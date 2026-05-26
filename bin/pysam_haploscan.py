#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
pysam_haploscan.py — Read-level haplotype scanning for ABO amplicons

Replaces the two-step  samtools mpileup → stats_from_pileup.py  pipeline with a
single pysam-based pass that works directly on aligned ONT reads.

Outputs
-------
{prefix}.AlignmentStatistics.tsv
    Identical column layout to stats_from_pileup.py — all downstream scripts
    (aggregate_abo_reports.py, predict_abo_phenotype.py) are unchanged.

{prefix}.Haplotypes.tsv
    Per-read haplotype table.  Each row records which diagnostic positions on
    that read carry a non-reference allele.  Because an ONT read can span an
    entire amplicon, two variants on the same row are CONFIRMED on the same
    molecule — enabling phasing without any wet-lab changes.

ABOReadPolymorphisms.txt
    Same polymorphic-position summary as before.

Phasing note
------------
c.297G (exon6 pos58) and c.1061del (exon7 pos685) lie on DIFFERENT amplicons
and therefore different BAM files.  Within-amplicon phasing (e.g. confirming
c.1032A + c.1061del on the same exon7 read) is fully supported here.
Cross-amplicon phasing requires a post-processing step that correlates exon6
and exon7 haplotype tables by read name (reads that span both amplicons).

Usage
-----
    pysam_haploscan.py -b sample.bam -f reference.fasta -o prefix
    pysam_haploscan.py -b sample.bam -f reference.fasta -o prefix -s ABOReadPolymorphisms.txt -v
"""

import argparse
import logging
import sys
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Dict, List, Optional

import pysam

__author__ = "Fredrick Mobegi"
__copyright__ = (
    "Copyright 2024, ABO blood group typing using third-generation sequencing (TGS) technology"
)
__credits__ = ["Fredrick Mobegi", "Benedict Matern", "Mathijs Groeneweg"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Fredrick Mobegi"
__email__ = "fredrick.mobegi@health.wa.gov.au"
__status__ = "Production"


# ---------------------------------------------------------------------------
# Exon classification
# ---------------------------------------------------------------------------

class ExonType(Enum):
    EXON6   = "exon6"
    EXON7   = "exon7"
    UNKNOWN = "unknown"


EXON6_LENGTH_RANGE = (130, 140)
EXON7_LENGTH_RANGE = (800, 830)


def determine_exon_type(ref_length: int) -> ExonType:
    if EXON6_LENGTH_RANGE[0] <= ref_length <= EXON6_LENGTH_RANGE[1]:
        return ExonType.EXON6
    if EXON7_LENGTH_RANGE[0] <= ref_length <= EXON7_LENGTH_RANGE[1]:
        return ExonType.EXON7
    return ExonType.UNKNOWN


# ---------------------------------------------------------------------------
# Diagnostic positions  (1-based, relative to each trimmed amplicon reference)
# ---------------------------------------------------------------------------

# Positions where an indel IS the diagnostic variant.
# At these positions indel counts are included in the frequency denominator.
# All other positions treat indels as noise and use ATGC-only denominator.
INDEL_DIAGNOSTIC: Dict[ExonType, frozenset] = {
    ExonType.EXON6: frozenset({22}),        # c.261delG  → O1 marker
    ExonType.EXON7: frozenset({431, 685}),  # 431 = O/A boundary; 685 = c.1061delC → A2.01
}

# Mirror of stats_from_pileup.LOW_COVERAGE_THRESHOLD
_LOW_COVERAGE_THRESHOLD = 200

# All 25 diagnostic positions interrogated for haplotype tagging (1-based).
# These match the positions used by aggregate_abo_reports.py / predict_abo_phenotype.py.
HAPLOTYPE_POSITIONS: Dict[ExonType, frozenset] = {
    ExonType.EXON6: frozenset({22, 27, 29, 58}),
    ExonType.EXON7: frozenset({
        93,                               # c.467 A1.02/A2 discriminator
        152, 153, 165, 283, 329, 348,     # A2 subtype panel
        368, 397, 404,                    # A2 subtype panel
        422, 428, 429, 431,               # Primary typing positions
        455, 533, 635, 658, 680,          # A2 subtype panel
        685,                              # c.1059 (A2 panel)
        687,                              # c.1061delC → A2.01 (KEY indel)
    }),
}

NUCLEOTIDES = frozenset({"A", "G", "C", "T"})

TSV_HEADER = "\t".join([
    "Ref_Position_1based", "Ref_Base", "Match_Percent", "Mismatch_Percent",
    "Insertion_Percent", "Deletion_Percent", "A_Percent", "G_Percent",
    "C_Percent", "T_Percent", "Depth",
])


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class PositionStats:
    pos:               int
    ref_base:          str
    depth:             int = 0
    match_percent:     int = 0
    mismatch_percent:  int = 0
    insertion_percent: int = 0
    deletion_percent:  int = 0
    A_percent:         int = 0
    G_percent:         int = 0
    C_percent:         int = 0
    T_percent:         int = 0


@dataclass
class ReadHaplotype:
    """Alleles observed at diagnostic positions on a single read."""
    read_name: str
    exon:      str
    alleles:   Dict[int, str] = field(default_factory=dict)

    def to_string(self) -> str:
        """Return e.g.  p22=del;p687=del  or  REF  for an all-reference read."""
        if not self.alleles:
            return "REF"
        return ";".join(f"p{pos}={allele}" for pos, allele in sorted(self.alleles.items()))


# ---------------------------------------------------------------------------
# Per-position frequency computation  (→ AlignmentStatistics.tsv)
# ---------------------------------------------------------------------------

def compute_position_stats(
    bam:              pysam.AlignmentFile,
    ref_name:         str,
    ref_seq:          str,
    exon_type:        ExonType,
    min_base_quality: int = 0,
    min_map_quality:  int = 0,
) -> List[PositionStats]:
    """
    Compute ATGC + indel frequencies at every reference position using the
    pysam pileup engine.  Logic mirrors stats_from_pileup.py exactly so that
    the TSV output is bit-for-bit compatible with the existing downstream code.
    """
    ref_length      = len(ref_seq)
    indel_positions = INDEL_DIAGNOSTIC.get(exon_type, frozenset())
    results: List[PositionStats] = []

    for col in bam.pileup(
        ref_name,
        0,
        ref_length,
        min_base_quality=min_base_quality,
        min_mapping_quality=min_map_quality,
        truncate=True,
        ignore_overlaps=False,
        stepper="nofilter",
    ):
        pos0     = col.reference_pos
        pos1     = pos0 + 1                                   # 1-based
        ref_base = ref_seq[pos0].upper() if pos0 < ref_length else "N"
        depth    = col.nsegments

        if depth == 0:
            results.append(PositionStats(pos=pos1, ref_base=ref_base))
            continue

        counts: Dict[str, int] = {"A": 0, "G": 0, "C": 0, "T": 0, "ins": 0, "del": 0}

        for pr in col.pileups:
            if pr.is_refskip:
                continue
            if pr.is_del:
                counts["del"] += 1
            else:
                qpos = pr.query_position
                if qpos is not None and pr.alignment.query_sequence:
                    base = pr.alignment.query_sequence[qpos].upper()
                    if base in NUCLEOTIDES:
                        counts[base] += 1
                    # An insertion FOLLOWS this position in this read
                    if pr.indel > 0:
                        counts["ins"] += 1

        # Mirror stats_from_pileup._should_include_indels():
        # Always include at explicit diagnostic positions; for all other positions
        # default to True unless coverage is low (same logic as original pipeline).
        if pos1 in indel_positions:
            include_indels = True
        elif depth < _LOW_COVERAGE_THRESHOLD:
            # Low coverage: exclude indels for exon6 or long-reference exon7
            # (mirrors stats_from_pileup: exon6 → False, total_rows > 140 → False)
            include_indels = not (exon_type == ExonType.EXON6 or ref_length > 140)
        else:
            include_indels = True  # High coverage: always include (original default)
        total_bases    = sum(counts[b] for b in NUCLEOTIDES)
        total_all      = total_bases + counts["ins"] + counts["del"]

        if include_indels:
            denom     = total_all if total_all > 0 else 1
            ins_pct   = int(counts["ins"] / denom * 100)
            del_pct   = int(counts["del"] / denom * 100)
            base_pcts = {b: int(counts[b] / denom * 100) for b in NUCLEOTIDES}
        else:
            denom     = total_bases if total_bases > 0 else 1
            ins_pct   = 0
            del_pct   = 0
            base_pcts = {b: int(counts[b] / denom * 100) for b in NUCLEOTIDES}
            # Normalise so ATGC sums exactly to 100 (same rounding fix as original)
            atgc_sum = sum(base_pcts.values())
            if atgc_sum != 100 and atgc_sum > 0:
                base_pcts[ref_base] = base_pcts.get(ref_base, 0) + (100 - atgc_sum)

        match_pct    = base_pcts.get(ref_base, 0)
        mismatch_pct = sum(v for k, v in base_pcts.items() if k != ref_base)

        results.append(PositionStats(
            pos=pos1,
            ref_base=ref_base,
            depth=depth,
            match_percent=match_pct,
            mismatch_percent=mismatch_pct,
            insertion_percent=ins_pct,
            deletion_percent=del_pct,
            A_percent=base_pcts.get("A", 0),
            G_percent=base_pcts.get("G", 0),
            C_percent=base_pcts.get("C", 0),
            T_percent=base_pcts.get("T", 0),
        ))

    return results


# ---------------------------------------------------------------------------
# Per-read haplotype computation  (→ Haplotypes.tsv)
# ---------------------------------------------------------------------------

def compute_read_haplotypes(
    bam:             pysam.AlignmentFile,
    ref_name:        str,
    ref_seq:         str,
    exon_type:       ExonType,
    min_map_quality: int = 0,
) -> List[ReadHaplotype]:
    """
    Iterate every primary aligned read and record which diagnostic positions
    carry a non-reference allele.

    For ONT amplicons the entire amplicon fits on a single read, so two
    variants in the same row are CONFIRMED co-occurring on the same molecule.

    Deletion at a diagnostic indel position  → allele recorded as "del"
    Insertion following a diagnostic position → allele recorded as "ins"
    SNP at any diagnostic SNP position        → allele recorded as the alt base
    Reference base                            → position omitted from alleles dict
    Position not covered by read              → position omitted from alleles dict
    """
    diag_positions  = HAPLOTYPE_POSITIONS.get(exon_type, frozenset())
    indel_positions = INDEL_DIAGNOSTIC.get(exon_type, frozenset())

    if not diag_positions:
        return []

    haplotypes: List[ReadHaplotype] = []

    for read in bam.fetch(ref_name):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < min_map_quality:
            continue
        if read.cigartuples is None or read.query_sequence is None:
            continue

        # Build ref_pos (0-based) → query_pos mapping from aligned pairs.
        # query_pos is None where the read has a deletion spanning that ref pos.
        ref_to_qpos: Dict[int, Optional[int]] = {}
        # Track ref positions followed by an insertion in this read
        ref_followed_by_ins: set = set()

        prev_rpos: Optional[int] = None
        for qpos, rpos in read.get_aligned_pairs(matches_only=False, with_seq=False):
            if rpos is not None:
                ref_to_qpos[rpos] = qpos   # qpos is None ↔ deletion
                prev_rpos = rpos
            elif qpos is not None and prev_rpos is not None:
                # Insertion in the read (no ref position)
                ref_followed_by_ins.add(prev_rpos)

        alleles: Dict[int, str] = {}

        for pos1 in diag_positions:
            rpos0 = pos1 - 1   # 0-based ref coordinate

            if rpos0 not in ref_to_qpos:
                continue    # read does not cover this position

            qpos = ref_to_qpos[rpos0]
            ref_base = ref_seq[rpos0].upper() if rpos0 < len(ref_seq) else "N"

            if pos1 in indel_positions:
                # Deletion spanning this position
                if qpos is None:
                    alleles[pos1] = "del"
                # Insertion following this position
                elif rpos0 in ref_followed_by_ins:
                    alleles[pos1] = "ins"
                # Reference base at indel position → no event (omit)
            else:
                # SNP position
                if qpos is None:
                    alleles[pos1] = "del"   # unexpected deletion
                else:
                    base = read.query_sequence[qpos].upper()
                    if base != ref_base and base in NUCLEOTIDES:
                        alleles[pos1] = base

        haplotypes.append(ReadHaplotype(
            read_name=read.query_name or "unknown",
            exon=exon_type.value,
            alleles=alleles,
        ))

    return haplotypes


# ---------------------------------------------------------------------------
# Co-occurrence / phasing analysis
# ---------------------------------------------------------------------------

def analyse_cooccurrence(
    haplotypes: List[ReadHaplotype],
    pos_a: int,
    allele_a: str,
    pos_b: int,
    allele_b: str,
    logger: logging.Logger,
) -> None:
    """
    Log how often two alleles co-occur on the same read.

    Example: pos685=A + pos687=del → confirms c.1032A and c.1061del are on
    the SAME molecule (same allele), ruling out trans configuration.

    Cross-amplicon phasing (e.g. exon6 pos58 + exon7 pos687) requires merging
    the two Haplotypes.tsv files by read name in a separate step.
    """
    reads_a  = sum(1 for h in haplotypes if h.alleles.get(pos_a) == allele_a)
    reads_b  = sum(1 for h in haplotypes if h.alleles.get(pos_b) == allele_b)
    reads_ab = sum(
        1 for h in haplotypes
        if h.alleles.get(pos_a) == allele_a and h.alleles.get(pos_b) == allele_b
    )
    total = len(haplotypes) or 1
    pct   = 100.0 * reads_ab / total

    logger.info(
        f"Phasing  p{pos_a}={allele_a} ∩ p{pos_b}={allele_b}: "
        f"{reads_a} reads carry p{pos_a}={allele_a}, "
        f"{reads_b} carry p{pos_b}={allele_b}, "
        f"{reads_ab} carry BOTH ({pct:.1f}% of all reads)"
    )

    if reads_a > 0:
        frac = reads_ab / reads_a
        if frac >= 0.8:
            logger.info(
                f"  → PHASED: {100*frac:.0f}% of p{pos_a}={allele_a} reads "
                f"also carry p{pos_b}={allele_b}  (same allele confirmed)"
            )
        elif frac <= 0.2:
            logger.info(
                f"  → TRANS: only {100*frac:.0f}% of p{pos_a}={allele_a} reads "
                f"carry p{pos_b}={allele_b}  (likely on different alleles)"
            )
        else:
            logger.warning(
                f"  → AMBIGUOUS: {100*frac:.0f}% co-occurrence — manual review recommended"
            )


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_stats_tsv(stats: List[PositionStats], path: Path) -> None:
    with open(path, "w") as fh:
        fh.write(TSV_HEADER + "\n")
        for s in stats:
            fh.write("\t".join([
                str(s.pos), s.ref_base,
                str(s.match_percent), str(s.mismatch_percent),
                str(s.insertion_percent), str(s.deletion_percent),
                str(s.A_percent), str(s.G_percent),
                str(s.C_percent), str(s.T_percent),
                str(s.depth),
            ]) + "\n")


def write_haplotypes_tsv(haplotypes: List[ReadHaplotype], path: Path) -> None:
    with open(path, "w") as fh:
        fh.write("Read_Name\tExon\tVariant_Positions\tHaplotype\n")
        for h in haplotypes:
            var_pos = (
                ",".join(str(p) for p in sorted(h.alleles.keys()))
                if h.alleles else "."
            )
            fh.write(f"{h.read_name}\t{h.exon}\t{var_pos}\t{h.to_string()}\n")


def write_summary(
    stats: List[PositionStats],
    path: Path,
    threshold: int = 10,
) -> None:
    """Write ABOReadPolymorphisms.txt in the same format as stats_from_pileup.py."""
    logger = logging.getLogger(__name__)
    polymorphic = 0
    with open(path, "w") as fh:
        for s in stats:
            if (
                s.mismatch_percent  >= threshold
                or s.insertion_percent >= threshold
                or s.deletion_percent  >= threshold
            ):
                fh.write(f"(1-based) Position:{s.pos}, Reference Base={s.ref_base}\n")
                fh.write(f"Aligned Read Count:{s.depth}\n")
                fh.write("Mat\tMis\tIns\tDel\tA\tG\tC\tT\n")
                fh.write(
                    f"{s.match_percent}\t{s.mismatch_percent}\t"
                    f"{s.insertion_percent}\t{s.deletion_percent}\t"
                    f"{s.A_percent}\t{s.G_percent}\t{s.C_percent}\t{s.T_percent}\n\n"
                )
                polymorphic += 1
    logger.info(f"Found {polymorphic} polymorphic positions (threshold={threshold}%)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def setup_logging(verbose: bool) -> logging.Logger:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    return logging.getLogger(__name__)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Read-level ABO haplotype scanner (replaces samtools mpileup + stats_from_pileup).\n"
            "Produces identical AlignmentStatistics.tsv output PLUS a per-read Haplotypes.tsv."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -b sample.bam -f reference.fasta -o sample_prefix
  %(prog)s -b sample.bam -f reference.fasta -o sample_prefix -s ABOReadPolymorphisms.txt -v
        """,
    )
    parser.add_argument("-b", "--bam",      required=True,
                        help="Sorted, indexed BAM file")
    parser.add_argument("-f", "--fasta",    required=True,
                        help="Reference FASTA (must have .fai index)")
    parser.add_argument("-o", "--output",   required=True,
                        help="Output prefix — produces <prefix>.AlignmentStatistics.tsv "
                             "and <prefix>.Haplotypes.tsv")
    parser.add_argument("-s", "--summary",  default="ABOReadPolymorphisms.txt",
                        help="Polymorphic positions summary (default: ABOReadPolymorphisms.txt)")
    parser.add_argument("-t", "--threshold", type=int, default=10,
                        help="Polymorphism threshold %% (default: 10)")
    parser.add_argument("-q", "--min-mapq",  type=int, default=0,
                        help="Minimum mapping quality (default: 0)")
    parser.add_argument("-Q", "--min-baseq", type=int, default=0,
                        help="Minimum base quality for pileup (default: 0)")
    parser.add_argument("-v", "--verbose",  action="store_true",
                        help="Verbose logging")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    args   = parser.parse_args()
    logger = setup_logging(args.verbose)

    bam_path      = Path(args.bam)
    fasta_path    = Path(args.fasta)
    prefix        = args.output
    summary_path  = Path(args.summary)
    stats_path    = Path(f"{prefix}.AlignmentStatistics.tsv")
    haplo_path    = Path(f"{prefix}.Haplotypes.tsv")

    # ── Load reference ──────────────────────────────────────────────────────
    try:
        with pysam.FastaFile(str(fasta_path)) as fa:
            ref_names = fa.references
            if not ref_names:
                logger.error("No sequences found in reference FASTA")
                return 1
            ref_name = ref_names[0]
            ref_seq  = fa.fetch(ref_name).upper()
    except Exception as exc:
        logger.error(f"Cannot load reference FASTA: {exc}")
        return 1

    ref_length = len(ref_seq)
    exon_type  = determine_exon_type(ref_length)
    logger.info(f"Reference: {ref_name}  length={ref_length}  exon={exon_type.value}")

    # ── Open BAM ────────────────────────────────────────────────────────────
    try:
        bam = pysam.AlignmentFile(str(bam_path), "rb")
    except Exception as exc:
        logger.error(f"Cannot open BAM: {exc}")
        return 1

    try:
        # ── 1. Per-position frequencies ─────────────────────────────────────
        logger.info("Computing per-position allele frequencies …")
        stats = compute_position_stats(
            bam, ref_name, ref_seq, exon_type,
            min_base_quality=args.min_baseq,
            min_map_quality=args.min_mapq,
        )
        write_stats_tsv(stats, stats_path)
        logger.info(f"✓  AlignmentStatistics  → {stats_path}")

        write_summary(stats, summary_path, threshold=args.threshold)
        logger.info(f"✓  Polymorphisms summary → {summary_path}")

        # ── 2. Per-read haplotypes ───────────────────────────────────────────
        logger.info("Computing per-read haplotypes …")
        haplotypes = compute_read_haplotypes(
            bam, ref_name, ref_seq, exon_type,
            min_map_quality=args.min_mapq,
        )
        write_haplotypes_tsv(haplotypes, haplo_path)
        logger.info(f"✓  Haplotypes            → {haplo_path}  ({len(haplotypes)} reads)")

        # ── 3. Within-amplicon phasing checks ────────────────────────────────
        if haplotypes and exon_type == ExonType.EXON7:
            # c.1032G>A (pos685) + c.1061del (pos687) — both on exon7
            analyse_cooccurrence(haplotypes, 685, "A",   687, "del", logger)
            # c.907A (pos422 area — note: check actual position) + c.1061del
            analyse_cooccurrence(haplotypes, 422, "A",   687, "del", logger)

        if haplotypes and exon_type == ExonType.EXON6:
            # c.261del (pos22) — O1 marker; confirm it is homozygous or het
            dels  = sum(1 for h in haplotypes if "del" in h.alleles.get(22, ""))
            total = len(haplotypes) or 1
            logger.info(
                f"Exon6 pos22 deletion: {dels}/{total} reads "
                f"({100*dels/total:.1f}%)  — "
                + ("homozygous O1" if dels/total > 0.8 else
                   "heterozygous O1/non-O1" if dels/total > 0.2 else
                   "non-O1")
            )

    finally:
        bam.close()

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\n! Interrupted by user")
        sys.exit(130)
    except Exception as exc:
        print(f"! CRITICAL: unhandled exception: {exc}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
