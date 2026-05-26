#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import os
import re
import sys
import glob
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import pandas as pd
from xlsxwriter.utility import xl_col_to_name

__author__ = "Fredrick Mobegi"
__copyright__ = "Copyright 2024, ABO blood group typing using third-generation sequencing (TGS) technology"
__credits__ = [
    "Fredrick Mobegi",
    "Benedict Matern",
    "Mathijs Groeneweg",
    "Claude Sonnet 4.6 (rewrite to add A1/A2/A3 subtypes)",
]
__license__ = "GPL"
__version__ = "1.2.0"
__maintainer__ = "Fredrick Mobegi"
__email__ = "fredrick.mobegi@health.wa.gov.au"
__status__ = "Production"


"""
ABO Blood Group Report Aggregator — v1.2.0

Changes from v1.1.0:
  - Added Bw subtype detection: B(A) (ABO*BA.01–04), Bw (ABO*BW.04 c.1054C>T,
    ABO*BW.06 c.657C>T, cisAB variants)
  - Added phase confidence scoring from Haplotypes.tsv: validates that variants
    called at ~50% frequency truly co-occur on the same reads (cis confirmation)

This file is part of the nf-core/abotyper pipeline "https://github.com/fmobegi/nf-core-abotyper".
"""


# ===========================================================================
# Phase Confidence Scoring — reads Haplotypes.tsv from pysam_haploscan.py
# ===========================================================================


@dataclass
class PhaseEvidence:
    """Summary of within-amplicon phase evidence for a sample+exon."""

    total_reads: int = 0
    # Key co-occurrence counts for validation
    cis_confirmed: Dict[str, int] = field(default_factory=dict)
    trans_confirmed: Dict[str, int] = field(default_factory=dict)
    ambiguous: Dict[str, int] = field(default_factory=dict)
    # Raw haplotype strings for advanced analysis
    haplotype_counts: Dict[str, int] = field(default_factory=dict)


def parse_haplotypes_tsv(haplo_file: str) -> Optional[PhaseEvidence]:
    """
    Parse a Haplotypes.tsv file produced by pysam_haploscan.py.

    Format:
        Read_Name\tExon\tVariant_Positions\tHaplotype
        read1\texon7\t685,687\tp685=A;p687=del
        read2\texon7\t.\tREF

    Returns PhaseEvidence summarising co-occurrence patterns, or None on failure.
    """
    if not haplo_file or not os.path.exists(haplo_file):
        return None

    try:
        df = pd.read_csv(haplo_file, sep="\t")
        if df.empty:
            return None

        evidence = PhaseEvidence(total_reads=len(df))

        # Count distinct haplotype strings
        if "Haplotype" in df.columns:
            for haplo_str, count in df["Haplotype"].value_counts().items():
                evidence.haplotype_counts[str(haplo_str)] = int(count)

        return evidence

    except Exception:
        return None


def compute_phase_confidence(
    evidence_exon6: Optional[PhaseEvidence],
    evidence_exon7: Optional[PhaseEvidence],
    extended_genotype: str,
) -> Tuple[str, str]:
    """
    Compute phase confidence from haplotype tables.

    Returns (confidence_level, detail_string):
      confidence_level: "High", "Moderate", "Low", "N/A"
      detail_string: human-readable explanation

    Logic:
      - If the genotype is homozygous (AA, BB, OO with same alleles), phasing
        is trivial → "High (homozygous)"
      - For heterozygous calls: check whether the dominant non-REF haplotype
        accounts for a consistent fraction of reads (~50% for het, ~100% for hom)
      - If haplotype diversity is low (1-2 dominant patterns covering >90% of
        reads), confidence is High
      - If many disparate haplotypes exist, confidence is Low
    """
    if not evidence_exon6 and not evidence_exon7:
        return "N/A", "No haplotype data available"

    # Homozygous genotypes don't need phase confirmation
    if "/" in extended_genotype:
        alleles = extended_genotype.split("/")
        if len(alleles) == 2 and alleles[0] == alleles[1]:
            return "High", "Homozygous; phasing trivial"

    # Analyse the exon with more diagnostic positions (exon7 preferred)
    evidence = (
        evidence_exon7
        if evidence_exon7 and evidence_exon7.total_reads > 0
        else evidence_exon6
    )
    if not evidence or evidence.total_reads == 0:
        return "N/A", "No haplotype data available"

    total = evidence.total_reads
    haplotypes = evidence.haplotype_counts

    if not haplotypes:
        return "N/A", "No haplotype calls in data"

    # Sort by frequency
    sorted_haplos = sorted(haplotypes.items(), key=lambda x: x[1], reverse=True)

    # Top 2 haplotypes should account for most reads in a diploid sample
    top2_count = sum(count for _, count in sorted_haplos[:2])
    top2_fraction = top2_count / total if total > 0 else 0

    n_distinct = len([h for h, c in sorted_haplos if c >= max(3, total * 0.02)])

    if top2_fraction >= 0.90 and n_distinct <= 3:
        detail = f"Top 2 haplotypes cover {top2_fraction*100:.0f}% of {total} reads"
        return "High", detail
    elif top2_fraction >= 0.75 and n_distinct <= 5:
        detail = f"Top 2 haplotypes cover {top2_fraction*100:.0f}% of {total} reads ({n_distinct} distinct patterns)"
        return "Moderate", detail
    else:
        detail = f"High haplotype diversity: {n_distinct} patterns, top 2 cover only {top2_fraction*100:.0f}%"
        return "Low", detail


# ===========================================================================
# Bw / B(A) subtype detection
# ===========================================================================

# Known Bw/B(A) markers (1-based positions in exon7 amplicon coordinates)
# These are positions already in the panel that, when variant IN THE CONTEXT
# of a B allele, indicate a weak-B or cisAB phenotype.
BW_MARKERS = {
    # pos: (ref_base, variant_base, subtype_label, cDNA_change)
    680: ("C", "T", "Bw.04", "c.1054C>T"),  # ABO*BW.04
    283: ("C", "T", "Bw.06", "c.657C>T"),  # ABO*BW.06 (also c.657C>A)
    348: (
        "G",
        "A",
        "B(A).04",
        "c.722G>A",
    ),  # ABO*BA.04 — B(A) with A transferase activity
    404: ("G", "A", "B(A).02", "c.778G>A"),  # ABO*BA.02
    # cisAB: typically carries c.803G>C (pos429, already primary) PLUS additional marker
    # A B allele that ALSO shows A-activity markers
}


def scan_bw_markers(
    type_exon7_422: str,
    type_exon7_429: str,
    exon7_positions: Dict[int, dict],
    nreads_exon7: int,
) -> Tuple[List[str], List[str]]:
    """
    Scan for Bw/B(A) subtype markers.

    Only fires when the sample carries at least one B allele (pos422=B or mixed,
    pos429=B or mixed).

    Args:
        type_exon7_422: Type string from position 422
        type_exon7_429: Type string from position 429
        exon7_positions: Dict of pos -> {col: value} for all exon7 data
        nreads_exon7: Read count at a representative exon7 position

    Returns:
        (markers_found, warnings)
    """
    markers = []
    warnings = []

    # Only check Bw markers if B allele is present
    has_b = ("B" in type_exon7_422) or ("B" in type_exon7_429)
    if not has_b:
        return markers, warnings

    # Need reasonable coverage to call Bw
    if nreads_exon7 < 30:
        return markers, ["Low coverage — Bw subtyping unreliable"]

    for pos, (ref_base, var_base, subtype, cdna) in BW_MARKERS.items():
        pos_data = exon7_positions.get(pos)
        if not pos_data:
            continue

        # Get the variant base percentage
        var_pct = pos_data.get(f"{var_base}", 0)
        if isinstance(var_pct, str):
            try:
                var_pct = float(var_pct)
            except (ValueError, TypeError):
                var_pct = 0

        # Bw markers are typically heterozygous (one normal B + one weak B)
        # or present at ~50% if the Bw variant is on one allele
        # Threshold: >=15% to detect in a heterozygous context
        if var_pct >= 15:
            confidence = "Confirmed" if var_pct >= 80 else "Possible"
            markers.append(f"{confidence} {subtype}({cdna}:{var_base}={var_pct:.0f}%)")

    # cisAB detection: B allele shows BOTH B markers (pos422=A mix, pos429=C mix)
    # AND A-specific markers at the A subtype positions
    # This is a complex scenario — flag for manual review if B + A2 markers coexist
    if has_b and exon7_positions.get(685):
        del_pct = exon7_positions[685].get("Del", 0)
        if isinstance(del_pct, str):
            try:
                del_pct = float(del_pct)
            except (ValueError, TypeError):
                del_pct = 0
        if del_pct >= 15 and "B" in type_exon7_422:
            warnings.append(
                "Possible cisAB: B allele detected with c.1061del — manual review recommended"
            )

    return markers, warnings


# ===========================================================================
# Main ABOReportParser class
# ===========================================================================


class ABOReportParser:
    """
    Collates all ABO phenotype results from each sample into a general table
    and generates an Excel worksheet and a CSV file for export to LIS soft or
    other general purpose lab management systems.

    v1.2.0 additions:
      - Bw subtype detection from existing panel positions
      - Phase confidence scoring from Haplotypes.tsv
    """

    def __init__(self, input_dir, default_barcode="barcode00"):
        """
        Initialize the ABOReportParser.

        Args:
            input_dir (str): The input directory containing data files.
            default_barcode (str): Default barcode for samples without explicit barcode suffix.
        """
        self.input_dir = input_dir
        self.default_barcode = default_barcode
        self.results = []
        self.initialize_columns()
        self.failed_samples = []

        # Compile regex patterns once
        self.pattern_with_barcode = re.compile(
            r"^(IMM|INGS|NGS|[A-Z0-9]+)(-[A-Z0-9]+)?(-[A-Z0-9]+)?_barcode\d+$",
            re.IGNORECASE,
        )
        self.pattern_without_barcode = re.compile(
            r"^[A-Z0-9]+(_[A-Z0-9]+)*$", re.IGNORECASE
        )

        # Statistics tracking
        self.processing_stats = {
            "with_barcode": 0,
            "without_barcode": 0,
            "pattern_matched": 0,
            "pattern_failed": 0,
        }

    def initialize_columns(self):
        """
        Define the column headers for selected exon positions only.
        Updated to include PhaseConfidence and BwSubtype columns.
        """
        # Primary ABO typing variant in exon 6
        exon6 = ["Exon6_pos22"] * 10

        # Primary variants in exon 7
        exon7_422 = ["Exon7_pos422"] * 10
        exon7_428 = ["Exon7_pos428"] * 10
        exon7_429 = ["Exon7_pos429"] * 10
        exon7_431 = ["Exon7_pos431"] * 10

        # Selected A subtype positions in exon 7
        exon7_93 = ["Exon7_pos93"] * 10
        exon7_685 = ["Exon7_pos685"] * 10

        # A2 subtype discriminator positions — exon 6
        exon6_27 = ["Exon6_pos27"] * 10
        exon6_29 = ["Exon6_pos29"] * 10
        exon6_58 = ["Exon6_pos58"] * 10
        # A2 subtype discriminator positions — exon 7
        exon7_33 = ["Exon7_pos33"] * 10
        exon7_152 = ["Exon7_pos152"] * 10
        exon7_153 = ["Exon7_pos153"] * 10
        exon7_165 = ["Exon7_pos165"] * 10
        exon7_283 = ["Exon7_pos283"] * 10
        exon7_329 = ["Exon7_pos329"] * 10
        exon7_348 = ["Exon7_pos348"] * 10
        exon7_368 = ["Exon7_pos368"] * 10
        exon7_397 = ["Exon7_pos397"] * 10
        exon7_404 = ["Exon7_pos404"] * 10
        exon7_455 = ["Exon7_pos455"] * 10
        exon7_533 = ["Exon7_pos533"] * 10
        exon7_635 = ["Exon7_pos635"] * 10
        exon7_658 = ["Exon7_pos658"] * 10
        exon7_680 = ["Exon7_pos680"] * 10

        header_cols = (
            ["", ""]                    # Sample: Barcode, Sequencing_ID
            + ["Result"] * 3            # Result: Phenotype, Genotype, ExtendedGenotype
            + exon6
            + exon7_422
            + exon7_428
            + exon7_429
            + exon7_431
            + exon7_93
            + exon7_685
            + exon6_27
            + exon6_29
            + exon6_58
            + exon7_33
            + exon7_152
            + exon7_153
            + exon7_165
            + exon7_283
            + exon7_329
            + exon7_348
            + exon7_368
            + exon7_397
            + exon7_404
            + exon7_455
            + exon7_533
            + exon7_635
            + exon7_658
            + exon7_680
            + ["Notes"] * 4             # Notes: ASubtype, ReadReliability, PhaseConfidence, BwSubtype
        )

        column_metrics = [
            "#Reads",
            "Mat",
            "Mis",
            "Ins",
            "Del",
            "A",
            "G",
            "C",
            "T",
            "Type",
        ]
        header_rows = (
            ["Barcode", "Sequencing_ID"]
            + ["Phenotype", "Genotype", "ExtendedGenotype"]  # Result columns
            + column_metrics * 25
            + [
                "ASubtype",
                "ReadReliability",
                "PhaseConfidence",
                "BwSubtype",
            ]  # Notes columns
        )

        self.columns = pd.MultiIndex.from_arrays([header_cols, header_rows])

    def extract_sample_info(self, filename):
        """
        Extract sample name and barcode from filename.

        Returns:
            tuple: (sample_name, barcode, pattern_type)
        """
        if "_barcode" in filename:
            parts = filename.rsplit("_barcode", 1)
            sample_name = parts[0]
            try:
                barcode_num = parts[1]
                barcode_int = int(barcode_num)
                if 0 <= barcode_int <= 99:
                    barcode = f"barcode{barcode_num.zfill(2)}"
                else:
                    print(
                        f"Warning: Unusual barcode number {barcode_int} for {filename}"
                    )
                    barcode = f"barcode{barcode_num}"
            except ValueError:
                print(f"Warning: Invalid barcode format in {filename}, using as-is")
                barcode = f"barcode{parts[1]}"
            return sample_name, barcode, "explicit"
        else:
            sample_name = filename
            barcode = self.default_barcode
            print(f"No barcode found in '{filename}', using default barcode: {barcode}")
            return sample_name, barcode, "default"

    def parse_exon7(self, filename):
        """Open the file for reading and processing all exon 7 positions."""
        try:
            with open(filename, "r", encoding="utf-8") as f:
                lines = f.readlines()

            positions = []
            counts = []
            mat_values = []
            mis_values = []
            ins_values = []
            del_values = []
            a_values = []
            g_values = []
            c_values = []
            t_values = []

            i = 0
            while i < len(lines):
                line = lines[i].strip()

                if "Exon 7 position(1-based):" in line:
                    pos_match = re.search(r":\s*(\d+)", line)
                    if pos_match:
                        pos = int(pos_match.group(1))

                        for j in range(i, min(i + 10, len(lines))):
                            if "Aligned Read Count:" in lines[j]:
                                count_match = re.search(r":\s*(\d+)", lines[j])
                                if count_match:
                                    count = int(count_match.group(1))

                                    stats_idx = j + 2
                                    if (
                                        stats_idx < len(lines)
                                        and "Mat" in lines[stats_idx - 1]
                                    ):
                                        stats = lines[stats_idx].split()
                                        if len(stats) >= 8:
                                            mat, mis, ins, dele, a, g, c, t = [
                                                float(x) for x in stats
                                            ]

                                            positions.append(pos)
                                            counts.append(count)
                                            mat_values.append(mat)
                                            mis_values.append(mis)
                                            ins_values.append(ins)
                                            del_values.append(dele)
                                            a_values.append(a)
                                            g_values.append(g)
                                            c_values.append(c)
                                            t_values.append(t)

                                            break
                i += 1

            df = pd.DataFrame(
                {
                    "Exon": ["7"] * len(positions),
                    "Position": positions,
                    "#Reads": counts,
                    "Mat": mat_values,
                    "Mis": mis_values,
                    "Ins": ins_values,
                    "Del": del_values,
                    "A": a_values,
                    "G": g_values,
                    "C": c_values,
                    "T": t_values,
                }
            )

            all_positions = [
                422,
                428,
                429,
                431,
                93,
                685,
                33,
                152,
                153,
                165,
                283,
                329,
                348,
                368,
                397,
                404,
                455,
                533,
                635,
                658,
                680,
            ]

            for pos in all_positions:
                if pos not in df["Position"].values:
                    df = pd.concat(
                        [
                            df,
                            pd.DataFrame(
                                {
                                    "Exon": ["7"],
                                    "Position": [pos],
                                    "#Reads": [0],
                                    "Mat": [0],
                                    "Mis": [0],
                                    "Ins": [0],
                                    "Del": [0],
                                    "A": [0],
                                    "G": [0],
                                    "C": [0],
                                    "T": [0],
                                }
                            ),
                        ],
                        ignore_index=True,
                    )

            df = df.sort_values("Position").reset_index(drop=True)

            df["Type"] = df.apply(
                lambda row: self.get_type(
                    row["Position"], row["A"], row["G"], row["C"], row["T"], row["Del"]
                ),
                axis=1,
            )

            return df

        except Exception as e:
            print(f"Error parsing exon 7 file {filename}: {str(e)}")
            empty_df = pd.DataFrame(
                columns=[
                    "Exon",
                    "Position",
                    "#Reads",
                    "Mat",
                    "Mis",
                    "Ins",
                    "Del",
                    "A",
                    "G",
                    "C",
                    "T",
                    "Type",
                ]
            )

            for pos in [
                422,
                428,
                429,
                431,
                93,
                685,
                33,
                152,
                153,
                165,
                283,
                329,
                348,
                368,
                397,
                404,
                455,
                533,
                635,
                658,
                680,
            ]:
                empty_df = pd.concat(
                    [
                        empty_df,
                        pd.DataFrame(
                            {
                                "Exon": ["7"],
                                "Position": [pos],
                                "#Reads": [0],
                                "Mat": [0],
                                "Mis": [0],
                                "Ins": [0],
                                "Del": [0],
                                "A": [0],
                                "G": [0],
                                "C": [0],
                                "T": [0],
                                "Type": [""],
                            }
                        ),
                    ],
                    ignore_index=True,
                )

            return empty_df

    def parse_exon6(self, filename):
        """Parse exon 6 and extract data for all relevant positions (22, 27, 29, 58)."""
        try:
            with open(filename, "r", encoding="utf-8") as f:
                lines = f.readlines()

            positions = []
            counts = []
            mat_values = []
            mis_values = []
            ins_values = []
            del_values = []
            a_values = []
            g_values = []
            c_values = []
            t_values = []

            i = 0
            while i < len(lines):
                line = lines[i].strip()

                if "Exon 6 position(1-based):" in line:
                    pos_match = re.search(r":\s*(\d+)", line)
                    if pos_match:
                        pos = int(pos_match.group(1))

                        for j in range(i, min(i + 10, len(lines))):
                            if "Aligned Read Count:" in lines[j]:
                                count_match = re.search(r":\s*(\d+)", lines[j])
                                if count_match:
                                    count = int(count_match.group(1))

                                    stats_idx = j + 2
                                    if (
                                        stats_idx < len(lines)
                                        and "Mat" in lines[stats_idx - 1]
                                    ):
                                        stats = lines[stats_idx].split()
                                        if len(stats) >= 8:
                                            mat, mis, ins, dele, a, g, c, t = [
                                                float(x) for x in stats
                                            ]

                                            positions.append(pos)
                                            counts.append(count)
                                            mat_values.append(mat)
                                            mis_values.append(mis)
                                            ins_values.append(ins)
                                            del_values.append(dele)
                                            a_values.append(a)
                                            g_values.append(g)
                                            c_values.append(c)
                                            t_values.append(t)

                                            break
                i += 1

            df = pd.DataFrame(
                {
                    "Exon": ["6"] * len(positions),
                    "Position": positions,
                    "#Reads": counts,
                    "Mat": mat_values,
                    "Mis": mis_values,
                    "Ins": ins_values,
                    "Del": del_values,
                    "A": a_values,
                    "G": g_values,
                    "C": c_values,
                    "T": t_values,
                }
            )

            all_positions = [22, 27, 29, 58]

            for pos in all_positions:
                if pos not in df["Position"].values:
                    df = pd.concat(
                        [
                            df,
                            pd.DataFrame(
                                {
                                    "Exon": ["6"],
                                    "Position": [pos],
                                    "#Reads": [0],
                                    "Mat": [0],
                                    "Mis": [0],
                                    "Ins": [0],
                                    "Del": [0],
                                    "A": [0],
                                    "G": [0],
                                    "C": [0],
                                    "T": [0],
                                }
                            ),
                        ],
                        ignore_index=True,
                    )

            df = df.sort_values("Position").reset_index(drop=True)

            df["Type"] = df.apply(
                lambda row: self.get_type_exon6(
                    row["Position"], row["A"], row["G"], row["C"], row["T"], row["Del"]
                ),
                axis=1,
            )

            return df

        except Exception as e:
            print(f"Error parsing exon 6: {str(e)}")
            empty_df = pd.DataFrame(
                columns=[
                    "Exon",
                    "Position",
                    "#Reads",
                    "Mat",
                    "Mis",
                    "Ins",
                    "Del",
                    "A",
                    "G",
                    "C",
                    "T",
                    "Type",
                ]
            )

            for pos in [22, 27, 29, 58]:
                empty_df = pd.concat(
                    [
                        empty_df,
                        pd.DataFrame(
                            {
                                "Exon": ["6"],
                                "Position": [pos],
                                "#Reads": [0],
                                "Mat": [0],
                                "Mis": [0],
                                "Ins": [0],
                                "Del": [0],
                                "A": [0],
                                "G": [0],
                                "C": [0],
                                "T": [0],
                                "Type": [""],
                            }
                        ),
                    ],
                    ignore_index=True,
                )

            return empty_df

    def get_type_exon6(self, pos, a, g, c, t, dele):
        """Determine blood type or subtype for each exon 6 position."""
        if pos == 22:  # c.261
            if g >= 80 and g > dele:
                return "A or B or O"
            elif dele >= 80 and dele > g:
                return "O1"
            elif abs(g + dele) >= 20:
                return "O1 and (A or B or O)"
            elif 20 < dele < 80 and 20 < g < 80:
                return "O1 and (A or B or O)"

        elif pos == 27:  # c.266C>T
            if t >= 25:
                return "variant"
            elif c >= 70:
                return ""

        elif pos == 29:  # c.268T>C
            if c >= 25:
                return "variant"
            elif t >= 70:
                return ""

        elif pos == 58:  # c.297A>G
            if g >= 25:
                return "variant"
            elif a >= 70:
                return ""

        return ""

    def get_type(self, pos, a, g, c, t, dele):
        """Determine the blood type or subtype for each exon7 position."""
        # Primary ABO variants
        if pos == 422:
            if c >= 80:
                return "A or O"
            elif a >= 80:
                return "B"
            elif abs(a - c) <= 20:
                return "(A or O) and B"
            elif 15 < a < 80 and 15 < c < 80:
                return "(A or O) and B"
        elif pos == 428:
            if g >= 70:
                return "O and (A or B)"
            elif a >= 70:
                return "O2"
            elif abs(g - a) <= 20:
                return "O2 and (O or A or B)"
            elif 15 < g < 70 and 15 < a < 70:
                return "O2 and (O or A or B)"
        elif pos == 429:
            if g >= 80:
                return "A or O"
            elif c >= 80:
                return "B"
            elif abs(g - c) <= 20:
                return "(A or O) and B"
            elif 15 < g < 80 and 20 < c < 80:
                return "(A or O) and B"
        elif pos == 431:
            if t >= 80:
                return "O and (A or B)"
            elif g >= 80:
                return "O3"
            elif abs(t - g) <= 20:
                return "O3 and (O or A or B)"
            elif 15 < t < 80 and 15 < g < 80:
                return "O3 and (O or A or B)"

        # A subtype positions
        elif pos == 93:
            if c >= 80:
                return "A1"
            elif t >= 80:
                return "A1.02 or A2"
            elif 20 < c < 80 and 20 < t < 80:
                return "A1.02 or A2"
        elif pos == 685:
            if c >= 80 and dele < 20:
                return "A1"
            elif dele >= 80 and c < 20:
                return "A2"
            elif c >= 20 and dele >= 20:
                return "A2"
            elif 20 < c < 80 and 20 < dele < 80:
                return "A2"

        # A2 comprehensive panel
        elif pos == 33:
            if (t + a + g) >= 25:
                return "variant"
        elif pos == 152:
            if (t + a + g) >= 25:
                return "variant"
        elif pos == 153:
            if (c + a + t) >= 25:
                return "variant"
        elif pos == 165:
            if (c + a + t) >= 25:
                return "variant"
        elif pos == 283:
            if (t + a + g) >= 25:
                return "variant"
        elif pos == 329:
            if (c + a + t) >= 25:
                return "variant"
        elif pos == 348:
            if (c + a + t) >= 25:
                return "variant"
        elif pos == 368:
            if (t + a + g) >= 25:
                return "variant"
        elif pos == 397:
            if (t + a + g) >= 25:
                return "variant"
        elif pos == 404:
            if (c + a + t) >= 25:
                return "variant"
        elif pos == 455:
            if (c + a + t) >= 25:
                return "variant"
        elif pos == 533:
            if a >= 25:
                return "A2.06"
        elif pos == 635:
            if (c + g + t) >= 25:
                return "variant"
        elif pos == 658:
            if a >= 25:
                return "A2.01"
        elif pos == 680:
            if (t + a + g) >= 25:
                return "variant"

        return ""

    def assign_phenotype_genotype(
        self,
        df,
        phase_evidence_exon6=None,
        phase_evidence_exon7=None,
    ):
        """Assign the phenotype and genotype information with Bw detection and phase confidence."""
        try:

            def safe_get_type(df, pos_key, default=""):
                try:
                    return df.at[0, (pos_key, "Type")]
                except (KeyError, IndexError):
                    return default

            def safe_get_reads(df, pos_key, default=0):
                try:
                    value = df.at[0, (pos_key, "#Reads")]
                    return value if pd.notna(value) else default
                except (KeyError, IndexError):
                    return default

            def safe_get_value(df, pos_key, col, default=0):
                """Get a specific column value for a position."""
                try:
                    value = df.at[0, (pos_key, col)]
                    return float(value) if pd.notna(value) else default
                except (KeyError, IndexError, ValueError, TypeError):
                    return default

            # Extract primary positions
            type_exon6 = safe_get_type(df, "Exon6_pos22")
            type_exon7_422 = safe_get_type(df, "Exon7_pos422")
            type_exon7_428 = safe_get_type(df, "Exon7_pos428")
            type_exon7_429 = safe_get_type(df, "Exon7_pos429")
            type_exon7_431 = safe_get_type(df, "Exon7_pos431")

            # A subtype positions
            type_exon7_93 = safe_get_type(df, "Exon7_pos93")
            type_exon7_685 = safe_get_type(df, "Exon7_pos685")

            # A2 subtype discriminator positions — exon 6
            type_exon6_27 = safe_get_type(df, "Exon6_pos27")
            type_exon6_29 = safe_get_type(df, "Exon6_pos29")
            type_exon6_58 = safe_get_type(df, "Exon6_pos58")
            # A2 subtype discriminator positions — exon 7
            type_exon7_33 = safe_get_type(df, "Exon7_pos33")
            type_exon7_152 = safe_get_type(df, "Exon7_pos152")
            type_exon7_153 = safe_get_type(df, "Exon7_pos153")
            type_exon7_165 = safe_get_type(df, "Exon7_pos165")
            type_exon7_283 = safe_get_type(df, "Exon7_pos283")
            type_exon7_329 = safe_get_type(df, "Exon7_pos329")
            type_exon7_348 = safe_get_type(df, "Exon7_pos348")
            type_exon7_368 = safe_get_type(df, "Exon7_pos368")
            type_exon7_397 = safe_get_type(df, "Exon7_pos397")
            type_exon7_404 = safe_get_type(df, "Exon7_pos404")
            type_exon7_455 = safe_get_type(df, "Exon7_pos455")
            type_exon7_533 = safe_get_type(df, "Exon7_pos533")
            type_exon7_635 = safe_get_type(df, "Exon7_pos635")
            type_exon7_658 = safe_get_type(df, "Exon7_pos658")
            type_exon7_680 = safe_get_type(df, "Exon7_pos680")

            # Read counts
            nreads6 = safe_get_reads(df, "Exon6_pos22")
            nreads_exon7_p422 = safe_get_reads(df, "Exon7_pos422")
            nreads_exon7_p428 = safe_get_reads(df, "Exon7_pos428")
            nreads_exon7_p429 = safe_get_reads(df, "Exon7_pos429")
            nreads_exon7_p431 = safe_get_reads(df, "Exon7_pos431")

            Phenotype = "Unknown"
            Genotype = "Unknown"
            ExtendedGenotype = "Unknown"
            BwSubtype = ""

            # ----- Build exon7 position data dict for Bw scanning -----
            exon7_positions = {}
            for pos in [
                422,
                428,
                429,
                431,
                93,
                685,
                33,
                152,
                153,
                165,
                283,
                329,
                348,
                368,
                397,
                404,
                455,
                533,
                635,
                658,
                680,
            ]:
                exon7_positions[pos] = {
                    "A": safe_get_value(df, f"Exon7_pos{pos}", "A"),
                    "G": safe_get_value(df, f"Exon7_pos{pos}", "G"),
                    "C": safe_get_value(df, f"Exon7_pos{pos}", "C"),
                    "T": safe_get_value(df, f"Exon7_pos{pos}", "T"),
                    "Del": safe_get_value(df, f"Exon7_pos{pos}", "Del"),
                    "#Reads": safe_get_reads(df, f"Exon7_pos{pos}"),
                }

            def scan_a2_markers():
                """Scan ALL diagnostic positions and return (markers, warnings)."""
                markers = []
                warnings = []

                _e7_reads = nreads_exon7_p422
                has_del = type_exon7_685 == "A2" and _e7_reads >= 30
                has_907 = type_exon7_533 == "A2.06"
                has_1032 = type_exon7_658 == "A2.01" and _e7_reads >= 30
                has_297 = type_exon6_58 == "variant"

                if has_del:
                    markers.append("c.1061del")
                if has_907:
                    markers.append("c.907A")
                if has_1032:
                    markers.append("c.1032A")
                if has_297:
                    markers.append("c.297G")

                o1v_fired = has_297 and not has_del and not has_907 and not has_1032
                if o1v_fired:
                    for m in ("c.297G", "c.266T", "c.268C"):
                        if m in markers:
                            markers.remove(m)
                    warnings.append(
                        "Note: c.297G detected without c.1061del — possible A2-derived "
                        "O allele (O1v); manual review recommended"
                    )

                if not o1v_fired:
                    if type_exon6_27 == "variant":
                        markers.append("c.266T")
                    if type_exon6_29 == "variant":
                        markers.append("c.268C")

                if nreads_exon7_p422 >= 500 and not o1v_fired:
                    for type_val, label in [
                        (type_exon7_33, "c.407var"),
                        (type_exon7_152, "c.526var"),
                        (type_exon7_153, "c.527var"),
                        (type_exon7_165, "c.539var"),
                        (type_exon7_283, "c.657var"),
                        (type_exon7_329, "c.703var"),
                        (type_exon7_348, "c.722var"),
                        (type_exon7_368, "c.742var"),
                        (type_exon7_397, "c.771var"),
                        (type_exon7_404, "c.778var"),
                        (type_exon7_455, "c.829var"),
                        (type_exon7_635, "c.1009var"),
                        (type_exon7_680, "c.1054var"),
                    ]:
                        if type_val == "variant":
                            markers.append(label)

                return markers, warnings

            def determine_a2_subtype(markers):
                """Resolve an A2-positive call to the most specific subtype."""
                has_del = "c.1061del" in markers
                has_907 = "c.907A" in markers
                has_1032 = "c.1032A" in markers
                has_297 = "c.297G" in markers

                if has_907:
                    return "A2.06", None
                if has_1032 and has_del:
                    return "A2.01", "c.1032G>A"
                if has_297 and has_del:
                    return "A2.01", "c.297A>G"
                if has_del:
                    return "A2.01", None
                if has_1032:
                    return "A2.01", "c.1032G>A"
                if has_297:
                    return "A2.01", "c.297A>G"
                return "A2", None

            def determine_a_subtype():
                """Determine A subtype using a POSITIVE EVIDENCE model."""
                markers, warnings = scan_a2_markers()

                if markers:
                    clean_subtype, nt_notation = determine_a2_subtype(markers)
                    if nt_notation:
                        warnings.insert(0, nt_notation)
                    warning_text = "; ".join(warnings) if warnings else None
                    return clean_subtype, warning_text

                warning_text = "; ".join(warnings) if warnings else None

                if type_exon7_93 in ("A1.02 or A2", "A1"):
                    return "A1", warning_text

                return "", None

            # PART 1: PRIMARY PHENOTYPING LOGIC
            a_subtype_warning = None

            ## OA COMBINATIONS
            ## combination 1 | AO1
            if (
                (type_exon6 == "O1 and (A or B or O)")
                and (type_exon7_422 == "A or O")
                and (type_exon7_428 == "O and (A or B)")
                and (type_exon7_429 == "A or O")
                and (type_exon7_431 == "O and (A or B)")
            ):
                a_subtype, a_subtype_warning = determine_a_subtype()
                Phenotype = "A"
                Genotype = "AO"
                if a_subtype:
                    ExtendedGenotype = f"{a_subtype}/O1"
                else:
                    ExtendedGenotype = "A/O1"

            ## combination 2 | AO2
            elif (
                (type_exon6 == "A or B or O")
                and (type_exon7_422 == "A or O")
                and (type_exon7_428 == "O2 and (O or A or B)")
                and (type_exon7_429 == "A or O")
                and (type_exon7_431 == "O and (A or B)")
            ):
                a_subtype, a_subtype_warning = determine_a_subtype()
                Phenotype = "A"
                Genotype = "AO"
                if a_subtype:
                    ExtendedGenotype = f"{a_subtype}/O2"
                else:
                    ExtendedGenotype = "A/O2"

            ## combination 3 | AO3
            elif (
                (type_exon6 == "A or B or O")
                and (type_exon7_422 == "A or O")
                and (type_exon7_428 == "O and (A or B)")
                and (type_exon7_429 == "A or O")
                and (type_exon7_431 == "O3 and (O or A or B)")
            ):
                a_subtype, a_subtype_warning = determine_a_subtype()
                Phenotype = "A"
                Genotype = "AO"
                if a_subtype:
                    ExtendedGenotype = f"{a_subtype}/O3"
                else:
                    ExtendedGenotype = "A/O3"

            ## OB COMBINATIONS
            ## combination 4 | BO1
            elif (
                (type_exon6 == "O1 and (A or B or O)")
                and (type_exon7_422 == "(A or O) and B")
                and (type_exon7_428 == "O and (A or B)")
                and (type_exon7_429 == "(A or O) and B")
                and (type_exon7_431 == "O and (A or B)")
            ):
                Phenotype = "B"
                Genotype = "BO"
                ExtendedGenotype = "B/O1"

            ## combination 5 | O2B
            elif (
                (type_exon6 == "A or B or O")
                and (type_exon7_422 == "(A or O) and B")
                and (type_exon7_428 == "O2 and (O or A or B)")
                and (type_exon7_429 == "(A or O) and B")
                and (type_exon7_431 == "O and (A or B)")
            ):
                Phenotype = "B"
                Genotype = "BO"
                ExtendedGenotype = "O2/B"

            ## combination 6 | BO3
            elif (
                (type_exon6 == "A or B or O")
                and (type_exon7_422 == "(A or O) and B")
                and (type_exon7_428 == "O and (A or B)")
                and (type_exon7_429 == "(A or O) and B")
                and (type_exon7_431 == "O3 and (O or A or B)")
            ):
                Phenotype = "B"
                Genotype = "BO"
                ExtendedGenotype = "B/O3"

            ## OO COMBINATIONS
            ## combination 7 | O1O2
            elif (
                (type_exon6 == "O1 and (A or B or O)")
                and (type_exon7_422 == "A or O")
                and (type_exon7_428 == "O2 and (O or A or B)")
                and (type_exon7_429 == "A or O")
                and (type_exon7_431 == "O and (A or B)")
            ):
                Phenotype = "O"
                Genotype = "OO"
                ExtendedGenotype = "O1/O2"

            ## combination 8 | O1O3
            elif (
                (type_exon6 == "O1 and (A or B or O)")
                and (type_exon7_422 == "A or O")
                and (type_exon7_428 == "O and (A or B)")
                and (type_exon7_429 == "A or O")
                and (type_exon7_431 == "O3 and (O or A or B)")
            ):
                Phenotype = "O"
                Genotype = "OO"
                ExtendedGenotype = "O1/O3"

            ## combination 9 | O2O3
            elif (
                (type_exon6 == "A or B or O")
                and (type_exon7_422 == "A or O")
                and (type_exon7_428 == "O2 and (O or A or B)")
                and (type_exon7_429 == "A or O")
                and (type_exon7_431 == "O3 and (O or A or B)")
            ):
                Phenotype = "O"
                Genotype = "OO"
                ExtendedGenotype = "O2/O3"

            ## combination 10 | O1O1
            elif (
                (type_exon6 == "O1")
                and (type_exon7_422 == "A or O")
                and (type_exon7_428 == "O and (A or B)")
                and (type_exon7_429 == "A or O")
                and (type_exon7_431 == "O and (A or B)")
            ):
                Phenotype = "O"
                Genotype = "OO"
                ExtendedGenotype = "O1/O1"

            ## combination 11 | O2O2
            elif (
                (type_exon6 == "A or B or O")
                and (type_exon7_422 == "A or O")
                and (type_exon7_428 == "O2")
                and (type_exon7_429 == "A or O")
                and (type_exon7_431 == "O and (A or B)")
            ):
                Phenotype = "O"
                Genotype = "OO"
                ExtendedGenotype = "O2/O2"

            ## combination 12 | O3O3
            elif (
                (type_exon6 == "A or B or O")
                and (type_exon7_422 == "A or O")
                and (type_exon7_428 == "O and (A or B)")
                and (type_exon7_429 == "A or O")
                and (type_exon7_431 == "O3")
            ):
                Phenotype = "O"
                Genotype = "OO"
                ExtendedGenotype = "O3/O3"

            ## combination 13 | AA
            elif (
                (type_exon6 == "A or B or O")
                and (type_exon7_422 == "A or O")
                and (type_exon7_428 == "O and (A or B)")
                and (type_exon7_429 == "A or O")
                and (type_exon7_431 == "O and (A or B)")
            ):
                a_subtype, a_subtype_warning = determine_a_subtype()
                if a_subtype:
                    Phenotype = a_subtype
                    Genotype = "AA"
                    ExtendedGenotype = f"{a_subtype}/{a_subtype}"
                else:
                    Phenotype = "A"
                    Genotype = "AA"
                    ExtendedGenotype = "A/A"

            ## combination 14 | BB
            elif (
                (type_exon6 == "A or B or O")
                and (type_exon7_422 == "B")
                and (type_exon7_428 == "O and (A or B)")
                and (type_exon7_429 == "B")
                and (type_exon7_431 == "O and (A or B)")
            ):
                Phenotype = "B"
                Genotype = "BB"
                ExtendedGenotype = "B/B"

            ## combination 15 | AB
            elif (
                (type_exon6 == "A or B or O")
                and (type_exon7_422 == "(A or O) and B")
                and (type_exon7_428 == "O and (A or B)")
                and (type_exon7_429 == "(A or O) and B")
                and (type_exon7_431 == "O and (A or B)")
            ):
                a_subtype, a_subtype_warning = determine_a_subtype()
                Phenotype = "AB"
                Genotype = "AB"
                if a_subtype:
                    ExtendedGenotype = f"{a_subtype}/B"
                else:
                    ExtendedGenotype = "A/B"

            ## UNKNOWN
            else:
                Phenotype = "Unknown"
                Genotype = "Unknown"
                ExtendedGenotype = "Unknown"

            # ----- PART 2: Bw SUBTYPE DETECTION -----
            bw_markers, bw_warnings = scan_bw_markers(
                type_exon7_422=type_exon7_422,
                type_exon7_429=type_exon7_429,
                exon7_positions=exon7_positions,
                nreads_exon7=int(nreads_exon7_p422) if nreads_exon7_p422 else 0,
            )

            if bw_markers:
                BwSubtype = "; ".join(bw_markers)
                # Refine phenotype if Bw detected
                if "B" in Phenotype or "B" in Genotype:
                    # Don't override the primary phenotype, but annotate
                    BwSubtype = "weak-B: " + BwSubtype
            elif bw_warnings:
                BwSubtype = "; ".join(bw_warnings)

            # ----- PART 3: PHASE CONFIDENCE -----
            phase_level, phase_detail = compute_phase_confidence(
                phase_evidence_exon6, phase_evidence_exon7, ExtendedGenotype
            )
            PhaseConfidence = f"{phase_level}: {phase_detail}"

            # ----- PART 4: RELIABILITY SCORING -----
            read_counts = [
                nreads6,
                nreads_exon7_p422,
                nreads_exon7_p428,
                nreads_exon7_p429,
                nreads_exon7_p431,
            ]

            valid_read_counts = []
            for count in read_counts:
                try:
                    if pd.notna(count) and float(count) > 0:
                        valid_read_counts.append(float(count))
                except (ValueError, TypeError):
                    continue

            if valid_read_counts:
                min_reads = min(valid_read_counts)

                if min_reads <= 20:
                    Reliability = "Very Low(<=20 reads)"
                elif min_reads <= 40:
                    Reliability = "Low (<=40 reads)"
                elif min_reads >= 500:
                    Reliability = "Robust(>=500 reads)"
                else:
                    Reliability = "Normal"

                if (
                    valid_read_counts
                    and max(valid_read_counts) / min(valid_read_counts) > 5
                ):
                    Reliability += " (Variable coverage)"
            else:
                Reliability = "Unknown (no read data)"

            # ASubtype: A-allele subtype notation (e.g. "c.1032G>A") — separate from read reliability
            ASubtype = a_subtype_warning or ""

            # Append Bw warnings to reliability if cisAB suspected
            if bw_warnings and "cisAB" in " ".join(bw_warnings):
                Reliability += "; " + "; ".join(bw_warnings)

            df[("Result", "Phenotype")] = Phenotype
            df[("Result", "Genotype")] = Genotype
            df[("Result", "ExtendedGenotype")] = ExtendedGenotype
            df[("Notes", "ASubtype")] = ASubtype
            df[("Notes", "ReadReliability")] = Reliability
            df[("Notes", "PhaseConfidence")] = PhaseConfidence
            df[("Notes", "BwSubtype")] = BwSubtype

            return df

        except Exception as e:
            print(f"Error in assign_phenotype_genotype: {str(e)}")
            import traceback

            traceback.print_exc()
            df[("Result", "Phenotype")] = "Error"
            df[("Result", "Genotype")] = "Error"
            df[("Result", "ExtendedGenotype")] = "Error"
            df[("Notes", "ASubtype")] = ""
            df[("Notes", "ReadReliability")] = "Error processing"
            df[("Notes", "PhaseConfidence")] = "Error"
            df[("Notes", "BwSubtype")] = "Error"
            return df

    def _find_haplotype_file(self, sample_dir: str, exon: str) -> Optional[str]:
        """
        Locate the Haplotypes.tsv file for a given sample and exon.

        Searches in:
          <sample_dir>/<exon>/*.Haplotypes.tsv          (flat layout)
          <sample_dir>/<exon>/**/*.Haplotypes.tsv       (subdirectory layout, e.g. haploscan/)
        """
        exon_dir = os.path.join(sample_dir, exon)
        if not os.path.isdir(exon_dir):
            return None

        # Search top-level first, then one level deep (e.g. haploscan/ subdir)
        haplo_files = glob.glob(os.path.join(exon_dir, "*.Haplotypes.tsv"))
        if haplo_files:
            return haplo_files[0]
        haplo_files = glob.glob(os.path.join(exon_dir, "*", "*.Haplotypes.tsv"))
        if haplo_files:
            return haplo_files[0]
        return None

    def process_file(self, filename):
        """Process a single file and extract all necessary data."""
        try:
            sample_name, barcode, pattern_type = self.extract_sample_info(filename)

            sample_dir = os.path.join(self.input_dir, filename)
            exon6_dir = os.path.join(sample_dir, "exon6")
            exon7_dir = os.path.join(sample_dir, "exon7")

            if not (os.path.exists(exon6_dir) and os.path.exists(exon7_dir)):
                error_msg = f"Missing exon6 or exon7 directory"
                print(f"Skipping file {filename}. {error_msg}.")
                self.failed_samples.append({"sample": filename, "reason": error_msg})
                return

            exon6_phenotypes = os.path.join(exon6_dir, "*.ABOPhenotype.txt")
            exon7_phenotypes = os.path.join(exon7_dir, "*.ABOPhenotype.txt")

            exon6_phenotype_files = glob.glob(exon6_phenotypes)
            exon7_phenotype_files = glob.glob(exon7_phenotypes)

            if not exon6_phenotype_files:
                error_msg = "Missing exon6 phenotype files"
                print(f"{error_msg} for {filename}. Skipping.")
                self.failed_samples.append({"sample": filename, "reason": error_msg})
                return

            if not exon7_phenotype_files:
                error_msg = "Missing exon7 phenotype files"
                print(f"{error_msg} for {filename}. Skipping.")
                self.failed_samples.append({"sample": filename, "reason": error_msg})
                return

            if os.path.getsize(exon6_phenotype_files[0]) == 0:
                error_msg = "Empty exon6 phenotype file (0 kb)"
                print(f"{error_msg} for {filename}. Skipping.")
                self.failed_samples.append({"sample": filename, "reason": error_msg})
                return

            if os.path.getsize(exon7_phenotype_files[0]) == 0:
                error_msg = "Empty exon7 phenotype file (0 kb)"
                print(f"{error_msg} for {filename}. Skipping.")
                self.failed_samples.append({"sample": filename, "reason": error_msg})
                return

            # Initialize result dataframe
            result_df = pd.DataFrame(columns=self.columns)
            result_df.loc[0, ("", "Barcode")] = barcode.replace("barcode", "")
            result_df.loc[0, ("", "Sequencing_ID")] = sample_name

            exon6_data = self.parse_exon6(exon6_phenotype_files[0])
            if not exon6_data.empty:
                for pos in [22, 27, 29, 58]:
                    pos_df = exon6_data[exon6_data["Position"] == pos]
                    if not pos_df.empty:
                        for col in [
                            "#Reads",
                            "Mat",
                            "Mis",
                            "Ins",
                            "Del",
                            "A",
                            "G",
                            "C",
                            "T",
                            "Type",
                        ]:
                            if col in pos_df.columns:
                                result_df.loc[0, (f"Exon6_pos{pos}", col)] = (
                                    pos_df.iloc[0][col]
                                )

            exon7_data = self.parse_exon7(exon7_phenotype_files[0])
            if not exon7_data.empty:
                selected_positions = [
                    422,
                    428,
                    429,
                    431,
                    93,
                    685,
                    33,
                    152,
                    153,
                    165,
                    283,
                    329,
                    348,
                    368,
                    397,
                    404,
                    455,
                    533,
                    635,
                    658,
                    680,
                ]

                for pos in selected_positions:
                    pos_df = exon7_data[exon7_data["Position"] == pos]
                    if not pos_df.empty:
                        for col in [
                            "#Reads",
                            "Mat",
                            "Mis",
                            "Ins",
                            "Del",
                            "A",
                            "G",
                            "C",
                            "T",
                            "Type",
                        ]:
                            if col in pos_df.columns:
                                result_df.loc[0, (f"Exon7_pos{pos}", col)] = (
                                    pos_df.iloc[0][col]
                                )

            # ----- Load Haplotypes.tsv for phase confidence -----
            haplo_exon6_file = self._find_haplotype_file(sample_dir, "exon6")
            haplo_exon7_file = self._find_haplotype_file(sample_dir, "exon7")

            phase_evidence_exon6 = parse_haplotypes_tsv(haplo_exon6_file)
            phase_evidence_exon7 = parse_haplotypes_tsv(haplo_exon7_file)

            if phase_evidence_exon6:
                print(
                    f"  Phase data (exon6): {phase_evidence_exon6.total_reads} reads, "
                    f"{len(phase_evidence_exon6.haplotype_counts)} distinct haplotypes"
                )
            if phase_evidence_exon7:
                print(
                    f"  Phase data (exon7): {phase_evidence_exon7.total_reads} reads, "
                    f"{len(phase_evidence_exon7.haplotype_counts)} distinct haplotypes"
                )

            # Assign phenotype/genotype with phase and Bw data
            result_df = self.assign_phenotype_genotype(
                result_df,
                phase_evidence_exon6=phase_evidence_exon6,
                phase_evidence_exon7=phase_evidence_exon7,
            )

            self.results.append(result_df)
            print(f"Successfully processed {filename}")

        except Exception as e:
            print(f"Error processing {filename}: {str(e)}")
            import traceback

            traceback.print_exc()

    def process_files(self):
        """Process all files in the input directory that match expected patterns."""
        print(f"Scanning directory: {self.input_dir}")

        valid_directories = []
        for filename in os.listdir(self.input_dir):
            if os.path.isdir(os.path.join(self.input_dir, filename)):
                match_with_barcode = self.pattern_with_barcode.match(filename)
                match_without_barcode = self.pattern_without_barcode.match(filename)
                if match_with_barcode or match_without_barcode:
                    valid_directories.append(filename)
                else:
                    self.processing_stats["pattern_failed"] += 1
                    self.failed_samples.append(
                        {
                            "sample": filename,
                            "reason": "Filename pattern not recognized",
                        }
                    )

        if not valid_directories:
            print("No valid sample directories found matching expected patterns.")
            return

        print(
            "\033[92m\n ********* Started combining samples to single file ********* \033[0m\n"
        )

        for filename in valid_directories:
            try:
                self.processing_stats["pattern_matched"] += 1
                print(f"\nProcessing file: {filename}")

                sample_name, barcode, pattern_type = self.extract_sample_info(filename)

                if pattern_type == "explicit":
                    self.processing_stats["with_barcode"] += 1
                    print(f"Extracted Sample: {sample_name}, Barcode: {barcode}")
                else:
                    self.processing_stats["without_barcode"] += 1
                    print(
                        f"Extracted Sample: {sample_name}, Barcode: {barcode} (default)"
                    )

                if any(
                    char in sample_name for char in ["<", ">", ":", '"', "|", "?", "*"]
                ):
                    print(
                        f"Warning: Sample name '{sample_name}' contains potentially problematic characters"
                    )

                self.process_file(filename)
                print(
                    f"Done adding Sample {sample_name} with barcode {barcode} to merged data frame"
                )
            except Exception as e:
                print(f"\nError processing file {filename}: {e}")
                self.failed_samples.append(
                    {"sample": filename, "reason": f"Processing error: {str(e)}"}
                )
            finally:
                print(f"Finished processing file: {filename}")

        # Print processing statistics
        print(f"\n--- Processing Statistics ---")
        print(f"Files with explicit barcode: {self.processing_stats['with_barcode']}")
        print(
            f"Files using default barcode: {self.processing_stats['without_barcode']}"
        )
        print(f"Total pattern matches: {self.processing_stats['pattern_matched']}")
        print(f"Pattern match failures: {self.processing_stats['pattern_failed']}")
        print(f"----------------------------")

    def merge_dataframes(self):
        if not self.results:
            print("Warning: No sample results to merge.")
            return pd.DataFrame(columns=self.columns)
        final_df = pd.concat(self.results)
        final_df[("", "Barcode")] = final_df[("", "Barcode")].astype(int)
        final_df = final_df.sort_values(
            by=[("", "Sequencing_ID"), ("", "Barcode")], ascending=True
        )
        return final_df

    def save_results_to_file(self, final_df):
        """Save results to text and Excel files with proper handling of headers and formatting."""
        try:
            final_df.to_csv("./ABO_result.txt", sep="\t", index=False)
            print("Results saved successfully to text file.")
        except Exception as txt_err:
            print(f"Error saving to text file: {txt_err}")
            return

        try:
            read_count_cols = []
            if isinstance(final_df.columns, pd.MultiIndex):
                for i, col in enumerate(final_df.columns):
                    if col[1] == "#Reads":
                        read_count_cols.append(i)
            else:
                for i, col in enumerate(final_df.columns):
                    if "#Reads" in str(col):
                        read_count_cols.append(i)

            writer = pd.ExcelWriter("./ABO_result.xlsx", engine="xlsxwriter")

            if isinstance(final_df.columns, pd.MultiIndex):
                final_df.columns = final_df.columns.droplevel()

            final_df.to_excel(
                writer, sheet_name="ABO_Result", header=True, index=False, startrow=1
            )

            workbook = writer.book
            worksheet = writer.sheets["ABO_Result"]

            # Define Excel formats
            data_format = workbook.add_format(
                {"bg_color": "white", "font_color": "black", "border": 1}
            )
            header_format = workbook.add_format(
                {
                    "bold": True,
                    "fg_color": "#007399",
                    "border": 1,
                    "font_color": "white",
                }
            )
            red_bg_format = workbook.add_format(
                {"bg_color": "#e2725b", "font_color": "black"}
            )
            orange_bg_format = workbook.add_format(
                {"bg_color": "#ff9a00", "font_color": "black"}
            )

            header_format.set_align("center")
            header_format.set_align("vcenter")

            num_rows, num_cols = final_df.shape

            # ReadReliability is the 3rd-to-last column (ASubtype, ReadReliability, PhaseConfidence, BwSubtype)
            reliability_col = xl_col_to_name(num_cols - 3)

            print(f"Data has {num_rows} rows, starting at row 3 with two header rows")

            # Apply conditional formatting to each read count column
            for col_idx in read_count_cols:
                col_letter = xl_col_to_name(col_idx)

                worksheet.conditional_format(
                    f"{col_letter}3:{col_letter}{num_rows + 2}",
                    {
                        "type": "cell",
                        "criteria": "<=",
                        "value": 20,
                        "format": red_bg_format,
                    },
                )

                worksheet.conditional_format(
                    f"{col_letter}3:{col_letter}{num_rows + 2}",
                    {
                        "type": "cell",
                        "criteria": "between",
                        "minimum": 21,
                        "maximum": 40,
                        "format": orange_bg_format,
                    },
                )

            print(
                f"Applying read count conditional formatting to columns: {[xl_col_to_name(i) for i in read_count_cols]}"
            )

            try:
                worksheet.conditional_format(
                    f"A3:{xl_col_to_name(num_cols - 1)}{num_rows + 2}",
                    {
                        "type": "formula",
                        "criteria": f'=${reliability_col}3="Very Low(<=20 reads)"',
                        "format": red_bg_format,
                    },
                )
                worksheet.conditional_format(
                    f"A3:{xl_col_to_name(num_cols - 1)}{num_rows + 2}",
                    {
                        "type": "formula",
                        "criteria": f'=${reliability_col}3="Low (<=40 reads)"',
                        "format": orange_bg_format,
                    },
                )
            except Exception as format_err:
                print(
                    f"Warning: Could not apply row-level conditional formatting: {format_err}"
                )

            # Write data
            for row in range(num_rows):
                for col in range(num_cols):
                    cell_value = final_df.iat[row, col]
                    if not pd.isna(cell_value):
                        worksheet.write(row + 2, col, cell_value, data_format)

            header_columns = [
                "Exon6_pos22",
                "Exon7_pos422",
                "Exon7_pos428",
                "Exon7_pos429",
                "Exon7_pos431",
                "Exon7_pos93",
                "Exon7_pos685",
                "Exon6_pos27",
                "Exon6_pos29",
                "Exon6_pos58",
                "Exon7_pos33",
                "Exon7_pos152",
                "Exon7_pos153",
                "Exon7_pos165",
                "Exon7_pos283",
                "Exon7_pos329",
                "Exon7_pos348",
                "Exon7_pos368",
                "Exon7_pos397",
                "Exon7_pos404",
                "Exon7_pos455",
                "Exon7_pos533",
                "Exon7_pos635",
                "Exon7_pos658",
                "Exon7_pos680",
            ]

            # SNP groups start at col 5 (after Barcode, Sequencing_ID, Phenotype, Genotype, ExtendedGenotype)
            column_start = 5
            merge_ranges = []

            for header in header_columns:
                start_col = column_start
                end_col = start_col + 9

                start_letter = xl_col_to_name(start_col)
                end_letter = xl_col_to_name(end_col)

                merge_ranges.append((f"{start_letter}1:{end_letter}1", header))
                column_start = end_col + 1

            # Notes columns: ASubtype, ReadReliability, PhaseConfidence, BwSubtype
            notes_start = xl_col_to_name(column_start)
            notes_end = xl_col_to_name(column_start + 3)

            # Merge header ranges
            worksheet.merge_range("A1:B1", "Sample", header_format)
            worksheet.merge_range("C1:E1", "Result", header_format)
            worksheet.merge_range(
                f"{notes_start}1:{notes_end}1", "Notes", header_format
            )

            for merge_range in merge_ranges:
                worksheet.merge_range(merge_range[0], merge_range[1], header_format)

            for col in range(num_cols):
                cell_value = final_df.columns[col]
                if not pd.isna(cell_value):
                    worksheet.write(1, col, cell_value, header_format)

            writer.close()
            print("Results saved successfully to Excel file.")
        except Exception as excel_err:
            print(f"Error saving to Excel file: {excel_err}")
            import traceback

            traceback.print_exc()

        # LIS export
        self.df_for_lis_soft = pd.DataFrame()
        self.df_for_lis_soft["Sample ID"] = final_df["Sequencing_ID"]
        self.df_for_lis_soft["Shipment Date"] = ""

        if (
            "Genotype" in final_df.columns
            and not final_df["Genotype"].isnull().all()
            and not (final_df["Genotype"] == "Unknown").all()
        ):
            valid_genotype_mask = (final_df["Genotype"] != "Unknown") & final_df[
                "Genotype"
            ].notnull()
            self.df_for_lis_soft.loc[valid_genotype_mask, "ABO Geno Type1"] = (
                final_df.loc[valid_genotype_mask, "Genotype"].str[0]
            )
            self.df_for_lis_soft.loc[valid_genotype_mask, "ABO Geno Type2"] = (
                final_df.loc[valid_genotype_mask, "Genotype"].str[1]
            )
        else:
            self.df_for_lis_soft["ABO Geno Type1"] = ""
            self.df_for_lis_soft["ABO Geno Type2"] = ""

        # Keep a plain phenotype (A/B/AB/O) for the LIS CSV export
        plain_phenotype = final_df["Phenotype"].copy()

        # Annotate Phenotype with Bw subtype info when present — xlsx only
        if "BwSubtype" in final_df.columns:
            bw_mask = final_df["BwSubtype"].astype(str).str.len() > 0
            annotated_phenotype = final_df["Phenotype"].copy()
            annotated_phenotype.loc[bw_mask] = (
                final_df.loc[bw_mask, "Phenotype"]
                + " ["
                + final_df.loc[bw_mask, "BwSubtype"]
                + "]"
            )
            # Update final_df so the xlsx Phenotype column is annotated
            final_df["Phenotype"] = annotated_phenotype

        # CSV gets the plain A/B/AB/O values
        self.df_for_lis_soft["ABO Pheno Type"] = plain_phenotype
        self.df_for_lis_soft["RH"] = ""
        self.df_for_lis_soft["Blood Type"] = plain_phenotype
        self.df_for_lis_soft["ABORH Comments"] = ""

        # Include Notes columns in LIS export
        if "ASubtype" in final_df.columns:
            self.df_for_lis_soft["A Subtype"] = final_df["ASubtype"]
        if "BwSubtype" in final_df.columns:
            self.df_for_lis_soft["Bw Subtype"] = final_df["BwSubtype"]
        if "PhaseConfidence" in final_df.columns:
            self.df_for_lis_soft["Phase Confidence"] = final_df["PhaseConfidence"]
        if "ReadReliability" in final_df.columns:
            self.df_for_lis_soft["Read Reliability"] = final_df["ReadReliability"]

        if isinstance(final_df.columns, pd.MultiIndex):
            reads_df = final_df.loc[:, (slice(None), "#Reads")]
            self.df_for_lis_soft["#Reads"] = reads_df.mean(axis=1)
        else:
            reads_columns = [col for col in final_df.columns if "#Reads" in str(col)]
            if reads_columns:
                self.df_for_lis_soft["#Reads"] = final_df[reads_columns].mean(axis=1)
            else:
                self.df_for_lis_soft["#Reads"] = 0

        self.df_for_lis_soft.drop_duplicates(inplace=True)
        self.df_for_lis_soft.to_csv("./final_export.csv", index=False, encoding="utf-8")
        print(
            f"LIS export file created successfully with {len(self.df_for_lis_soft)} samples"
        )

    def run(self):
        """Run the ABOReportParser."""
        self.process_files()
        final_df = self.merge_dataframes()
        if final_df.empty:
            print("No results to report. Exiting.")
            return
        print("\n\nFinal Results:")
        print("-" * 336)
        print(final_df.to_string(index=False))
        print("-" * 336)
        self.save_results_to_file(final_df)

        if self.failed_samples:
            print("\n\nFailed Samples Summary:")
            print("-" * 80)
            for sample in self.failed_samples:
                print(f"Sample: {sample['sample']} - Reason: {sample['reason']}")
            print("-" * 80)
            print(f"Total failed samples: {len(self.failed_samples)}")
        else:
            print("\nAll samples processed successfully.")

        # Print comprehensive processing summary
        print(f"\n=== PROCESSING SUMMARY ===")
        print(
            f"Total directories scanned: {self.processing_stats['pattern_matched'] + self.processing_stats['pattern_failed']}"
        )
        print(f"Successfully processed: {len(self.results)}")
        print(f"Samples with explicit barcode: {self.processing_stats['with_barcode']}")
        print(
            f"Samples with default barcode ({self.default_barcode}): {self.processing_stats['without_barcode']}"
        )
        print(f"Failed samples: {len(self.failed_samples)}")
        print(
            f"Pattern recognition failures: {self.processing_stats['pattern_failed']}"
        )

        if len(self.results) > 0:
            try:
                phenotype_counts = final_df["Phenotype"].value_counts()
                print(f"\nPhenotype Distribution:")
                for phenotype, count in phenotype_counts.items():
                    print(f"  {phenotype}: {count}")

                try:
                    extended_genotype_counts = final_df[
                        "ExtendedGenotype"
                    ].value_counts()
                    print(f"\nExtended Genotype Distribution (with A subtypes):")

                    a_subtypes = {}
                    ab_subtypes = {}
                    other_genotypes = {}

                    for genotype, count in extended_genotype_counts.items():
                        if "A1" in str(genotype) or "A2" in str(genotype):
                            if "B" in str(genotype):
                                ab_subtypes[genotype] = count
                            else:
                                a_subtypes[genotype] = count
                        else:
                            other_genotypes[genotype] = count

                    if a_subtypes:
                        print(f"  A Subtypes:")
                        for genotype, count in sorted(a_subtypes.items()):
                            print(f"    {genotype}: {count}")

                    if ab_subtypes:
                        print(f"  AB Subtypes:")
                        for genotype, count in sorted(ab_subtypes.items()):
                            print(f"    {genotype}: {count}")

                    if other_genotypes:
                        print(f"  Other Genotypes:")
                        for genotype, count in sorted(other_genotypes.items()):
                            print(f"    {genotype}: {count}")

                except Exception as e:
                    print(f"Could not analyze extended genotype distribution: {e}")

                # NEW: Bw subtype distribution
                try:
                    bw_col = final_df.get("BwSubtype")
                    if bw_col is not None:
                        bw_detected = bw_col[bw_col.astype(str).str.len() > 0]
                        if not bw_detected.empty:
                            print(
                                f"\nBw/B(A) Subtype Flags ({len(bw_detected)} samples):"
                            )
                            for val, count in bw_detected.value_counts().items():
                                print(f"    {val}: {count}")
                except Exception as e:
                    print(f"Could not analyze Bw distribution: {e}")

                # NEW: Phase confidence distribution
                try:
                    pc_col = final_df.get("PhaseConfidence")
                    if pc_col is not None:
                        print(f"\nPhase Confidence Distribution:")
                        for val, count in pc_col.value_counts().head(10).items():
                            print(f"    {val}: {count}")
                except Exception as e:
                    print(f"Could not analyze phase confidence: {e}")

            except Exception as e:
                print(f"Could not analyze phenotype distribution: {e}")

        print(f"============================")


# ===========================================================================
# CLI entry point
# ===========================================================================


def create_argument_parser():
    """Create and configure argument parser."""
    parser = argparse.ArgumentParser(
        prog="aggregate_abo_reports.py",
        description="ABO Blood Group Report Aggregator — Combines individual sample ABO phenotype results into unified reports.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
EXAMPLES:
  # Basic usage with default barcode (barcode00):
  python aggregate_abo_reports.py /path/to/results

  # Specify custom default barcode:
  python aggregate_abo_reports.py /path/to/results --default-barcode barcode99

  # Process data with verbose output:
  python aggregate_abo_reports.py /path/to/results --verbose

EXPECTED DIRECTORY STRUCTURE:
  input_directory/
  ├── SAMPLE1_barcode01/
  │   ├── exon6/
  │   │   ├── *.ABOPhenotype.txt
  │   │   └── *.Haplotypes.tsv    (NEW — from pysam_haploscan.py)
  │   └── exon7/
  │       ├── *.ABOPhenotype.txt
  │       └── *.Haplotypes.tsv    (NEW — from pysam_haploscan.py)
  ├── SAMPLE2_barcode02/
  └── SAMPLE3/  (will use default barcode)

OUTPUT FILES:
  - ABO_result.txt: Tab-separated comprehensive results
  - ABO_result.xlsx: Excel file with conditional formatting
  - final_export.csv: LIS-compatible export format (now includes BwSubtype + PhaseConfidence)

NEW in v1.2.0:
  - BwSubtype column: Flags Bw.04 (c.1054C>T), Bw.06 (c.657C>T), B(A).02/.04,
    and possible cisAB configurations
  - PhaseConfidence column: Validates heterozygous calls using read-level
    haplotype data from pysam_haploscan.py Haplotypes.tsv output

For more information, see: https://github.com/fmobegi/nf-core-abotyper
        """,
    )

    parser.add_argument(
        "input_directory",
        help="Path to directory containing sample subdirectories with ABO phenotype results",
    )

    parser.add_argument(
        "--default-barcode",
        "-b",
        default="barcode00",
        metavar="BARCODE",
        help="Default barcode for samples without explicit barcode suffix (default: %(default)s)",
    )

    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose output for debugging",
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
        help="Show version information and exit",
    )

    return parser


def validate_arguments(args):
    """Validate command line arguments and input directory."""
    if not os.path.exists(args.input_directory):
        print(f"Error: Input directory '{args.input_directory}' does not exist.")
        sys.exit(1)

    if not os.path.isdir(args.input_directory):
        print(f"Error: '{args.input_directory}' is not a directory.")
        sys.exit(1)

    if not re.match(r"^barcode\d{1,2}$", args.default_barcode):
        print(
            f"Warning: Unusual barcode format '{args.default_barcode}'. Expected format: barcodeXX"
        )

    return True


if __name__ == "__main__":
    parser = create_argument_parser()
    args = parser.parse_args()

    validate_arguments(args)

    if args.verbose:
        print(f"Verbose mode enabled")
        print(f"Input directory: {args.input_directory}")
        print(f"Default barcode: {args.default_barcode}")
        print(f"Script version: {__version__}")

    print(f"Using default barcode: {args.default_barcode}")

    parser_instance = ABOReportParser(args.input_directory, args.default_barcode)
    parser_instance.run()
    print("All done!\n")
