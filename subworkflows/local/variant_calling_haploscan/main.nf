/*
  SUBWORKFLOW: VARIANTS_QUANTIFICATION (haploscan)

  Drop-in replacement for variant_calling_mpileup.
  Uses pysam_haploscan.py to compute per-position allele frequencies AND
  per-read haplotypes in a single pass over the aligned BAM — no pileup file.

  Channel matching logic mirrors variant_calling_mpileup exactly so that
  the same metadata quirks (fai id has .fasta suffix, bam carries ref_id)
  are handled identically.
*/

include { HAPLOSCAN } from '../../../modules/local/haploscan/main'

workflow VARIANTS_QUANTIFICATION {
    take:
    ch_bam   // channel: [ val(meta),  path(bam)   ]
    ch_bai   // channel: [ val(meta),  path(bai)   ]
    ch_fasta // channel: [ val(meta1), path(fasta) ]
    ch_fai   // channel: [ val(meta1), path(fai)   ]
    ch_bed   // channel: [ val(meta1), path(bed)   ]  — accepted but not used

    main:

    // Join BAM and BAI (metadata matches exactly)
    ch_bam_bai = ch_bam.join(ch_bai, by: 0)

    // Match FASTA and FAI — fai metadata id has .fasta suffix, fasta does not
    ch_fasta_fai = ch_fasta
        .combine(ch_fai)
        .filter { fasta_meta, fasta, fai_meta, fai ->
            fasta_meta.exon == fai_meta.exon &&
            fasta_meta.id   == fai_meta.id.replace('.fasta', '')
        }
        .map { fasta_meta, fasta, fai_meta, fai ->
            [fasta_meta, fasta, fai]
        }

    // Match BAM+BAI with FASTA+FAI — bam carries ref_id, fasta carries id
    ch_haploscan_input = ch_bam_bai
        .combine(ch_fasta_fai)
        .filter { bam_meta, bam, bai, fasta_meta, fasta, fai ->
            bam_meta.exon   == fasta_meta.exon &&
            bam_meta.ref_id == fasta_meta.id
        }

    HAPLOSCAN(
        ch_haploscan_input.map { bam_meta, bam, bai, fasta_meta, fasta, fai ->
            [bam_meta, bam, bai]
        },
        ch_haploscan_input.map { bam_meta, bam, bai, fasta_meta, fasta, fai ->
            [fasta_meta, fasta, fai]
        },
    )

    emit:
    metrics    = HAPLOSCAN.out.tsv        // channel: [ val(meta), path(*.AlignmentStatistics.tsv) ]
    haplotypes = HAPLOSCAN.out.haplotypes // channel: [ val(meta), path(*.Haplotypes.tsv) ]
}

