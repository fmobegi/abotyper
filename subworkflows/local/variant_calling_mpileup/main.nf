/*
  SUBWORKFLOW: VARIANTS_QUANTIFICATION
*/

//   Description: 
//     Performs variant quantification by generating mpileup files from aligned BAM files
//     and computing nucleotide frequency statistics at specific genomic positions.
//   
//   Workflow Steps:
//     1. Match BAM/BAI files with BED interval files based on exon and reference ID
//     2. Run SAMTOOLS_MPILEUP to generate pileup data for variant positions
//     3. Analyze pileup output with MPILEUPSTATS to compute nucleotide frequencies

include { SAMTOOLS_MPILEUP } from '../../../modules/nf-core/samtools/mpileup/main'
include { MPILEUPSTATS     } from '../../../modules/local/mpileupstats/main'

workflow VARIANTS_QUANTIFICATION {
    take:
    ch_bam   // channel: [ val(meta), path(bam)]
    ch_bai   // channel: [ val(meta), path(bai)]
    ch_fasta // channel: [ val(meta1), path(fasta)]
    ch_fai   // channel: [ val(meta1), path(fai)]
    ch_bed   // channel: [ val(meta1), path(bed)]

    main:

    // Join BAM and BAI files (these have matching metadata)
    ch_bam_bai = ch_bam.join(ch_bai, by: 0)

    // Match BAM+BAI with BED files based on exon and reference ID
    ch_mpileup_input = ch_bam_bai
        .combine(ch_bed)
        .filter { bam_meta, bam, bai, bed_meta, bed ->
            bam_meta.exon == bed_meta.exon && bam_meta.ref_id == bed_meta.id.replace('.fasta', '')
        }
        .map { bam_meta, bam, bai, bed_meta, bed ->
            [bam_meta, bam, bai, bed]
        }

    // Match FASTA and FAI using combine+filter (because metadata IDs don't match)
    ch_fasta_fai = ch_fasta
        .combine(ch_fai)
        .filter { fasta_meta, fasta, fai_meta, fai ->
            // Match by exon and normalize the IDs by removing .fasta extension
            fasta_meta.exon == fai_meta.exon && 
            fasta_meta.id == fai_meta.id.replace('.fasta', '')
        }
        .map { fasta_meta, fasta, fai_meta, fai ->
            [fasta_meta, fasta, fai]  // Keep FASTA metadata (without extension)
        }

    // Match BAM+BAI+BED with FASTA+FAI based on exon and reference ID
    ch_fasta_matched = ch_mpileup_input
        .combine(ch_fasta_fai)
        .filter { bam_meta, bam, bai, bed, fasta_meta, fasta, fai ->
            bam_meta.exon == fasta_meta.exon && bam_meta.ref_id == fasta_meta.id
        }
        .map { bam_meta, bam, bai, bed, fasta_meta, fasta, fai ->
            [bam_meta, bam, bai, bed, fasta_meta, fasta, fai]
        }

    /*
    MODULE: SAMTOOLS_MPILEUP
    */
    SAMTOOLS_MPILEUP(
        ch_fasta_matched.map { bam_meta, bam, bai, bed, fasta_meta, fasta, fai -> 
            [bam_meta, bam, bai, bed] 
        },
        ch_fasta_matched.map { bam_meta, bam, bai, bed, fasta_meta, fasta, fai -> 
            [fasta_meta, fasta, fai] 
        }
    )

    // Match mpileup output with fasta for nucleotide frequency analysis
    ch_nucl_freq_input = SAMTOOLS_MPILEUP.out.mpileup
        .combine(ch_fasta)
        .filter { mpileup_meta, mpileup, fasta_meta, fasta ->
            mpileup_meta.exon == fasta_meta.exon && mpileup_meta.ref_id == fasta_meta.id
        }

    /*
    MODULE: MPILEUPSTATS
    */
    MPILEUPSTATS(
        ch_nucl_freq_input.map { mpileup_meta, mpileup, fasta_meta, fasta -> 
            [mpileup_meta, mpileup] 
        },
        ch_nucl_freq_input.map { mpileup_meta, mpileup, fasta_meta, fasta -> 
            [fasta_meta, fasta] 
        }
    )

    emit:
    metrics  = MPILEUPSTATS.out.tsv // channel: [ val(meta), path(tsv) ]
}
