include { SAMTOOLS_MPILEUP   } from '../../../modules/nf-core/samtools/mpileup/main'
include { MPILEUP_NUCL_FREQ  } from '../../../modules/local/mpileupstats/main'

workflow VARIANTS_QUANTIFICATION {

    take:
    ch_bam       // channel: [ val(meta), [ bam ] ] - with exon metadata
    ch_bai       // channel: [ val(meta), [ bai ] ] - with exon metadata
    ch_fasta     // channel: [ val(meta), [ fasta ] ] - with exon metadata
    ch_fai       // channel: [ val(meta), [ fai ] ] - with exon metadata
    ch_bed       // channel: [ val(meta), [ bed ] ] - combined bed files with exon metadata

    main:

    ch_versions = Channel.empty()

    // Prepare mpileup input by matching bam with bed files based on exon metadata
    ch_mpileup_input = ch_bam
        .combine(ch_bed)
        .filter { bam_meta, bam, bed_meta, bed ->
            bam_meta.exon == bed_meta.exon
        }
        .map { bam_meta, bam, bed_meta, bed ->
            [bam_meta, bam, bed]
        }

    // Match fasta files with samples based on exon
    ch_fasta_matched = ch_mpileup_input
        .combine(ch_fasta)
        .filter { bam_meta, bam, bed, fasta_meta, fasta ->
            bam_meta.exon == fasta_meta.exon
        }
        .map { bam_meta, bam, bed, fasta_meta, fasta ->
            [bam_meta, bam, bed, fasta_meta, fasta]
        }

    /*
    Some sanity check to ensure paths and metadata are as expected by the processes below.
    SAMtools is the main offender. Modules parameterization keeps fluctuating between samtools/subtool
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    // @TODO NOTE TO SELF:  UNCOMMENT TO EXECUTE FOR TESTING
    // ch_fasta_matched.view { "MPILEUP input: $it" }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    //
    // MODULE: samtools/mpileup
    // Input: tuple val(meta), path(input), path(intervals)
    //        tuple val(meta2), path(fasta)
    //
    SAMTOOLS_MPILEUP (
        ch_fasta_matched.map { bam_meta, bam, bed, fasta_meta, fasta -> [bam_meta, bam, bed] },
        ch_fasta_matched.map { bam_meta, bam, bed, fasta_meta, fasta -> [fasta_meta, fasta] }
    )

    ch_versions = ch_versions.mix(SAMTOOLS_MPILEUP.out.versions.first())

    // Match mpileup output with fasta for nucleotide frequency analysis
    ch_nucl_freq_input = SAMTOOLS_MPILEUP.out.mpileup
        .combine(ch_fasta)
        .filter { mpileup_meta, mpileup, fasta_meta, fasta ->
            mpileup_meta.exon == fasta_meta.exon
        }

    //
    // MODULE: mpileupmetrics
    //
    MPILEUP_NUCL_FREQ (
        ch_nucl_freq_input.map { mpileup_meta, mpileup, fasta_meta, fasta -> [mpileup_meta, mpileup] },
        ch_nucl_freq_input.map { mpileup_meta, mpileup, fasta_meta, fasta -> [fasta_meta, fasta] }
    )

    ch_versions = ch_versions.mix(MPILEUP_NUCL_FREQ.out.versions.first())

    emit:
    metrics       = MPILEUP_NUCL_FREQ.out.tsv              // channel: [ val(meta), [ tsv ] ] - with exon metadata
    versions      = ch_versions                            // channel: [ versions.yml ]
}
