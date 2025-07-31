include { SAMTOOLS_MPILEUP   } from '../../../modules/nf-core/samtools/mpileup/main'
include { MPILEUP_NUCL_FREQ  } from '../../../modules/local/mpileupstats/main'

workflow VARIANTS_QUANTIFICATION {

    take:
    ch_bam       // channel: [ val(meta), [ bam ] ] - with exon metadata
    ch_bai       // channel: [ val(meta), [ bai ] ] - with exon metadata
    ch_fasta     // channel: [ val(meta), [ fasta ] ] - with exon metadata
    ch_fai       // channel: [ val(meta), [ fai ] ] - with exon metadata
    ch_exon6_bed // channel: [ val(meta), [ bed ] ] - pre-computed bed file
    ch_exon7_bed // channel: [ val(meta), [ bed ] ] - pre-computed bed file

    main:

    ch_versions = Channel.empty()

    // Add exon metadata to bed files to match with samples
    ch_exon6_bed_with_meta = ch_exon6_bed
        .map { meta, bed -> 
            def new_meta = meta + [exon: 'exon6']
            [new_meta, bed]
        }
    
    ch_exon7_bed_with_meta = ch_exon7_bed
        .map { meta, bed -> 
            def new_meta = meta + [exon: 'exon7']
            [new_meta, bed]
        }

    ch_combined_bed = ch_exon6_bed_with_meta.mix(ch_exon7_bed_with_meta)

    // Prepare mpileup input by combining bam with bed files based on exon metadata
    ch_mpileup_input = ch_bam
        .combine(ch_combined_bed)
        .filter { bam_meta, bam, bed_meta, bed -> 
            bam_meta.exon == bed_meta.exon 
        }
        .map { bam_meta, bam, bed_meta, bed -> 
            [bam_meta, bam, bed]
        }

    
    /*
    Some sanity check to ensure paths and metadata are as expected by the processes below. 
    SAMtools is the main offender. Modules parameterization keeps fluctuating between samtools/subtool 
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    // @TODO NOTE TO SELF:  UNCOMMENT TO EXECUTE FOR TESTING
    // ch_mpileup_input.view { "MPILEUP input: $it" }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    //
    // MODULE: samtools/mpileup
    //
    SAMTOOLS_MPILEUP (
        ch_mpileup_input,
        ch_fasta
    )
    
    ch_versions = ch_versions.mix(SAMTOOLS_MPILEUP.out.versions.first())

    // 
    // MODULE: mpileupmetrics
    // 
    MPILEUP_NUCL_FREQ (
        SAMTOOLS_MPILEUP.out.mpileup,
        ch_fasta
    )

    ch_versions = ch_versions.mix(MPILEUP_NUCL_FREQ.out.versions.first())

    emit:
    metrics       = MPILEUP_NUCL_FREQ.out.tsv              // channel: [ val(meta), [ tsv ] ] - with exon metadata
    versions      = ch_versions                            // channel: [ versions.yml ]
}
