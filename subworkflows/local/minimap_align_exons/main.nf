/*
 * Subworkflow: minimap_align_exons
 * Description: Aligns each sample to two exon references using Minimap2, 
 * then runs Samtools coverage, flagstat, and stats.
 * Uses a pure metadata-driven approach with no aliasing or premature splitting.
 */

include { MINIMAP2_ALIGN    } from '../../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_COVERAGE } from '../../../modules/nf-core/samtools/coverage/main'
include { SAMTOOLS_FLAGSTAT } from '../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_STATS    } from '../../../modules/nf-core/samtools/stats/main'

workflow MINIMAP2_ALIGN_READS {

    take:
    ch_combined_input   // channel: [ val(meta), [ fastq ] ] - with exon metadata
    ch_combined_fasta   // channel: [ val(meta), [ fasta ] ] - with exon metadata  
    ch_combined_fai     // channel: [ val(meta), [ fai ] ] - with exon metadata

    main:

    ch_versions = Channel.empty()

    // 
    // MODULE: Minimap2/align
    // 
    MINIMAP2_ALIGN (
        ch_combined_input,
        ch_combined_fasta,
        bam_format="bam",
        bam_index_extension="bai",
        cigar_paf_format=false,
        cigar_bam=false
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    //
    // MODULE: Samtools/coverage 
    // 
    SAMTOOLS_COVERAGE (
        MINIMAP2_ALIGN.out.bam
            .join(MINIMAP2_ALIGN.out.index),
        ch_combined_fasta,
        ch_combined_fai
    )
    ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE.out.versions)

    //
    // MODULE: Samtools/flagstat 
    // 
    SAMTOOLS_FLAGSTAT (
        MINIMAP2_ALIGN.out.bam
            .join(MINIMAP2_ALIGN.out.index)
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)

    //
    // MODULE: Samtools/stats 
    // 
    SAMTOOLS_STATS (
        MINIMAP2_ALIGN.out.bam
            .join(MINIMAP2_ALIGN.out.index),
        ch_combined_fasta
    )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)

    emit:
    bam           = MINIMAP2_ALIGN.out.bam                  // channel: [ val(meta), [ bam ] ] - with exon metadata
    bai           = MINIMAP2_ALIGN.out.index                // channel: [ val(meta), [ bai ] ] - with exon metadata
    coverage      = SAMTOOLS_COVERAGE.out.coverage         // channel: [ val(meta), [ txt ] ] - with exon metadata
    fasta         = ch_combined_fasta                       // channel: [ val(meta), [ fasta ] ] - with exon metadata
    fai           = ch_combined_fai                         // channel: [ val(meta), [ fai ] ] - with exon metadata
    versions      = ch_versions                             // channel: [ versions.yml ]
}