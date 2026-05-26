/*
 * Subworkflow: minimap_align_exons
 * Description: Aligns each sample to two exon references using Minimap2,
 * then runs Samtools coverage, flagstat, and stats.
 * Uses a pure metadata-driven approach with no aliasing or premature splitting.
 */

include { MINIMAP2_ALIGN    } from '../../../modules/nf-core/minimap2/align'
include { SAMTOOLS_COVERAGE } from '../../../modules/nf-core/samtools/coverage'
include { SAMTOOLS_FLAGSTAT } from '../../../modules/nf-core/samtools/flagstat'
include { SAMTOOLS_STATS    } from '../../../modules/nf-core/samtools/stats'

workflow MINIMAP2_ALIGN_READS {
    take:
    ch_combined_input // channel: [ val(meta), [ fastq ] ] - with exon metadata
    ch_combined_fasta // channel: [ val(meta), [ fasta ] ] - with exon metadata
    ch_combined_fai // channel: [ val(meta), [ fai ] ] - with exon metadata

    main:

    ch_versions = channel.empty()

    // Simple combine and filter approach to match samples with correct exon references
    ch_minimap_ready = ch_combined_input
        .combine(ch_combined_fasta)
        .filter { sample_meta, fastq, ref_meta, fasta ->
            sample_meta.exon == ref_meta.exon
        }
        .map { sample_meta, fastq, ref_meta, fasta ->
            def combined_meta = sample_meta + [ref_id: ref_meta.id]
            [combined_meta, fastq, fasta]
        }

    /*
    MODULE: MINIMAP2_ALIGN
    */
    MINIMAP2_ALIGN(
        ch_minimap_ready.map { meta, fastq, fasta -> [meta, fastq] },
        ch_minimap_ready.map { meta, fastq, fasta -> [meta, fasta] },
        "bam",
        "bai",
        false,
        false,
    )

    // Join BAM and BAI files for processes that need both
    ch_bam_bai = MINIMAP2_ALIGN.out.bam.join(MINIMAP2_ALIGN.out.index, by: 0)

    /*
    MODULE: SAMTOOLS_COVERAGE
    */
    ch_coverage_input = ch_bam_bai
        .combine(ch_combined_fasta)
        .combine(ch_combined_fai)
        .filter { bam_meta, bam, bai, fasta_meta, fasta, fai_meta, fai ->
            bam_meta.exon == fasta_meta.exon && bam_meta.exon == fai_meta.exon
        }

    SAMTOOLS_COVERAGE(
        ch_coverage_input.map { bam_meta, bam, bai, fasta_meta, fasta, fai_meta, fai ->
            [bam_meta, bam, bai]
        },
        ch_coverage_input.map { bam_meta, bam, bai, fasta_meta, fasta, fai_meta, fai ->
            [fasta_meta, fasta, fai]
        }
    )

    /*
    MODULE: SAMTOOLS_FLAGSTAT
    */
    SAMTOOLS_FLAGSTAT(
        ch_bam_bai
    )

    /*
    MODULE: SAMTOOLS_STATS
    */
    ch_stats_input = ch_bam_bai
        .combine(ch_combined_fasta)
        .combine(ch_combined_fai)  // ADD THIS
        .filter { bam_meta, bam, bai, fasta_meta, fasta, fai_meta, fai ->
            bam_meta.exon == fasta_meta.exon && bam_meta.exon == fai_meta.exon
        }

    SAMTOOLS_STATS(
        ch_stats_input.map { bam_meta, bam, bai, fasta_meta, fasta, fai_meta, fai ->
            [bam_meta, bam, bai]
        },
        ch_stats_input.map { bam_meta, bam, bai, fasta_meta, fasta, fai_meta, fai ->
            [fasta_meta, fasta, fai]
        }
    )

    emit:
    bam      = MINIMAP2_ALIGN.out.bam // channel: [ val(meta), path(bam) ]
    bai      = MINIMAP2_ALIGN.out.index // channel: [ val(meta), path(bai) ]
    coverage = SAMTOOLS_COVERAGE.out.coverage // channel: [ val(meta), path(txt) ]
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat // channel: [ val(meta), path(flagstat) ]
    stats    = SAMTOOLS_STATS.out.stats // channel: [ val(meta), path(stats) ]
    fasta    = ch_combined_fasta // channel: [ val(meta1), path(fasta) ]
    fai      = ch_combined_fai // channel: [ val(meta1), path(fai) ]
    versions = ch_versions // channel: [ versions.yml ] | All modules here are using topics
}
