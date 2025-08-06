include { GETABOSNPS        } from '../../../modules/local/abo/abosnps/main'
include { ABOSNPS2PHENO     } from '../../../modules/local/abo/snps2pheno/main'

workflow PREDICTABOPHENOTYPE {

    take:
    ch_variants_freq    // channel: [ val(meta), [ freq ] ] - with exon metadata
    ch_bam_coverage     // channel: [ val(meta), [ cov ] ] - with exon metadata

    main:

    ch_versions = Channel.empty()

    // Join variants frequency with coverage based on matching metadata
    ch_combined_input = ch_variants_freq
        .join(ch_bam_coverage)
        .map { meta, freq, cov ->
            [meta, freq, cov, meta.exon]  // Pass exon as parameter
        }

    // Check channels for sanity
    // ch_combined_input.view { meta, freq, cov, exon -> "Combined: meta=$meta, freq=$freq, cov=$cov, exon=$exon" }

    GETABOSNPS ( ch_combined_input )

    ch_versions = ch_versions.mix(GETABOSNPS.out.versions)

    // Prepare SNP reports for ABOSNPS2PHENO
    ch_SNP_reports = GETABOSNPS.out.phenotype
        .map { meta, file ->
            [meta.id, [exon: meta.exon, file: file]]
        }
        .groupTuple()
        .map { id, files ->
            def sample_dir = file("${params.outdir}/per_sample_processing/${id}")
            sample_dir.mkdirs()
            files.each {
                def exon_dir = sample_dir.resolve(it.exon)
                exon_dir.mkdirs()
                it.file.copyTo(exon_dir.resolve(it.file.name))
            }
            return sample_dir
        }
        .collect()

    // Stage the existing per_sample_processing directory
    ch_per_sample_processing = Channel.fromPath("${params.outdir}/per_sample_processing", type: 'dir')

    ABOSNPS2PHENO (
        ch_SNP_reports,
        ch_per_sample_processing
    )
    ch_versions = ch_versions.mix(ABOSNPS2PHENO.out.versions.first())

    emit:
    versions = ch_versions    // channel: [ versions.yml ]
}
