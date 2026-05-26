/*
  SUBWORKFLOW: PREDICTABOPHENOTYPE
*/

//  Description: Predicts ABO phenotype by combining variant frequency data
//   with BAM coverage, extracting SNPs, and mapping them to phenotypes.
//   Uses metadata-driven joining and structured per-sample output.


include { ABO_GETABOSNPS } from '../../../modules/local/abo/getabosnps/main'
include { ABO_SNPS2PHENO } from '../../../modules/local/abo/snps2pheno/main'

workflow PREDICTABOPHENOTYPE {
    take:
    ch_variants_freq // channel: [ val(meta), [ freq ] ] - with exon metadata
    ch_bam_coverage // channel: [ val(meta), [ cov ] ] - with exon metadata

    main:

    // JOIN: Variant frequency with BAM coverage
    ch_combined_input = ch_variants_freq
        .join(ch_bam_coverage)
        .map { meta, freq, cov ->
            [meta, freq, cov, meta.exon]
        }

    /*
    MODULE: ABO_GETABOSNPS
    */
    ABO_GETABOSNPS(
        ch_combined_input
    )

    // PREP:Organize SNP reports by sample and exon
    ch_snp_reports = ABO_GETABOSNPS.out.phenotype
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

    // STAGE: Existing per_sample_processing directory in results
    ch_per_sample_processing = channel.fromPath("${params.outdir}/per_sample_processing", type: 'dir')

    /*
    MODULE: ABO_SNPS2PHENO
    */
    ABO_SNPS2PHENO(
        ch_snp_reports,
        ch_per_sample_processing,
    )

    emit:
    abo_results = ABO_SNPS2PHENO.out.txt
}
