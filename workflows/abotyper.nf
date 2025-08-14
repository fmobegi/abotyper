/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
// include { MULTIQC                } from '../modules/nf-core/multiqc/main'
// Using  MULTIQC v1.28 (broad compatibility)
// This is a temporary workaround until we can figure out the AVX2 compatibility issues with MultiQC v1.30
// Some python libraries are being compiled with AVX2 instructions under containerised environments, even though local versions behave differently.
include { MULTIQC                } from '../modules/local/multiqc/main'
include { MAKEINDEX              } from '../modules/local/makeindex/main'

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_abotyper_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MINIMAP2_ALIGN_READS        } from '../subworkflows/local/minimap_align_exons/main'
include { PREDICTABOPHENOTYPE         } from '../subworkflows/local/predictabophenotype/main'
include { VARIANTS_QUANTIFICATION     } from '../subworkflows/local/variant_calling_mpileup/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ABOTYPER {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    fai            // channel: fasta from params.fai (ABO database)
    fasta          // channel: fasta from params.fasta (ABO database)
    exon6fai       // channel: fasta from params.exon6fai
    exon6fasta     // channel: fasta from params.exon6fasta
    exon7fai       // channel: fasta from params.exon7fai
    exon7fasta     // channel: fasta from params.exon7fasta
    logo           // channel: png from params.logo (custom pathwest logo)

    main:

    ch_versions = Channel.empty()

    //
    // Prepare sample channels with exon metadata for mapping to each reference
    //
    ch_exon6_samples = ch_samplesheet
        .map { meta, fastq ->
            def new_meta = meta + [exon: 'exon6']
            [new_meta, fastq]
        }

    ch_exon7_samples = ch_samplesheet
        .map { meta, fastq ->
            def new_meta = meta + [exon: 'exon7']
            [new_meta, fastq]
        }

    // // DEBUG: Check individual channels
    // ch_exon6_samples.view { "EXON6 SAMPLE: $it" }
    // ch_exon7_samples.view { "EXON7 SAMPLE: $it" }

    // Combine sample channels
    ch_combined_input = ch_exon6_samples.mix(ch_exon7_samples)

    // // DEBUG: Check combined channel
    // ch_combined_input.view { "COMBINED INPUT: $it" }

    // Use reference channels as-is (they already have exon metadata)
    ch_combined_fasta = exon6fasta.mix(exon7fasta)
    ch_combined_fai = exon6fai.mix(exon7fai)

    //
    // MODULE: Create index files once for reference sequences
    //
    MAKEINDEX (
        exon6fai,
        exon7fai
    )
    ch_versions = ch_versions.mix(MAKEINDEX.out.versions)

    //
    // MODULE: FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // SUBWORKFLOW: minimap2/align (with pre-prepared inputs for better caching)
    //
    MINIMAP2_ALIGN_READS (
        ch_combined_input,
        ch_combined_fasta,
        ch_combined_fai
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN_READS.out.versions)

    //
    // SUBWORKFLOW: Run pileup for variants with properly combined bed files
    //
    ch_combined_bed = MAKEINDEX.out.exon6bed
        .map { meta, bed -> [meta + [exon: 'exon6'], bed] }
        .mix(
            MAKEINDEX.out.exon7bed.map { meta, bed -> [meta + [exon: 'exon7'], bed] }
        )

    VARIANTS_QUANTIFICATION(
        MINIMAP2_ALIGN_READS.out.bam,
        MINIMAP2_ALIGN_READS.out.bai,
        MINIMAP2_ALIGN_READS.out.fasta,
        MINIMAP2_ALIGN_READS.out.fai,
        ch_combined_bed  // Pass combined bed files with exon metadata
    )
    ch_versions = ch_versions.mix(VARIANTS_QUANTIFICATION.out.versions)

    //
    // SUBWORKFLOW: Run ABO prediction
    //
    PREDICTABOPHENOTYPE(
        VARIANTS_QUANTIFICATION.out.metrics,
        MINIMAP2_ALIGN_READS.out.coverage
    )
    ch_versions = ch_versions.mix(PREDICTABOPHENOTYPE.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'abotyper_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // Prepare all MultiQC inputs without staging
    //
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

    // Collect all MultiQC inputs directly without staging
    ch_multiqc_inputs = Channel.empty()
        .mix(FASTQC.out.zip.map { meta, files -> files }.flatten())
        .mix(MINIMAP2_ALIGN_READS.out.coverage.map { meta, cov -> cov })
        .mix(VARIANTS_QUANTIFICATION.out.metrics.map { meta, metrics -> metrics })
        .mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        .mix(ch_collated_versions)
        .mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    MULTIQC (
        ch_multiqc_inputs.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
