/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                  } from '../modules/nf-core/fastqc'
include { MULTIQC                 } from '../modules/nf-core/multiqc'
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText  } from '../subworkflows/local/utils_nfcore_abotyper_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MAKEINDEX               } from '../modules/local/makeindex'
include { MINIMAP2_ALIGN_READS    } from '../subworkflows/local/minimap_align_exons'
include { PREDICTABOPHENOTYPE     } from '../subworkflows/local/predictabophenotype'
// include { VARIANTS_QUANTIFICATION } from '../subworkflows/local/variant_calling_mpileup'
include { VARIANTS_QUANTIFICATION } from '../subworkflows/local/variant_calling_haploscan'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ABOTYPER {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    exon6fai // channel: fasta from params.exon6fai
    exon6fasta // channel: fasta from params.exon6fasta
    exon7fai // channel: fasta from params.exon7fai
    exon7fasta // channel: fasta from params.exon7fasta
    logo // channel: png from params.logo (custom pathwest logo)

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    // Prepare sample channels with exon metadata for mapping to each reference
    ch_exon6_samples = ch_samplesheet.map { meta, fastq ->
        def new_meta = meta + [exon: 'exon6']
        [new_meta, fastq]
    }

    ch_exon7_samples = ch_samplesheet.map { meta, fastq ->
        def new_meta = meta + [exon: 'exon7']
        [new_meta, fastq]
    }

    ch_combined_input = ch_exon6_samples.mix(ch_exon7_samples)
    ch_combined_fasta = exon6fasta.mix(exon7fasta)
    ch_combined_fai = exon6fai.mix(exon7fai)

    /*
    MODULE: MAKEINDEX
    */
    MAKEINDEX(
        exon6fai,
        exon7fai,
    )

    /*
    MODULE: FASTQC
    */
    FASTQC(
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})

    /*
    SUBWORKFLOW: MINIMAP2_ALIGN_READS
    */
    MINIMAP2_ALIGN_READS(
        ch_combined_input,
        ch_combined_fasta,
        ch_combined_fai,
    )

    /*
    SUBWORKFLOW: VARIANTS_QUANTIFICATION
    */
    VARIANTS_QUANTIFICATION(
        MINIMAP2_ALIGN_READS.out.bam,
        MINIMAP2_ALIGN_READS.out.bai,
        MINIMAP2_ALIGN_READS.out.fasta,
        MINIMAP2_ALIGN_READS.out.fai,
        MAKEINDEX.out.exon6bed.mix(MAKEINDEX.out.exon7bed),
    )

    /*
    SUBWORKFLOW: PREDICTABOPHENOTYPE
    */
    // Join metrics and coverage by metadata to ensure correct pairing
    ch_prediction_input = VARIANTS_QUANTIFICATION.out.metrics
        .join(MINIMAP2_ALIGN_READS.out.coverage)
        .multiMap { meta, metrics, coverage ->
            metrics: [meta, metrics]
            coverage: [meta, coverage]
        }

    PREDICTABOPHENOTYPE(
        ch_prediction_input.metrics,
        ch_prediction_input.coverage,
        VARIANTS_QUANTIFICATION.out.haplotypes,
    )

    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'abotyper_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)

    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))

    def ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'abotyper'],
                files,
                params.multiqc_config
                    ? file(params.multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                params.multiqc_logo ? file(params.multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList()
    versions       = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
