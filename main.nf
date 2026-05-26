#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/abotyper
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/abotyper
    Website: https://nf-co.re/abotyper
    Slack  : https://nfcore.slack.com/channels/abotyper
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ABOTYPER                } from './workflows/abotyper'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_abotyper_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_abotyper_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_abotyper_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.exon6fai   = getGenomeAttribute('exon6fai')
params.exon6fasta = getGenomeAttribute('exon6fasta')
params.exon7fai   = getGenomeAttribute('exon7fai')
params.exon7fasta = getGenomeAttribute('exon7fasta')
params.logo       = getGenomeAttribute('logo')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NFCORE_ABOTYPER {
    take:
    samplesheet // channel: samplesheet read in from --input

    main:
    exon6fai = params.exon6fai ? Channel.fromPath(params.exon6fai).map { it -> [[id: it.baseName, exon: 'exon6'], it] } : Channel.empty()
    exon6fasta = params.exon6fasta ? Channel.fromPath(params.exon6fasta).map { it -> [[id: it.baseName, exon: 'exon6'], it] } : Channel.empty()
    exon7fai = params.exon7fai ? Channel.fromPath(params.exon7fai).map { it -> [[id: it.baseName, exon: 'exon7'], it] } : Channel.empty()
    exon7fasta = params.exon7fasta ? Channel.fromPath(params.exon7fasta).map { it -> [[id: it.baseName, exon: 'exon7'], it] } : Channel.empty()
    logo = params.logo ? Channel.fromPath(params.logo).collect() : Channel.empty()

    ABOTYPER(
        samplesheet,
        exon6fai,
        exon6fasta,
        exon7fai,
        exon7fasta,
        logo
    )

    emit:
    multiqc_report = ABOTYPER.out.multiqc_report
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_ABOTYPER(
        PIPELINE_INITIALISATION.out.samplesheet
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_ABOTYPER.out.multiqc_report
    )
}
