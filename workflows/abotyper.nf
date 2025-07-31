/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
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
    ch_multiqc_files = Channel.empty()
    
    //
    // Prepare input channels with exon metadata ONCE (for better caching)
    //
    ch_exon6_input = ch_samplesheet
        .map { meta, fastq -> 
            def new_meta = meta + [exon: 'exon6']
            [new_meta, fastq]
        }
    
    ch_exon7_input = ch_samplesheet
        .map { meta, fastq -> 
            def new_meta = meta + [exon: 'exon7']
            [new_meta, fastq]
        }

    ch_combined_input = ch_exon6_input.mix(ch_exon7_input)

    // Prepare reference files with matching metadata ONCE
    ch_exon6_fasta_prepared = exon6fasta
        .map { ref_meta, fasta -> 
            [ref_meta + [exon: 'exon6'], fasta]
        }
    
    ch_exon7_fasta_prepared = exon7fasta
        .map { ref_meta, fasta -> 
            [ref_meta + [exon: 'exon7'], fasta]
        }

    ch_combined_fasta = ch_exon6_fasta_prepared.mix(ch_exon7_fasta_prepared)

    ch_exon6_fai_prepared = exon6fai
        .map { ref_meta, fai -> 
            [ref_meta + [exon: 'exon6'], fai]
        }
    
    ch_exon7_fai_prepared = exon7fai
        .map { ref_meta, fai -> 
            [ref_meta + [exon: 'exon7'], fai]
        }

    ch_combined_fai = ch_exon6_fai_prepared.mix(ch_exon7_fai_prepared)
    
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

    // Collect fastqc reports for multiqc
    ch_multiqc_files = ch_multiqc_files.mix(
        FASTQC.out.zip.map { meta, zip -> 
            def new_name = "${meta.id}_fastqc.zip"
            [zip, new_name]
        }
    )
    
    //
    // SUBWORKFLOW: minimap2/align (with pre-prepared inputs for better caching)
    //
    MINIMAP2_ALIGN_READS (
        ch_combined_input,
        ch_combined_fasta,
        ch_combined_fai
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN_READS.out.versions)
    
    // Collect alignment QC files for MultiQC
    ch_multiqc_files = ch_multiqc_files.mix(
        MINIMAP2_ALIGN_READS.out.coverage.map { meta, cov -> 
            def new_name = "${meta.id}_${meta.exon}.coverage.txt"
            [cov, new_name]
        }
    )
    
    //
    // SUBWORKFLOW: Run pileup for variants
    //
    VARIANTS_QUANTIFICATION(
        MINIMAP2_ALIGN_READS.out.bam,
        MINIMAP2_ALIGN_READS.out.bai,
        MINIMAP2_ALIGN_READS.out.fasta,
        MINIMAP2_ALIGN_READS.out.fai,
        MAKEINDEX.out.exon6bed,  // Pass pre-computed exon 6 bed files
        MAKEINDEX.out.exon7bed   // Pass pre-computed exon 7 bed files
    )
    ch_versions = ch_versions.mix(VARIANTS_QUANTIFICATION.out.versions)
    
    // Collect variant metrics for MultiQC
    ch_multiqc_files = ch_multiqc_files.mix(
        VARIANTS_QUANTIFICATION.out.metrics.map { meta, metrics -> 
            def new_name = "${meta.id}_${meta.exon}.freq.tsv"
            [metrics, new_name]
        }
    )
    
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

    // Stage all files for MultiQC with proper naming
    ch_staged_files = ch_multiqc_files
        .collectFile() { file, new_name ->
            [new_name, file]
        }
        .collect()

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

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    
    ch_multiqc_files_final = ch_staged_files
        .mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        .mix(ch_collated_versions)

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files_final = ch_multiqc_files_final.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files_final.collect(),
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