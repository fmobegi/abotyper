/*
Typically, in variant calling, a "variant" is defined as a position where the observed sequence
differs from the reference genome.
When REF and ALT are the same, it's not usually considered a variant in the traditional sense.
However, for ABO analysis, it is necessary to include all relevant REF positions in the decision-making tree.
The samtools/mpileup module output is processed using python3 to achieve this.
*/

process MPILEUP_NUCL_FREQ {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.9--1'
        : 'quay.io/biocontainers/python:3.9--1'}"

    input:
    tuple val(meta), path(pileup)
    tuple val(meta1), path(fasta)

    output:
    tuple val(meta), path("*.AlignmentStatistics.tsv"), emit: tsv
    path ("ABOReadPolymorphisms.txt"), emit: txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    stats_from_pileup.py \\
        -i ${pileup} \\
        -o ${prefix}.AlignmentStatistics.tsv \\
        -s ABOReadPolymorphisms.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        re: \$(python3 -c "import re; print(re.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.AlignmentStatistics.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        re: \$(python3 -c "import re; print(re.__version__)")
    END_VERSIONS
    """
}
