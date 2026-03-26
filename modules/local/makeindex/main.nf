process MAKEINDEX {
    tag "FAI to BED"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(exon6fai)
    tuple val(meta1), path(exon7fai)

    output:
    tuple val(meta), path("*_exon6.bed"), emit: exon6bed
    tuple val(meta1), path("*_exon7.bed"), emit: exon7bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id.toString().tokenize('.')[0]}"
    def prefix1 = task.ext.prefix1 ?: "${meta1.id.toString().tokenize('.')[0]}"

    """
    awk -v FS="\\t" -v OFS="\\t" '{print \$1 FS "0" FS (\$2)-1}' ${exon6fai} > ${prefix}_exon6.bed
    awk -v FS="\\t" -v OFS="\\t" '{print \$1 FS "0" FS (\$2)-1}' ${exon7fai} > ${prefix1}_exon7.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1 | sed 's/.*Awk //; s/,.*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id.toString().tokenize('.')[0]}"
    def prefix1 = task.ext.prefix1 ?: "${meta1.id.toString().tokenize('.')[0]}"

    """
    touch ${prefix}_exon6.bed
    touch ${prefix1}_exon7.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1 | sed 's/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
