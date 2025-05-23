
process TUMOUR_FILTER {
    // Merge VCF files into a single VCF for downstream processing

    publishDir "${params.publishDir}/filtered", mode: 'copy'

    container 'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2'

    input:
        tuple val(id), path(vcf_file)
        path normal_panel
        val tool

    output:
        tuple val(id), path("${id}_${tool}_filtered.vcf"), emit: filtered_vcf

    script:
    """
    # Add header
    grep '^#' ${vcf_file} > ${id}_${tool}_filtered.vcf
    bedtools subtract -a ${vcf_file} -b ${normal_panel} >> ${id}_${tool}_filtered.vcf
    """
}
