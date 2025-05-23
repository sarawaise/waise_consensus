
process FILTER_TUMOUR {
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

process TUMOUR_CONSENSUS_VCF {
    // Merge VCF files into a single VCF for downstream processing

    publishDir "${params.publishDir}/somatic", mode: 'copy'

    container 'quay.io/biocontainers/jasminesv:1.1.5--hdfd78af_0'

    input:
        tuple val(id), path(vcf_files)

    output:
        val id, emit: id
        path "${id}_consensus_jasmine_out.vcf", emit: consensus_vcf
        path "consensus_filelist.txt", emit: filelist

    script:
    """
    # Create a file list for Jasmine
    ls -1 ${vcf_files} > consensus_filelist.txt

    # Run Jasmine
    jasmine file_list=consensus_filelist.txt out_file=${id}_consensus_jasmine_out.vcf
    """
}

process TUMOUR_CONSENSUS_CALL {
    // Merge VCF files into a single VCF for downstream processing

    publishDir "${params.publishDir}/somatic", mode: 'copy'

    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'

    input:
        val id
        path vcf

    output:
        path "${id}_tumor_cons.vcf", emit: panel_vcf
        val id, emit: id

    script:

    """
    awk '/^##INFO=<ID=SUPP/ { sub("Type=String", "Type=Integer"); } { print }' ${vcf} > tmp.vcf0
    bcftools filter -i "INFO/SUPP > 4" tmp.vcf0 > ${id}_tumor_cons.vcf
    """
}

