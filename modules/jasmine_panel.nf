// Process VCF files from Normal samples with Jasmine. 
// For each tool, the subworkflow JASMINE_NORMAL_PANEL create a VCF with SVs observed in multiple samples.


process MAKE {
    // Merge VCF files into a single VCF for downstream processing

    publishDir "${params.publishDir}/panel", mode: 'copy'

    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'

    input:
        path vcf
        val tool

    output:
        path "${tool}_jasmine_cons.vcf", emit: panel_vcf

    script:

    """
    awk '/^##INFO=<tool=SUPP/ { sub("Type=String", "Type=Integer"); } { print }' ${vcf} > tmp.vcf0
    bcftools filter -i "INFO/SUPP > 1" tmp.vcf0 > ${tool}_jasmine_cons.vcf
    """
}

process MERGE {
    // Merge VCF files into a single VCF for downstream processing

    publishDir "${params.publishDir}/jasmine", mode: 'copy'

    container 'quay.io/biocontainers/jasminesv:1.1.5--hdfd78af_0'

    input:
        path vcf_files
        val tool

    output:
        path "${tool}_jasmine_out.vcf", emit: merged_vcf
        path "${tool}_jasmine_filelist.txt", emit: filelist

    script:
    """
    # Create a file list for Jasmine
    ls -1 ${vcf_files} > ${tool}_jasmine_filelist.txt

    # Run Jasmine
    jasmine file_list=${tool}_jasmine_filelist.txt out_file=${tool}_jasmine_out.vcf
    """
}

workflow JASMINE_PANEL {
    

    take:
    vcf_files 
    tool // identifier

    main:
    (merged_vcf, file_list) = MERGE(vcf_files, tool)
    panel_vcf  = MAKE(merged_vcf, tool)

    emit:
    panel_vcf  // channel: vcf
}