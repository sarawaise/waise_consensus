#!/usr/bin/env nextflow
// Startup: source .startup.sh
// Command: nextflow run main.nf -profile singularity 
// Test: SVABA[13m]; DELLY[1m]; LUMPY[3m]; MANTA[5m]; GRIDSSPL[36m]

nextflow.enable.dsl=2

include { JASMINE_PANEL as MANTA_NORMAL_PANEL } from './modules/jasmine_panel.nf'
include { JASMINE_PANEL as LUMPY_NORMAL_PANEL } from './modules/jasmine_panel.nf'
include { JASMINE_PANEL as SVABA_NORMAL_PANEL } from './modules/jasmine_panel.nf'
include { JASMINE_PANEL as DELLY_NORMAL_PANEL } from './modules/jasmine_panel.nf'
include { JASMINE_PANEL as GRIDSS_NORMAL_PANEL } from './modules/jasmine_panel.nf'

include { TUMOUR_FILTER as MANTA_TUMOUR_FILTER } from './modules/tumour_filter.nf'
include { TUMOUR_FILTER as LUMPY_TUMOUR_FILTER } from './modules/tumour_filter.nf'
include { TUMOUR_FILTER as SVABA_TUMOUR_FILTER } from './modules/tumour_filter.nf'
include { TUMOUR_FILTER as DELLY_TUMOUR_FILTER } from './modules/tumour_filter.nf'
include { TUMOUR_FILTER as GRIDSS_TUMOUR_FILTER } from './modules/tumour_filter.nf'

include { TUMOUR_CONSENSUS_VCF } from './modules/tumour_consensus.nf'
include { TUMOUR_CONSENSUS_CALL } from './modules/tumour_consensus.nf'

process MANTA {
    
    publishDir "${params.publishDir}", mode: 'copy'

    container 'quay.io/biocontainers/manta:1.6.0--py27_0'

    input:
        tuple val(id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
        path reference_fasta
        path reference_index

    output:
        path "manta/${id}/", emit: outdir
        path "manta/${id}/results/variants/${id}*diploidSV.vcf", emit: vcf_n
        tuple val(id), path("manta/${id}/results/variants/${id}*somaticSV.vcf"), emit: vcf_t

    script:
    def manta_outdir = "manta/${id}/"
    """
    configManta.py \
    --tumorBam ${tumor_bam} \
    --normalBam ${normal_bam} \
    --reference ${params.reference_fasta} \
    --runDir ${manta_outdir}

    # Run manta
    ${manta_outdir}/runWorkflow.py

    # Unzip for jasmine compatibility
    gunzip -c ${manta_outdir}/results/variants/somaticSV.vcf.gz > ${manta_outdir}/results/variants/${id}_somaticSV.vcf
    gunzip -c ${manta_outdir}/results/variants/diploidSV.vcf.gz > ${manta_outdir}/results/variants/${id}_diploidSV.vcf
    """

}

process LUMPY {

    publishDir "${params.publishDir}", mode: 'copy'
    
    container 'quay.io/biocontainers/smoove:0.2.8--h9ee0642_1'

    input:
        tuple val(id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)

    output:
        path "lumpy/${id}/", emit: outdir
        path "lumpy/${id}/*_normal.vcf", emit: vcf_n
        tuple val(id), path("lumpy/${id}/*_tumor.vcf"), emit: vcf_t

    script:
    def lumpy_outdir = "lumpy/${id}/"

    """
    smoove call \
        --name ${id} \
        --outdir ${lumpy_outdir} \
        -f ${params.reference_fasta} \
        -processes 1 \
        --removepr \
        --support 3 \
        ${normal_bam} \
        ${tumor_bam}

    # Get sample IDs from VCF header
    echo "Getting sample IDs from joint VCF"
    lumpy_file="${lumpy_outdir}/${id}-smoove.vcf"
    gunzip -k \${lumpy_file}.gz
    sids=`vcfutils.pl listsam \$lumpy_file`
    normal_id=`echo \$sids | cut -f1 -d' '`
    tumor_id=`echo \$sids | cut -f2 -d' '`

    # Split VCF into tumor and normal
    echo "Splitting joint output VCF"
    vcfutils.pl subsam \$lumpy_file \$normal_id > ${lumpy_outdir}/${id}_normal.vcf
    vcfutils.pl subsam \$lumpy_file \$tumor_id > ${lumpy_outdir}/${id}_tumor.vcf
    """

}

process SVABA {

    publishDir "${params.publishDir}", mode: 'copy'

    container 'quay.io/biocontainers/svaba:1.1.0--h468198e_3'

    input:
        tuple val(id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)

    output:
        path "svaba/${id}/", emit: outdir
        path "svaba/${id}/*svaba.germline.sv.vcf", emit: vcf_n
        tuple val(id), path("svaba/${id}/*svaba.unfiltered.somatic.sv.vcf"), emit: vcf_t

    script:
    def svaba_outdir = "svaba/${id}/"
    """
    mkdir -p ${svaba_outdir}
    tbam=\$(pwd)/${tumor_bam}
    nbam=\$(pwd)/${normal_bam}

    cd ${svaba_outdir}
    svaba run \
        -t \$tbam \
        -n \$nbam \
        -a ${id} \
        --threads $params.threads \
        --reference-genome $params.reference_fasta
    """
}

process DELLY {

    publishDir "${params.publishDir}", mode: 'copy'

    container 'quay.io/biocontainers/delly:1.0.3--h358d541_4'

    input:
        tuple val(id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)

    output:
        val "${id}", emit: id
        path "delly/${id}/", emit: outdir
        path "delly/${id}/*_delly.bcf", emit: bcf

    script:
    def delly_outdir = "delly/${id}/"
    
    """
    mkdir -p ${delly_outdir}
    delly call \
    -g $params.reference_fasta \
    -o ${delly_outdir}/${id}_delly.bcf \
    ${tumor_bam} \
    ${normal_bam}
    """
}

process DELLY_SPLIT {

    publishDir "${params.publishDir}", mode: 'copy'

    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'

    input:
        val id
        path delly_bcf

    output:
        path "delly/${id}/*_normal.vcf", emit: vcf_n
        tuple val(id), path("delly/${id}/*_tumor.vcf"), emit: vcf_t

    script:
    def delly_outdir = "delly/${id}/"

    """
    mkdir -p ${delly_outdir} # For publishing
    sids=`bcftools query -l ${delly_bcf}`
    tumor_id=`echo \$sids | cut -f1 -d' '` # ID order is as given to DELLY_CALL
    normal_id=`echo \$sids | cut -f2 -d' '`

    # Split VCF into tumor and normal
    echo "Splitting joint output VCF"
    bcftools view -c1 -Ov -s \$normal_id -o ${delly_outdir}/${id}_normal.vcf ${delly_bcf}
    bcftools view -c1 -Ov -s \$tumor_id -o ${delly_outdir}/${id}_tumor.vcf ${delly_bcf}
    """
}
    
process GRIDSS {

    publishDir "${params.publishDir}", mode: 'copy'

    container 'gridss/gridss:2.13.2'

    input:
        tuple val(id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)

    output:
        val "${id}", emit: id
        path "gridsspl/${id}/", emit: outdir
        path "gridsspl/${id}/${id}_normal.vcf", emit: vcf_n
        tuple val(id), path("gridsspl/${id}/${id}_tumor.vcf"), emit: vcf_t

    script:
    def gridsspl_outdir = "gridsspl/${id}/"
    def vcf_pre = "gridsspl/${id}/${id}"
    def gridss_args = ""
    """
    # Run gridss
    gridss \
            -o ${gridsspl_outdir}/gridss.vcf \
            -r $params.reference_fasta \
            -b $params.gridss_blacklist \
            -c $params.gridss_properties \
            --threads $params.threads \
            --jvmheap $params.jvmheap \
            $gridss_args \
            $normal_bam $tumor_bam

    # Split BAMs
    sids=`bcftools query -l ${gridsspl_outdir}/gridss.vcf`
    normal_id=`echo \$sids | cut -f1 -d' '` # ID order is as given to GRIDSSPL
    tumor_id=`echo \$sids | cut -f2 -d' '`

    # Split VCF into tumor and normal
    echo "Splitting joint output VCF"
    bcftools view -c1 -Ov -s \$normal_id -o ${vcf_pre}_normal.vcf ${gridsspl_outdir}/gridss.vcf
    bcftools view -c1 -Ov -s \$tumor_id -o ${vcf_pre}_tumor.vcf ${gridsspl_outdir}/gridss.vcf
    """
}

workflow {

    // Load inputs from csv
    bam_channel = Channel
            .fromPath(params.csv)
            .splitCsv(header: true, sep: ',', strip: true)
            .map { row -> tuple(
                row.id,
                file(row.tumor_bam),
                file(row.tumor_bai),
                file(row.normal_bam),
                file(row.normal_bai)
                ) } 
    
    // Run callers on BAM
    MANTA(bam_channel, file(params.reference_fasta), file(params.reference_index))
    LUMPY(bam_channel)
    SVABA(bam_channel)
    DELLY(bam_channel)
    DELLY_SPLIT(DELLY.out.id, DELLY.out.bcf)
    GRIDSS(bam_channel)

    // Combine VCFs from normal outputs
    manta_normals = MANTA.out.vcf_n.collect()
    manta_panel = MANTA_NORMAL_PANEL(manta_normals, 'manta')
    lumpy_normals = LUMPY.out.vcf_n.collect()
    lumpy_panel = LUMPY_NORMAL_PANEL(lumpy_normals, 'lumpy')
    svaba_normals = SVABA.out.vcf_n.collect()
    svaba_panel = SVABA_NORMAL_PANEL(svaba_normals, 'svaba')
    delly_normals = DELLY_SPLIT.out.vcf_n.collect()
    delly_panel = DELLY_NORMAL_PANEL(delly_normals, 'delly')
    gridss_normals = GRIDSS.out.vcf_n.collect()
    gridss_panel = GRIDSS_NORMAL_PANEL(gridss_normals, 'gridss')

    // Filter tumor VCFs using panels
    MANTA_TUMOUR_FILTER(MANTA.out.vcf_t, manta_panel, 'manta')
    LUMPY_TUMOUR_FILTER(LUMPY.out.vcf_t, lumpy_panel, 'lumpy')
    SVABA_TUMOUR_FILTER(SVABA.out.vcf_t, svaba_panel, 'svaba')
    DELLY_TUMOUR_FILTER(DELLY_SPLIT.out.vcf_t, delly_panel, 'delly')
    GRIDSS_TUMOUR_FILTER(GRIDSS.out.vcf_t, gridss_panel, 'gridss')
    
    // Combine all tumor channels by Sample ID
    TUMOR_COMBINED = MANTA_TUMOUR_FILTER.out.filtered_vcf
        .combine(SVABA_TUMOUR_FILTER.out.filtered_vcf, by: 0)
        .combine(LUMPY_TUMOUR_FILTER.out.filtered_vcf, by: 0)
        .combine(DELLY_TUMOUR_FILTER.out.filtered_vcf, by: 0)
        .combine(GRIDSS_TUMOUR_FILTER.out.filtered_vcf, by: 0)
    
    // Aggregate VCF files for each sample
    //  Create channel with (val(id), [path(vcf), path(vcf), ...])
    TUMOR_FORMATTED = TUMOR_COMBINED.map { v -> tuple(v[0], v.subList(1, v.size())) } // `*vcfs` captures all elements except the first
    //  Create single VCF file for each sample
    TUMOUR_CONSENSUS_VCF(TUMOR_FORMATTED)

    // Call consensus SVs
    TUMOUR_CONSENSUS_CALL(TUMOUR_CONSENSUS_VCF.out.id, TUMOUR_CONSENSUS_VCF.out.consensus_vcf)

}
