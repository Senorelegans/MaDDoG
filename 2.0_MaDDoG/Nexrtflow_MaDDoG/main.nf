#!/usr/bin/env nextflow
/*


*/

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.name = false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.plaintext_email = false
params.bbmap_adapters = "$baseDir/bin/adapters.fa"
params.bedGraphToBigWig = "$baseDir/bin/bedGraphToBigWig"
params.rcc = "$baseDir/bin/rcc.py"
params.maddog = "$baseDir/bin/nb_MaDDoG.py"
params.rsubreadpy = "$baseDir/bin/40_Create_R_subread_DESeq2_NEXTFLOW.py"
params.rsubreadRtemplate = "$baseDir/bin/R_subread_DEseq2_TEMPLATE.R"
multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")



// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}



process MADDOG {
    memory '30 GB'
    validExitStatus 0
    tag "$chr"
    cpus 4
    maxForks 24
    time '10h'
    publishDir "${params.outdir}/MaDDoG/gtf/${chr}", mode: 'copy', pattern: "*gtf*"        
    publishDir "${params.outdir}/MaDDoG/Chrom_out/${chr}", mode: 'copy', pattern: "*"


    input:
    set val(chr) from file(params.CHR_NAMES).splitCsv(header: true, sep:'\t')

    output:
    file("*ALL_BestModel_*") into bestmodel
    file("*") into maddog_out  
    set val(chr), file("*gtf") into maddog_gtf       

    script:
    """

    python ${params.maddog} \
            --chrom ${chr} \
            --DOGOUT_PATH ${params.outdir} \
            --MADDOG_ONLY True \
            --max_num_states ${params.max_num_states} \
            --bayes_factor_thresh ${params.bayes_factor_thresh} \
            --iterations ${params.iterations} \
            --convolution_window ${params.convolution_window}
    """
 }


process RSUBREAD_DESEQ {
    memory '50 GB'
    validExitStatus 0
    cpus 18
    maxForks 12
    time '14h'
    publishDir "${params.outdir}/MaDDoG/rsubread_deseq/rsubread_deseq_${chr}", mode: 'copy', pattern: "*"

    input:
    set val(chr), file(gtf) from maddog_gtf

    output:
    file("*") into rsubread_out

    script:
    """
    python ${params.rsubreadpy} \
    --annotation_file_with_path ${gtf} \
    --control_BAM_list \
    ${params.BAMPATH}/CTL1.bam \
    ${params.BAMPATH}/CTL2.bam \
    ${params.BAMPATH}/CTL3.bam \
    --treatment_BAM_list \
    ${params.BAMPATH}/KCL1.bam \
    ${params.BAMPATH}/KCL2.bam \
    ${params.BAMPATH}/KCL3.bam \
    --input_R_template_file ${params.rsubreadRtemplate} \
    --cpus 18 \
    --padj 0.05

    R CMD BATCH Rsubread_DESeq2_initial.R Rsubread_DESeq2_initial.R.out
    """
 }







// process Flatten {
//     memory '50 GB'
//     publishDir "${params.outdir}/gtf", mode: 'copy', pattern: "*"

//     input:
//     file(ANNOTATION) from annotation_in

//     output:
//     file("*CHROMNAMES.txt") into flat_gtf_chrnames
//     file("*_flat.gtf") into flat_gtf

//     script:
//     """

//     python ${params.flatten} \
//               --annotation_in ${annotation_in} \
//               --file_type ${params.annotation_type}
//     """
//  }





// process trim {
//     validExitStatus 0
//     tag "$name"
//     cpus 20
//     maxForks 6
//     time '24h'
//     memory '20 GB'
//     //publishDir "${params.outdir}/qc/trimstats", mode: 'copy', pattern: "*stats.txt"
//     publishDir "${params.outdir}/qc/fastq_count", mode: 'copy', pattern: "*fastq_count.txt"
//     publishDir "${params.outdir}/qc/fastqtrim_count", mode: 'copy', pattern: "*fastqtrim_count.txt"    
//     if (params.saveTrim || params.saveAllfq) {
//         publishDir "${params.outdir}/fastq_trimmed", mode: 'copy', pattern: "*.fastq.gz"
//     }      

//     input:
//     set val(name), val(condition) from fastq_reads_trim.splitCsv(header: ["name", "condition"],sep:'\t')

 



 