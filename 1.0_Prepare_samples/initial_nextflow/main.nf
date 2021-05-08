#!/usr/bin/env nextflow
/*
========================================================================================
                         SteadyFlow - Steady Sate Transcription PIPELINE
========================================================================================
Steady Sate Analysis Pipeline. Started 2018-06-21.
 #### Homepage / Documentation
 https://github.com/Dowell-Lab/RNA-seq-Flow
 #### Authors
 Margaret Gruca <magr0763@colorado.edu>
========================================================================================
========================================================================================

Pipeline steps:

    1. Pre-processing sra/fastq
        1a. SRA tools -- fasterq-dump sra to generate fastq file
        1b. FastQC (pre-trim) -- perform pre-trim FastQC on fastq files

    2. Trimming
        2a. BBDuk -- trim fastq files for quality and adapters
        2b. FastQC (post-trim) -- perform post-trim FastQC on fastq files (ensure trimming performs as expected)

    3. Mapping w/ HISAT2 -- map to genome reference file

    4. SAMtools -- convert SAM file to BAM, index BAM, flagstat BAM

    5. Quality control
        5a. preseq -- estimate library complexity
        5b. RSeQC -- calculate genomic coverage relative to a reference file, infer experiement (single- v. paired-end), read duplication
        5c. Pileup.sh : BBMap Suite -- genomic coverage by chromosome, GC content, pos/min reads, intron/exon ratio

    6. Coverage files
        6d. BEDTools : non-normalized & nornmalized bedgraphs
        6b. BEDTools and kentUtils : 5' bigwigs for dREG & normalized bigwigs for genome browser
        
    7. Normalizing bigwigs for Genome Browser use

    8. IGV Tools : bedGraph --> tdf
    
    9. RSeQC Counts : Read counts for a provided RefSeq annotation file

    10. MultiQC : generate QC report for pipeline

*/


def helpMessage() {
    log.info"""
    =========================================
     SteadyFlow v${params.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf -profile slurm --fastqs '/project/*_{R1,R2}*.fastq' --outdir '/project/'

    Required arguments:
         -profile                      Configuration profile to use. <base, slurm>
         --fastqs                      Directory pattern for fastq files: /project/*{R1,R2}*.fastq (Required if --sras not specified)
         --sras                        Directory pattern for SRA files: /project/*.sras (Required if --fastqs not specified)
         --workdir                     Nextflow working directory where all intermediate files are saved.
         --email                       Where to send workflow report email.

    Performance options:
        --threadfqdump                 Runs multi-threading for fastq-dump for sra processing.

    Input File options:
        --singleEnd                    Specifies that the input files are not paired reads (default is paired-end).
        --flip                         Reverse complements each strand. Necessary for some library preps.
        --flipR2                       Reverse complements R2 only.       

    Save options:
        --outdir                       Specifies where to save the output from the nextflow run.
        --savefq                       Compresses and saves raw fastq reads.
        --saveTrim                     Compresses and saves trimmed fastq reads.
        --skipBAM                      Skip saving BAM files. Only CRAM files will be saved with this option.
        --saveAll                      Compresses and saves all fastq reads.

    QC Options:
        --skipMultiQC                  Skip running MultiQC.
        --skipRSeQC                    Skip running RSeQC.
        
    Analysis Options:
        --count                        Run RSeQC FPKM count over RefSeq annotated genes.

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
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
genome_refseq = params.genome_refseq
WINDOW = params.WINDOW
multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")

// Validate inputs

if ( params.genome ){
    genome = file(params.genome)
    if( !genome.exists() ) exit 1, "Genome directory not found: ${params.genome}"
}

if ( params.chrom_sizes ){
    chrom_sizes = file(params.chrom_sizes)
    if( !chrom_sizes.exists() ) exit 1, "Genome chrom sizes file not found: ${params.chrom_sizes}"
 }

if ( params.hisat2_indices ){
    hisat2_indices = file("${params.hisat2_indices}")
}

if ( params.genome_refseq ){
    genome_refseq = file("${params.genome_refseq}")
}

if ( params.bbmap_adapters){
    bbmap_adapters = file("${params.bbmap_adapters}")
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


if (params.forward_stranded && !params.unStranded){
    rnastrandness = params.singleEnd ? '--rna-strandness F' : '--rna-strandness FR'
} else if (params.reverse_stranded && !params.unStranded){
    rnastrandness = params.singleEnd ? '--rna-strandness R' : '--rna-strandness RF'
}


/*
 * Create a channel for input read files
 */
if (params.fastqs) {
    if (params.singleEnd) {
        fastq_reads_qc = Channel
                            .fromPath(params.fastqs)
                            .map { file -> tuple(file.simpleName, file) }
        fastq_reads_trim = Channel
                            .fromPath(params.fastqs)
                            .map { file -> tuple(file.simpleName, file) }
    } else {
        Channel
            .fromFilePairs( params.fastqs, size: params.singleEnd ? 1 : 2 )
            .ifEmpty { exit 1, "Cannot find any reads matching\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
            .into { fastq_reads_qc; fastq_reads_trim; fastq_reads_gzip }
    }
}

else {
    Channel
        .empty()
        .into { fastq_reads_qc; fastq_reads_trim; fastq_reads_gzip }
}

if (params.sras) {
    println("pattern for SRAs provided")
    read_files_sra = Channel
                        .fromPath(params.sras)
                        .map { file -> tuple(file.baseName, file) }
}

else {
    read_files_sra = Channel.empty()
}


// Header log info
log.info """=======================================================
NascentFlow v${params.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']    = 'NascentFlow'
summary['Help Message']     = params.help
summary['Pipeline Version'] = params.version
summary['Run Name']         = custom_runName ?: workflow.runName
if(params.reads) summary['Reads']     = params.reads
if(params.fastqs) summary['Fastqs']   = params.fastqs
if(params.sras) summary['SRAs']       = params.sras
summary['Genome Ref']       = params.genome
summary['Thread fqdump']    = params.threadfqdump ? 'YES' : 'NO'
summary['Data Type']        = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Save All fastq']   = params.saveAllfq ? 'YES' : 'NO'
summary['Save BAM']         = params.skipBAM ? 'NO' : 'YES'
summary['Save fastq']       = params.savefq ? 'YES' : 'NO'
summary['Save Trimmed']     = params.saveTrim ? 'YES' : 'NO'
summary['Reverse Comp']     = params.flip ? 'YES' : 'NO'
summary['Reverse Comp R2']  = params.flipR2 ? 'YES' : 'NO'
summary['Run RSeQC']        = params.skipRSeQC ? 'NO' : 'YES'
summary['Run Count']        = params.count ? 'NO' : 'YES'
summary['Run MultiQC']      = params.skipMultiQC ? 'NO' : 'YES'
summary['Max Memory']       = params.max_memory
summary['Max CPUs']         = params.max_cpus
summary['Max Time']         = params.max_time
summary['Output dir']       = params.outdir
summary['Working dir']      = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']     = "$HOME"
summary['Current user']     = "$USER"
summary['Current path']     = "$PWD"
summary['Output dir']       = params.outdir
summary['Script dir']       = workflow.projectDir
summary['Config Profile']   = workflow.profile
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "======================================================="

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
try {
    if( ! nextflow.version.matches(">= $params.nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}


/*
 * Parse software version numbers
 */
process get_software_versions {
    validExitStatus 0,1,127
    publishDir "${params.outdir}/software_versions/", mode: 'copy', pattern: '*.txt'

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml
    file '*.txt' into software_versions_text

    script:
    """
    echo $params.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    bbversion.sh --version > v_bbduk.txt
    hisat2 --version > v_hisat2.txt
    samtools --version > v_samtools.txt
    fastq-dump --version > v_fastq-dump.txt
    preseq > v_preseq.txt
    bedtools --version > v_bedtools.txt
    igvtools version > v_igv-tools.txt

    for X in `ls *.txt`; do
        cat \$X >> all_versions.txt;
    done
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}

/*
 * Step 1a -- get fastq files from downloaded sras
 */

process sra_dump {
    tag "$prefix"
    if (params.threadfqdump) {
        cpus 8 }
    else {
        cpus 1
    }
    if (params.savefq || params.saveAllfq) {
        publishDir "${params.outdir}/fastq", mode: 'copy'
    }
    
    input:
    set val(prefix), file(reads) from read_files_sra

    output:
    set val(prefix), file("*.fastq.gz") into fastq_reads_qc_sra, fastq_reads_trim_sra, fastq_reads_gzip_sra
   

    script:
    prefix = reads.baseName
    if (!params.threadfqdump) {
        """
        echo ${prefix}

        fastq-dump ${reads} --gzip
        """
    } else if (!params.singleEnd) {
         """
        export PATH=~/.local/bin:$PATH

        parallel-fastq-dump \
            --threads 8 \
            --gzip \
            --split-3 \
            --sra-id ${reads}
        """
    } else if (!params.threadfqdump && !params.singleEnd) {
        """
        echo ${prefix}

        fastq-dump --split-3 ${reads} --gzip
        """
    } else {
        """
        export PATH=~/.local/bin:$PATH

        parallel-fastq-dump \
            --threads 8 \
            --gzip \
            --sra-id ${reads}
        """
    }
}

/*
 * STEP 1b - FastQC
 */

process fastQC {
    validExitStatus 0,1
    tag "$prefix"
    memory '8 GB'
    publishDir "${params.outdir}/qc/fastqc/", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(prefix), file(reads) from fastq_reads_qc.mix(fastq_reads_qc_sra)

    output:
    file "*.{zip,html,txt}" into fastqc_results

    script:
    """
    fastqc $reads
    """
}


/*
 * STEP 2a - Trimming
 */

process bbduk {
    validExitStatus 0,1
    tag "$name"
    cpus 32
    time '4h'
    memory '40 GB'
    publishDir "${params.outdir}/qc/trimstats", mode: 'copy', pattern: "*.txt"
    if (params.saveTrim || params.saveAllfq) {
        publishDir "${params.outdir}/fastq_trimmed", mode: 'copy', pattern: "*.fastq.gz"
    }      

    input:
    set val(name), file(reads) from fastq_reads_trim.mix(fastq_reads_trim_sra)

    output:
    set val(name), file ("*.trim.fastq.gz") into trimmed_reads_fastqc, trimmed_reads_hisat2
    file "*.txt" into trim_stats

    script:
    if (!params.singleEnd && params.flip) {
        """
        echo ${name}

        reformat.sh -Xmx20g \
                t=16 \
                in=${name}_R1.fastq.gz \
                in2=${name}_R2.fastq.gz \
                out=${name}_R1.flip.fastq.gz \
                out2=${name}_R2.flip.fastq.gz \
                rcomp=t

        bbduk.sh -Xmx20g \
                t=16 \
                in=${name}_R1.flip.fastq.gz \
                in2=${name}_R2.flip.fastq.gz \
                out=${name}_R1.flip.trim.fastq.gz \
                out2=${name}_R2.flip.trim.fastq.gz \
                ref=${bbmap_adapters} \
                ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
                nullifybrokenquality=t \
                maq=10 minlen=25 \
                tpe tbo \
                literal=AAAAAAAAAAAAAAAAAAAAAAA \
                stats=${name}.trimstats.txt \
                refstats=${name}.refstats.txt
        """
    } else if (!params.singleEnd && params.flipR2) {
                """
        echo ${name}

        reformat.sh -Xmx20g \
                t=16 \
                in=${name}_R1.fastq.gz \
                in2=${name}_R2.fastq.gz \
                out=${name}_R1.flip.fastq.gz \
                out2=${name}_R2.flip.fastq.gz \
                rcompmate=t

        bbduk.sh -Xmx20g \
                t=16 \
                in=${name}_R1.flip.fastq.gz \
                in2=${name}_R2.flip.fastq.gz \
                out=${name}_R1.flip.trim.fastq.gz \
                out2=${name}_R2.flip.trim.fastq.gz \
                ref=${bbmap_adapters} \
                ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
                nullifybrokenquality=t \
                maq=10 minlen=25 \
                tpe tbo \
                literal=AAAAAAAAAAAAAAAAAAAAAAA \
                stats=${name}.trimstats.txt \
                refstats=${name}.refstats.txt
        """
    } else if (params.flip) {
        """
        echo ${name}


        reformat.sh -Xmx20g \
                t=16 \
                in=${name}.fastq.gz \
                out=${name}.flip.fastq.gz \
                rcomp=t

        
        bbduk.sh -Xmx20g \
                  t=16 \
                  in=${name}.flip.fastq.gz \
                  out=${name}.flip.trim.fastq.gz \
                  ref=${bbmap_adapters} \
                  ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
                  nullifybrokenquality=t \
                  maq=10 minlen=25 \
                  tpe tbo \
                  literal=AAAAAAAAAAAAAAAAAAAAAAA \
                  stats=${name}.trimstats.txt \
                  refstats=${name}.refstats.txt
        """
    }
        else if (!params.singleEnd) {
        """
        echo ${name}      

        bbduk.sh -Xmx20g \
                  t=16 \
                  in=${name}_R1.fastq.gz \
                  in2=${name}_R2.fastq.gz \
                  out=${name}_R1.trim.fastq.gz \
                  out2=${name}_R2.trim.fastq.gz \
                  ref=${bbmap_adapters} \
                  ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
                  nullifybrokenquality=t \
                  maq=10 minlen=25 \
                  tpe tbo \
                  literal=AAAAAAAAAAAAAAAAAAAAAAA \
                  stats=${name}.trimstats.txt \
                  refstats=${name}.refstats.txt
        """
    } else {
        """
        echo ${name}
        
        bbduk.sh -Xmx20g \
                  t=16 \
                  in=${name}.fastq.gz \
                  out=${name}.trim.fastq.gz \
                  ref=${bbmap_adapters} \
                  ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
                  nullifybrokenquality=t \
                  maq=10 minlen=25 \
                  tpe tbo \
                  literal=AAAAAAAAAAAAAAAAAAAAAAA \
                  stats=${name}.trimstats.txt \
                  refstats=${name}.refstats.txt
        """
    }
}


/*
 * STEP 2b - Trimmed FastQC
 */

process fastqc_trimmed {
    validExitStatus 0,1
    tag "$name"
    memory '4 GB'
    publishDir "${params.outdir}/qc/fastqc/", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(trimmed_reads) from trimmed_reads_fastqc

    output:
    file "*_fastqc.{zip,html,txt}" into trimmed_fastqc_results

    script:
    """
    echo ${name}

    fastqc $trimmed_reads
    """
}


/*
 * STEP 3 - Map reads to reference genome
 */

process hisat2 {
    tag "$name"
    validExitStatus 0
    cpus 32
    memory '100 GB'
    time '2h'
    publishDir "${params.outdir}/qc/hisat2_mapstats", mode: 'copy', pattern: "*.txt"

    input:
    //file(indices) from hisat2_indices
    val(indices_path) from hisat2_indices
    set val(name), file(trimmed_reads) from trimmed_reads_hisat2

    output:
    set val(name), file("*.sam") into hisat2_sam
    file("*.txt") into hisat2_mapstats    

  

    script:
    if (!params.singleEnd) {
        """
        echo ${name}
    
        hisat2  -p 32 \
                --very-sensitive \
                $rnastrandness \
                -x /scratch/Shares/dowell/genomes/hg38/HISAT2/genome \
                --pen-noncansplice 14 \
                --mp 1,0 \
                --sp 3,1 \
                -1 ${name}_R1.trim.fastq.gz \
                -2 ${name}_R2.trim.fastq.gz \
                --new-summary \
                > ${name}.sam \
                2> ${name}.hisat2_mapstats.txt                
        """
    } else {
        """
        echo ${name}
    
        hisat2  -p 32 \
                --very-sensitive \
                $rnastrandness \
                --pen-noncansplice 14 \
                --mp 1,0 \
                --sp 3,1 \
                -x ${indices_path}\
                -U ${trimmed_reads} \
                --new-summary \
                > ${name}.sam \
                2> ${name}.hisat2_mapstats.txt                
        """
    }
}


/*
 * STEP 4 - Convert to BAM format and sort
 */

process samtools {
    tag "$name"
    memory '100 GB'
    cpus 16
    publishDir "${params.outdir}" , mode: 'copy',
    saveAs: {filename ->
             if ((filename.indexOf(".bam") > 0) & !params.skipBAM)                                                                                                                             "mapped/bams/$filename"
        else if ((filename.indexOf(".bam.bai") > 0) & !params.skipBAM)                                                                                                                         "mapped/bams/$filename"
        else if (filename.indexOf("flagstat") > 0)                    "qc/mapstats/$filename"
        else if (filename.indexOf("millionsmapped") > 0)              "qc/mapstats/$filename"
    }

    input:
    set val(name), file(mapped_sam) from hisat2_sam

    output:
    set val(name), file("${name}.bam") into sorted_bam_ch
    set val(name), file("${name}.bam.bai") into sorted_bam_indices_ch
    set val(name), file("${name}.flagstat") into bam_flagstat
    set val(name), file("${name}.millionsmapped") into bam_milmapped_bedgraph

    script:
    if (!params.singleEnd) {
    """

    samtools view -@ 16 -bS -o ${name}.bamtmp ${mapped_sam}
    samtools sort -@ 16 ${name}.bamtmp > ${name}.bam
    samtools flagstat ${name}.bam > ${name}.flagstat
    samtools view -@ 16 -F 0x40 ${name}.bam | cut -f1 | sort | uniq | wc -l > ${name}.millionsmapped
    samtools index ${name}.bam ${name}.bam.bai
    #samtools view -@ 16 -C -T ${genome} -o ${name}.cram ${name}.bam
    #samtools sort -@ 16 -O cram ${name}.cramtmp > ${name}.cram
    #samtools index -c ${name}.cram ${name}.cram.crai
    """
    } else {
    """

    samtools view -@ 16 -bS -o ${name}.bam ${mapped_sam}
    samtools sort -@ 16 ${name}.bam > ${name}.sorted.bam
    samtools flagstat ${name}.sorted.bam > ${name}.flagstat
    samtools view -@ 16 -F 0x904 -c ${name}.sorted.bam > ${name}.millionsmapped
    samtools index ${name}.sorted.bam ${name}.sorted.bam.bai
    samtools view -@ 16 -C -T ${genome} -o ${name}.cram ${name}.sorted.bam
    samtools sort -@ 16 -O cram ${name}.cram > ${name}.sorted.cram
    samtools index -c ${name}.sorted.cram ${name}.sorted.cram.crai
    """
    }
}

sorted_bam_ch
   .into {sorted_bams_for_bedtools_bedgraph; sorted_bams_for_preseq; sorted_bams_for_rseqc; sorted_bams_for_dreg_prep; sorted_bams_for_pileup; sorted_bams_for_rseqc_count; sorted_bams_for_deeptools_bedgraph}

sorted_bam_indices_ch
    .into {sorted_bam_indices_for_bedtools_bedgraph; sorted_bam_indices_for_bedtools_normalized_bedgraph; sorted_bam_indicies_for_pileup; sorted_bam_indices_for_preseq; sorted_bam_indices_for_rseqc; sorted_bam_indices_for_rseqc_count; sorted_bams_indices_for_deeptools_bedgraph}

/*
 *STEP 5a - Plot the estimated complexity of a sample, and estimate future yields
 *         for complexity if the sample is sequenced at higher read depths.
 */

process preseq {
    tag "$name"
    memory '20 GB'
    time '8h'
    errorStrategy 'ignore'
    publishDir "${params.outdir}/qc/preseq/", mode: 'copy', pattern: "*.txt"

    input:
    set val(name), file(bam_file) from sorted_bams_for_preseq
    file(bam_indices) from sorted_bam_indices_for_preseq

    output:
    file("*.txt") into preseq_results

    script:
    """
    preseq c_curve -B -o ${name}.c_curve.txt \
           ${bam_file}

    preseq lc_extrap -B -o ${name}.lc_extrap.txt \
           ${bam_file}
    """
 }


/*
 *STEP 5b - Analyze read distributions using RSeQC
 */



/*
 *STEP 5c - Analyze coverage using pileup.sh
 */

process pileup {
    tag "$name"
    memory '50 GB'
    publishDir "${params.outdir}/qc/pileup", mode: 'copy', pattern: "*.txt"

    input:
    set val(name), file(bam_file) from sorted_bams_for_pileup
    file(bam_indices) from sorted_bam_indicies_for_pileup

    output:
    file("*.txt") into pileup_results

    script:
    """

    pileup.sh -Xmx20g \
              in=${bam_file} \
              out=${name}.coverage.stats.txt \
              hist=${name}.coverage.hist.txt
    """
 }


// process deeptools_bedgraph {
//     validExitStatus 0,143
//     echo true
//     tag "$name"
//     memory '80 GB'
//     time '4h'
//     cpus 32
//     publishDir "${params.outdir}/mapped/deep_tools/stranded_bedgraphs", mode: 'copy', pattern: "*{min,plu}.BedGraph"
//     publishDir "${params.outdir}/mapped/deep_tools/combined_bedgraphs", mode: 'copy', pattern: "*.BedGraph"
//     publishDir "${params.outdir}/mapped/deep_tools/rcc_bedgraphs", mode: 'copy', pattern: "${name}_rcc.BedGraph"
    
//     input:
//     set val(name), file(bam_file), file(millions_mapped) from bam_milmapped_deeptools_bedgraph
//     set val(name), file(bam_indices) from sorted_bams_indices_for_deeptools_bedgraph

//     output:
//     set val(name), file("*plu.BedGraph") into blank1
//     set val(name), file("*min.BedGraph") into blank2
//     set val(name), file("${name}.BedGraph") into blank3
//     set val(name), file("${name}_rcc.BedGraph") into tdf_deeptools
//     set val(name), file("${name}_plu_rcc.BedGraph") into blank5
//     set val(name), file("${name}_min_rcc.BedGraph") into blank6


//     script:

//     //////// From deep tools documentation
//     // Selects RNA-seq reads (single-end or paired-end) originating from genes on the given strand. 
//     // This option assumes a standard dUTP-based library preparation 
//     // (that is, –filterRNAstrand=forward keeps minus-strand reads, which originally came from genes
//     //  on the forward strand using a dUTP-based method). 
//     //  Consider using –samExcludeFlag instead for filtering by strand in other contexts.
//     // --extendReads not good for splicing, usually used on chipseq. but for dogs it might be good
//     """
//     bamCoverage -b ${bam_file} -o ${name}_plu.BedGraph -of "bedgraph" -p 32 --filterRNAstrand reverse --extendReads
//     bamCoverage -b ${bam_file} -o ${name}_tpm.min.BedGraph -of "bedgraph" -p 32 --filterRNAstrand forward --extendReads

//     awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' ${name}_tpm.min.BedGraph > ${name}_min.BedGraph
//     cat ${name}_plu.BedGraph ${name}_min.BedGraph > ${name}_unsorted.BedGraph

//     sortBed -i ${name}_unsorted.BedGraph > ${name}.BedGraph

//     python ${params.rcc} ${name}.BedGraph ${millions_mapped} ${name}_rcc.BedGraph 
//     python ${params.rcc} ${name}_plu.BedGraph ${millions_mapped} ${name}_plu_rcc.BedGraph 
//     python ${params.rcc} ${name}_min.BedGraph ${millions_mapped} ${name}_min_rcc.BedGraph 
//     """
// }

// process deep_tools_tdf {
//     tag "$name"
//     memory '200 GB'
//     time '1h'
//     // This often blows up due to a ceiling in memory usage, so we can ignore
//     // and re-run later as it's non-essential.
    
//     publishDir "${params.outdir}/mapped/deeptools/tdfs", mode: 'copy', pattern: "*.tdf"
//     input:
//     set val(name), file(normalized_bg) from tdf_deeptools
//     file(chrom_sizes) from chrom_sizes

//     output:
//     set val(name), file("*.tdf") into tiled_data_ch_deeptools

//     script:
//     """
//     igvtools toTDF ${normalized_bg} ${name}_rcc.tdf ${chrom_sizes}
//     """
//  }



process bedgraphs {
    validExitStatus 0,143
    tag "$name"
    memory '80 GB'
    time '8h'
    maxForks 20
    cpus 20    
    publishDir "${params.outdir}/mapped/bedgraphs_stranded", mode: 'copy', pattern: "*{min,plu}.BedGraph"
    publishDir "${params.outdir}/mapped/bam_stranded", mode: 'copy', pattern: "*{plu,min}.bam*"
    publishDir "${params.outdir}/mapped/wgs_stranded", mode: 'copy', pattern: "*.wgs.regions.bed.*"    
    //publishDir "${params.outdir}/mapped/pileup", mode: 'copy', pattern: "*pileup*"
    // publishDir "${params.outdir}/mapped/bedgraphs", mode: 'copy', pattern: "${name}.BedGraph"
    // publishDir "${params.outdir}/mapped/rcc_bedgraphs", mode: 'copy', pattern: "${name}_rcc.BedGraph"
    publishDir "${params.outdir}/mapped/rcc_bedgraphs_stranded", mode: 'copy', pattern: "*{min,plu}_rcc.BedGraph"

    input:
    set val(name), file(bam_file) from sorted_bams_for_bedtools_bedgraph
    set val(name), file(bam_indices) from sorted_bam_indices_for_bedtools_bedgraph
    set val(name), file(millions_mapped) from bam_milmapped_bedgraph

    output:
    //set val(name), file("*pileup.txt") into pileup_out
    set val(name), file("*plu.BedGraph") into plu_non_normalized_bedgraphs
    set val(name), file("*min.BedGraph") into min_non_normalized_bedgraphs
    set val(name), file("*.bam") into bamouts
    set val(name), file("*.bam.bai") into bamoutsbai
    // set val(name), file("${name}.BedGraph") into non_normalized_bedgraphs
    set val(name), file("${name}_rcc.BedGraph") into bedgraph_tdf
    set val(name), file("${name}_plu_rcc.BedGraph") into bedgraph_bigwig_plu
    set val(name), file("${name}_min_rcc.BedGraph") into bedgraph_bigwig_min


    //#http://davetang.org/wiki/tiki-index.php?page=SAMTools#Extracting_only_the_first_read_from_paired_end_BAM_files
    ///opt/samtools/0.1.19/samtools view -h -b -f 0x0040 ${od}${bamroot}.sorted.bam > ${bedgraphfortdfdir}${bamroot}.pairfirst.bam
    //#0x0040 is hexadecimal for 64 (i.e. 16 * 4), which is binary for 1000000, corresponding to the read in the first read pair.

    //#need to know the flag for the second strand
    //https://broadinstitute.github.io/picard/explain-flags.html
    //#128 means second in pair
    //#128 in hexadecimal is 0x0080
    //16 read minerse strand

    script:
    if (!params.singleEnd) {

        if (params.forward_stranded && !params.unStranded){

        """
            #https://www.biostars.org/p/92935/

        ########### FORWARD STRAND
        # 1. alignments of the second in pair if they map to the forward strand
        # 2. alignments of the first in pair if they map to the minerse  strand

        # Keep second in pair remove read minerse strand

        samtools view -b -f 128 -F 16 ${bam_file} > ${name}_plu1.bam
        samtools index ${name}_plu1.bam
        samtools view -b -f 80 ${bam_file} > ${name}_plu2.bam
        samtools index ${name}_plu2.bam

        samtools merge -f ${name}_plu.bam ${name}_plu1.bam ${name}_plu2.bam
        samtools index ${name}_plu.bam

        ########### REVERSE STRAND
        # 1. alignments of the second in pair if they map to the minerse strand
        # 2. alignments of the first in pair if they map to the forward strand

        # Flag 144 - keep
            #read minerse strand
            #second in pair

        samtools view -b -f 144 ${bam_file} > ${name}_min1.bam
        samtools index ${name}_min1.bam

        # Flag 64 first in pair       - keep
        # Flag 16 read minerse strand - remove
        samtools view -b -f 64 -F 16 ${bam_file} > ${name}_min2.bam
        samtools index ${name}_min2.bam

        samtools merge -f ${name}_min.bam ${name}_min1.bam ${name}_min2.bam
        samtools index ${name}_min.bam


        ########## SPLIT BAM #################
        genomeCoverageBed -bg -split -strand + -ibam ${name}_plu.bam > ${name}_pos.bed.tmp
        genomeCoverageBed -bg -split -strand - -ibam ${name}_min.bam > ${name}_neg.bed.tmp1

        awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' ${name}_neg.bed.tmp1 > ${name}_neg.bed.tmp2
        sortBed -i ${name}_pos.bed.tmp > ${name}_plu.BedGraph
        sortBed -i ${name}_neg.bed.tmp2 > ${name}_min.BedGraph



        python ${params.rcc} ${name}_plu.BedGraph ${millions_mapped} ${name}_plu_rcc.BedGraph
        python ${params.rcc} ${name}_min.BedGraph ${millions_mapped} ${name}_min_rcc.BedGraph


        cat ${name}_plu.BedGraph ${name}_min.BedGraph > ${name}.BedGraph.tmp
        sortBed -i ${name}.BedGraph.tmp > ${name}.BedGraph
        python ${params.rcc} ${name}.BedGraph ${millions_mapped} ${name}_rcc.BedGraph



        mosdepth -n -t 20 --fast-mode --by ${WINDOW} ${name}_plu.wgs ${name}_plu.bam
        mosdepth -n -t 20 --fast-mode --by ${WINDOW} ${name}_min.wgs ${name}_min.bam



        rm ${name}_min1.bam
        rm ${name}_min2.bam
        rm ${name}_plu2.bam
        rm ${name}_plu1.bam
        rm ${name}_pos.bed.tmp
        rm ${name}_neg.bed.tmp1
        rm ${name}_neg.bed.tmp2
        rm ${name}.BedGraph.tmp
        rm ${name}.BedGraph
        """

        } else if (params.reverse_stranded && !params.unStranded){

        """
            #https://www.biostars.org/p/92935/

        ########### FORWARD STRAND
        # 1. alignments of the second in pair if they map to the forward strand
        # 2. alignments of the first in pair if they map to the minerse  strand

        # Keep second in pair remove read minerse strand

        samtools view -b -f 64 -F 16 ${bam_file} > ${name}_plu1.bam
        samtools index ${name}_plu1.bam
        samtools view -b -f 144 ${bam_file} > ${name}_plu2.bam
        samtools index ${name}_plu2.bam

        samtools merge -f ${name}_plu.bam ${name}_plu1.bam ${name}_plu2.bam
        samtools index ${name}_plu.bam

        ########### REVERSE STRAND
        # 1. alignments of the second in pair if they map to the minerse strand
        # 2. alignments of the first in pair if they map to the forward strand

        # Flag 144 - keep
            #read minerse strand
            #second in pair

        samtools view -b -f 80 ${bam_file} > ${name}_min1.bam
        samtools index ${name}_min1.bam

        # Flag 64 first in pair       - keep
        # Flag 16 read minerse strand - remove
        samtools view -b -f 128 -F 16 ${bam_file} > ${name}_min2.bam
        samtools index ${name}_min2.bam

        samtools merge -f ${name}_min.bam ${name}_min1.bam ${name}_min2.bam
        samtools index ${name}_min.bam


        # Get pileup
        ###samtools mpileup -f ${genome} ${name}_plu.bam > ${name}_plu_pileup.txt
        ###samtools mpileup -f ${genome} ${name}_min.bam > ${name}_min_pileup.txt



        ########## SPLIT BAM #################
        genomeCoverageBed -bg -split -strand + -ibam ${name}_plu.bam > ${name}_pos.bed.tmp
        genomeCoverageBed -bg -split -strand - -ibam ${name}_min.bam > ${name}_neg.bed.tmp1

        awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' ${name}_neg.bed.tmp1 > ${name}_neg.bed.tmp2
        sortBed -i ${name}_pos.bed.tmp > ${name}_plu.BedGraph
        sortBed -i ${name}_neg.bed.tmp2 > ${name}_min.BedGraph

        python ${params.rcc} ${name}_plu.BedGraph ${millions_mapped} ${name}_plu_rcc.BedGraph
        python ${params.rcc} ${name}_min.BedGraph ${millions_mapped} ${name}_min_rcc.BedGraph


        cat ${name}_plu.BedGraph ${name}_min.BedGraph > ${name}.BedGraph.tmp
        sortBed -i ${name}.BedGraph.tmp > ${name}.BedGraph
        python ${params.rcc} ${name}.BedGraph ${millions_mapped} ${name}_rcc.BedGraph

        rm ${name}_min1.bam
        rm ${name}_min2.bam
        rm ${name}_plu2.bam
        rm ${name}_plu1.bam
        rm ${name}_pos.bed.tmp
        rm ${name}_neg.bed.tmp1
        rm ${name}_neg.bed.tmp2
        rm ${name}.BedGraph.tmp
        rm ${name}.BedGraph
        """
        }

    }   
 }






// process bedgraphs {
//     validExitStatus 0,143
//     tag "$name"
//     memory '80 GB'
//     time '4h'
//     publishDir "${params.outdir}/mapped/bedgraphs_stranded", mode: 'copy', pattern: "*{min,plu}.BedGraph"
//     publishDir "${params.outdir}/mapped/bedgraphs", mode: 'copy', pattern: "${name}.BedGraph"
//     publishDir "${params.outdir}/mapped/rcc_bedgraphs", mode: 'copy', pattern: "${name}_rcc.BedGraph"
//     publishDir "${params.outdir}/mapped/rcc_bedgraphs_stranded", mode: 'copy', pattern: "*{min,plu}.BedGraph"
//     publishDir "${params.outdir}/mapped/bedgraphs_bp", mode: 'copy', pattern: "*_bp_*"  

//     input:
//     set val(name), file(bam_file) from sorted_bams_for_bedtools_bedgraph
//     set val(name), file(bam_indices) from sorted_bam_indices_for_bedtools_bedgraph
//     set val(name), file(millions_mapped) from bam_milmapped_bedgraph

//     output:
//     set val(name), file("*plu.BedGraph") into plu_non_normalized_bedgraphs
//     set val(name), file("*min.BedGraph") into min_non_normalized_bedgraphs
//     set val(name), file("${name}.BedGraph") into non_normalized_bedgraphs
//     set val(name), file("${name}_rcc.BedGraph") into bedgraph_tdf
//     set val(name), file("${name}_plu_rcc.BedGraph") into bedgraph_bigwig_plu
//     set val(name), file("${name}_min_rcc.BedGraph") into bedgraph_bigwig_min


//     //#http://davetang.org/wiki/tiki-index.php?page=SAMTools#Extracting_only_the_first_read_from_paired_end_BAM_files
//     ///opt/samtools/0.1.19/samtools view -h -b -f 0x0040 ${od}${bamroot}.sorted.bam > ${bedgraphfortdfdir}${bamroot}.pairfirst.bam
//     //#0x0040 is hexadecimal for 64 (i.e. 16 * 4), which is binary for 1000000, corresponding to the read in the first read pair.

//     //#need to know the flag for the second strand
//     //https://broadinstitute.github.io/picard/explain-flags.html
//     //#128 means second in pair
//     //#128 in hexadecimal is 0x0080
//     //16 read reverse strand

//     script:
//     """
//         #https://www.biostars.org/p/92935/

//     ########### FORWARD STRAND
//     # 1. alignments of the second in pair if they map to the forward strand
//     # 2. alignments of the first in pair if they map to the reverse  strand

//     # Keep second in pair remove read reverse strand
//     samtools view -b -f 128 -F 16 ${bam_file} > ${name}_fwd1.bam
//     samtools index ${name}_fwd1.bam
//     samtools view -b -f 80 ${bam_file} > ${name}_fwd2.bam
//     samtools index ${name}_fwd2.bam

//     samtools merge -f ${name}_fwd.bam ${name}_fwd1.bam ${name}_fwd2.bam
//     samtools index ${name}_fwd.bam

//     ########### REVERSE STRAND
//     # 1. alignments of the second in pair if they map to the reverse strand
//     # 2. alignments of the first in pair if they map to the forward strand

//     # Flag 144 - keep
//         #read reverse strand
//         #second in pair

//     samtools view -b -f 144 ${bam_file} > ${name}_rev1.bam
//     samtools index ${name}_rev1.bam

//     # Flag 64 first in pair       - keep
//     # Flag 16 read reverse strand - remove
//     samtools view -b -f 64 -F 16 ${bam_file} > ${name}_rev2.bam
//     samtools index ${name}_rev2.bam

//     samtools merge -f ${name}_rev.bam ${name}_rev1.bam ${name}_rev2.bam
//     samtools index ${name}_rev.bam


//     ########## SPLIT BAM #################
//     genomeCoverageBed -d -strand + -ibam ${name}_fwd.bam > ${name}_bp_plu.BedGraph
//     genomeCoverageBed -d -strand - -ibam ${name}_rev.bam > ${name}_bp_min.BedGraph 


//     genomeCoverageBed -bg -split -strand - -ibam ${name}_fwd.bam > ${name}_pos.bed.tmp
//     genomeCoverageBed -bg -split -strand - -ibam ${name}_rev.bam > ${name}_neg.bed.tmp1

//     awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' ${name}_neg.bed.tmp1 > ${name}_neg.bed.tmp2
//     sortBed -i ${name}_pos.bed.tmp > ${name}_plu.BedGraph
//     sortBed -i ${name}_neg.bed.tmp2 > ${name}_min.BedGraph

//     cat ${name}_plu.BedGraph ${name}_min.BedGraph > ${name}.bedGraph.tmp
//     sortBed -i ${name}.bedGraph.tmp > ${name}.BedGraph

//     python ${params.rcc} ${name}_plu.BedGraph ${millions_mapped} ${name}_plu_rcc.BedGraph
//     python ${params.rcc} ${name}_min.BedGraph ${millions_mapped} ${name}_min_rcc.BedGraph
//     python ${params.rcc} ${name}.BedGraph ${millions_mapped} ${name}_rcc.BedGraph
//     """
//  }

process igvtools {
    tag "$name"
    memory '200 GB'
    time '1h'
    // This often blows up due to a ceiling in memory usage, so we can ignore
    // and re-run later as it's non-essential.
    errorStrategy 'ignore'
    publishDir "${params.outdir}/mapped/tdfs", mode: 'copy', pattern: "*.tdf"

    input:
    set val(name), file(normalized_bg) from bedgraph_tdf
    file(chrom_sizes) from chrom_sizes

    output:
    set val(name), file("*.tdf") into tiled_data_ch

    script:
    """
    igvtools toTDF ${normalized_bg} ${name}_rcc.tdf ${chrom_sizes}
    """
 }



/*
 * STEP 9 - MultiQC
 */
process multiQC {
    validExitStatus 0,1,143
    errorStrategy 'ignore'
    publishDir "${params.outdir}/multiqc/", mode: 'copy', pattern: "multiqc_report.html"
    publishDir "${params.outdir}/multiqc/", mode: 'copy', pattern: "*_data"

    when:
    !params.skipMultiQC

    input:
    file multiqc_config
    file (fastqc:'qc/fastqc/*') from fastqc_results.collect()
    file ('qc/fastqc/*') from trimmed_fastqc_results.collect()
    file ('qc/trimstats/*') from trim_stats.collect()
    file ('qc/mapstats/*') from bam_flagstat.collect()
    // file ('qc/rseqc/*') from rseqc_results.collect()
    file ('qc/preseq/*') from preseq_results.collect()
    file ('software_versions/*') from software_versions_yaml
    file ('qc/hisat2_mapstats/*') from hisat2_mapstats.collect()

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data" into multiqc_report_files

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''

    """
    multiqc . -f $rtitle $rfilename --config $multiqc_config
    """
}

/*
 * STEP 10 - Output Description HTML
 */
//
//process output_documentation {
//    tag "$prefix"
//    publishDir "${params.outdir}/pipeline_info/", mode: 'copy'
//
//    input:
//    file output_docs
//
//    output:
//    file "results_description.html"
//
//    script:
//    """
//    markdown_to_html.r $output_docs results_description.html
//    """
//}
//


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[SteadyFlow] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[SteadyFlow] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = params.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[SteadyFlow] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[SteadyFlow] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[SteadyFlow] Pipeline Complete"

}
