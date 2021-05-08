#!/bin/bash
#SBATCH --job-name=covid_nextflow# Job name
#SBATCH --mail-type=FAIL,END # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=markomelnicksupercomputer@gmail.com # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=1 # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=12:05:00 # Time limit hrs:min:sec
#SBATCH --partition short
#SBATCH --mem=20gb # Memory limit
#SBATCH --output=nextflow_OUT.out
#SBATCH --error=nextflow_ERR.err



################### SET VARIABLES ######################################

export PATH=~:$PATH
export PATH=~/.local/bin:$PATH
export PATH=~/.local:$PATH

########################################################################
################### LOAD NECESSARY MODULES #############################

module load sra/2.8.0
module load samtools/1.8
module load hisat2/2.1.0
module load bedtools/2.25.0
module load gcc/7.1.0
module load seqkit/0.9.0
module load fastqc/0.11.8
module load bbmap/38.05
module load igvtools/2.3.75
module load preseq/2.0.3
module load mpich/3.2.1
ml gcc/7.1.0

source /Users/mame5141/.bashrc2

########################################################################
################## PRINT JOB INFO ######################################

#PROJECT='/scratch/Users/mame5141/2019/RepEditing/nextflow_out/'
PROJECT="/scratch/Users/mame5141/2021/Vilborg2015/Dogcatcher_nextflow"

mkdir -p $PROJECT

printf "Sample ID: $ROOTNAME"
printf "\nDirectory: $PROJECT"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date
printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

#######################################################################
MAIN=/scratch/Users/mame5141/2021/Vilborg2015/Dogcatcher_nextflow
BIN=${MAIN}/bin
nextflow run ${MAIN}/main.nf -resume -profile slurm \
        --bin = ${BIN} \
        --workdir ${PROJECT}/temp \
        --outdir ${PROJECT} \
        --fastqs "/Users/mame5141/2021/Vilborg2015/Dogcatcher_nextflow/fastqs/*_{R1,R2}.fastq.gz" \
        --reads "/Users/mame5141/2021/Vilborg2015/Dogcatcher_nextflow/fastqs/*_{R1,R2}.fastq.gz" \
        --genome "/scratch/Shares/dowell/genomes/hg38/hg38.fa" \
        --chrom_sizes "/scratch/Shares/dowell/genomes/hg38/hg38.chrom.sizes" \
        --hisat2_indices = "/scratch/Shares/dowell/genomes/hg38/HISAT2/genome" \
        --genome_refseq = "/scratch/Shares/dowell/genomes/hg38/hg38_refseq.bed" \
        --email "markomelnicksupercomputer@gmail.com" \
        --WINDOW 50 \
        --forward_stranded true \
        --reverse_stranded false \
        --unStranded false

#        --genome_refseq ${MAIN}"/gtf/hg38_refseq_EnsembleLIKE.gtf" \