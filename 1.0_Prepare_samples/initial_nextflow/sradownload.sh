#!/bin/bash
#SBATCH --job-name=ANTH_DOWNLOAD# Job name
#SBATCH --mail-type=FAIL,END # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=markomelnicksupercomputer@gmail.com # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=8 # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=18:05:00 # Time limit hrs:min:sec
#SBATCH --partition short
#SBATCH --mem=40gb # Memory limit
#SBATCH --output=OUT.out
#SBATCH --error=ERR.out


module load sra/2.8.0

#Fastq
fastq-dump --split-3 SRR2038198 &
fastq-dump --split-3 SRR2038259 &
fastq-dump --split-3 SRR2038310 &
fastq-dump --split-3 SRR2038322 &
fastq-dump --split-3 SRR2038440 &
fastq-dump --split-3 SRR2038441 &
wait


mv SRR2038198_1.fastq SRR2038198_R1.fastq
mv SRR2038259_1.fastq SRR2038259_R1.fastq
mv SRR2038310_1.fastq SRR2038310_R1.fastq
mv SRR2038322_1.fastq SRR2038322_R1.fastq
mv SRR2038440_1.fastq SRR2038440_R1.fastq
mv SRR2038441_1.fastq SRR2038441_R1.fastq
mv SRR2038198_2.fastq SRR2038198_R2.fastq
mv SRR2038259_2.fastq SRR2038259_R2.fastq
mv SRR2038310_2.fastq SRR2038310_R2.fastq
mv SRR2038322_2.fastq SRR2038322_R2.fastq
mv SRR2038440_2.fastq SRR2038440_R2.fastq
mv SRR2038441_2.fastq SRR2038441_R2.fastq



gzip SRR2038198_R1.fastq &
gzip SRR2038259_R1.fastq &
gzip SRR2038310_R1.fastq &
gzip SRR2038322_R1.fastq &
gzip SRR2038440_R1.fastq &
gzip SRR2038441_R1.fastq &
gzip SRR2038198_R2.fastq &
gzip SRR2038259_R2.fastq &
gzip SRR2038310_R2.fastq &
gzip SRR2038322_R2.fastq &
gzip SRR2038440_R2.fastq &
gzip SRR2038441_R2.fastq &  
wait