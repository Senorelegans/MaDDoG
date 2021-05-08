#!/bin/bash
#SBATCH --job-name=dogcatcher_nextflow# Job name
#SBATCH --mail-type=FAIL,END # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=markomelnicksupercomputer@gmail.com # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=12 # Number of CPU (processer cores i.e. tasks)
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


ml R/3.5.1


eval "$(conda shell.bash hook)"
conda activate /Users/mame5141/miniconda2/envs/tfp
########################################################################
################## PRINT JOB INFO ######################################


PROJECT=/scratch/Users/mame5141/2021/Vilborg2015/Dogcatcher2

MAIN=/scratch/Users/mame5141/2021/Vilborg2015/Dogcatcher2/NEXTFLOW_MaDDoG

GTF="/scratch/Users/mame5141/2021/Vilborg2015/Dogcatcher2/gtf/GCF_000001405.39_GRCh38.p13_genomic.sorted.gff"

annotation_type="REFSEQGFF"
CPUS=12
COMPARISON="VILBORG_cpu"
GENEREGION=WINDOW
SLIDINGWINDOW=5000
DOGGENEPERCCOV=1
TITLE=${COMPARISON}_Region_${GENEREGION}_SlidingWindow_${SLIDINGWINDOW}_DogGenePercCov_${DOGGENEPERCCOV}
mkdir -p $PROJECT
BIN=${MAIN}/bin


python ${BIN}/nb_10_DogcatcherFlatten.py --annotation_in ${GTF} --file_type ${annotation_type}

python ${BIN}/nb_MaDDoG.py \
--Comparison ${COMPARISON} \
--OutFolder OUT \
--annotation_in ${GTF} \
--wgs_path ${PROJECT}"/data/WGS/Vilborg2015/wgs_stranded/normalized" \
--sample_table ${MAIN}"/Vilborg_Sample_Table.tsv" \
--cpus ${CPUS} \
--DOGCATCHER_ONLY "True" \
--GeneRegion ${GENEREGION} \
--RunDOG True \
--DogGenePercCov 1 \
--GeneRegionCoverageMin 5 \
--InsideGeneSize 1000 \
--SlidingWindow 5000 \
--reads_window 50


## RUN MADDOG in parallel
#######################################################################

nextflow run ${MAIN}/main.nf -resume -profile slurm \
        --bin = ${BIN} \
        --workdir ${PROJECT}/temp \
        --outdir ${MAIN}/${TITLE} \
        --email "markomelnicksupercomputer@gmail.com" \
        --CHR_NAMES "/scratch/Users/mame5141/2021/Vilborg2015/Dogcatcher2/gtf/GCF_000001405.39_GRCh38.p13_genomic.sorted_flat_CHROMNAMES.txt" \
        --BAMPATH "/scratch/Users/mame5141/2021/Vilborg2015/Dogcatcher_nextflow/mapped/bams" \
        --max_num_states 6 \
        --bayes_factor_thresh 1.3 \
        --iterations 40 \
        --convolution_window 10