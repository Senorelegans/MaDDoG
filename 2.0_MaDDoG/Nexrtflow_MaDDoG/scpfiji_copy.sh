#!/usr/bin/env bash

#source /Users/mame5141/.bashrc


# From fiji to summit



##### HOME
######/Users/mame5141/2021/Vilborg2015/Dogcatcher_nextflow


SCRATCH=mame5141@fiji.colorado.edu:/scratch/Users/mame5141/2021/Vilborg2015/Dogcatcher2/NEXTFLOW_MaDDoG

# scp -r bin/*.py ${SCRATCH}/bin
# scp -r bin/*.R ${SCRATCH}/bin
# # scp ../*.py ${SCRATCH}/bin 
# scp -r main.nf ${SCRATCH}
# scp -r dogcatcher_nextflow_RNAseq.sh ${SCRATCH}

# scp -r ${SCRATCH}/VILBORG_cpu_Region_WINDOW_SlidingWindow_5000_DogGenePercCov_1/MaDDoG/gtf bin
scp -r ${SCRATCH}/VILBORG_cpu_Region_WINDOW_SlidingWindow_5000_DogGenePercCov_1/MaDDoG .