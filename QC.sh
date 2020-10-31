#! /bin/bash

#PBS -l nodes=1:ppn=8
#PBS -N Asm
#PBS -o log.out
#PBS -e err.out
#PBS -q four
cd $PBS_O_WORKDIR

/gpfshome/home/Changq/miniconda3/envs/rna/bin/fastqc -o ../QC_res/ *.gz
