#! /bin/bash

#PBS -l nodes=1:ppn=8
#PBS -N Asm
#PBS -o log.out
#PBS -e err.out
#PBS -q four
cd $PBS_O_WORKDIR

export MYPATH="/gpfshome/home/Changq/miniconda3/envs/rna/bin/"
cat file.txt|while read id;do $MYPATH/samtools view -bS ${id}.sam|$MYPATH/samtools sort -o ${id}_sort.bam;done
