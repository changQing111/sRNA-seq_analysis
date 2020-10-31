#! /bin/bash

#PBS -l nodes=1:ppn=8
#PBS -N Asm
#PBS -o log.out
#PBS -e err.out
#PBS -q four
cd $PBS_O_WORKDIR

#| /gpfshome/home/Changq/miniconda3/envs/rna/bin/samtools view -bS > ${id}.bam

#~/rice_virus/RRSV_index/RRSV


cat file.txt|while read id;do /gpfshome/home/Changq/miniconda3/envs/rna/bin/bowtie2 -x /gpfshome/home/Changq/MSU7/bt2_index/Os -U ../Cleandata/trim_${id}.fq.gz -S ${id}.sam;done
