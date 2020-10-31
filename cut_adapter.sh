#! /bin/bash

#PBS -l nodes=1:ppn=8
#PBS -N Asm
#PBS -o log.out
#PBS -e err.out
#PBS -q four
cd $PBS_O_WORKDIR


cat file_list.txt|while read id;do /gpfshome/home/Changq/miniconda3/bin/cutadapt -a AGATCGAAGAG \
--quality-base 33 -m 10 -q 20, 20 -o trim_${id} /gpfshome/home/Changq/Project/guoyi/raw_data/$id;done
