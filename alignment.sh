#! /bin/bash

#PBS -l nodes=1:ppn=8
#PBS -N Asm
#PBS -o log.out
#PBS -e err.out
#PBS -q four
cd $PBS_O_WORKDIR

# alignment reference genome


bowtie -v 0 -p 8 -x reference_index/ref -a -m 50 -q data.fastq -S data.sam 2>align_rate.txt
