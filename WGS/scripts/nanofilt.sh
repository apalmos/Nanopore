#!/bin/bash
#SBATCH -n 20
#SBATCH --mem-per-cpu=19G
#SBATCH -t 48:00:00
#SBATCH -o 1_QC

file=$1
directory=`cat directory.txt`

gunzip -c ${directory}/fastq_pass/${file} | NanoFilt \
-l 500 \
--headcrop 10 \
> ${directory}/QC/output/${file}_nanofilt_trimmed.fastq

cat ${directory}/QC/output/${file}_nanofilt_trimmed.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > ${directory}/QC/output/${file}_nanofilt_trimmed.fa
