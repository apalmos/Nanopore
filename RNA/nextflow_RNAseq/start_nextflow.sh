#!/bin/bash
#SBATCH -n 10
#SBATCH --mem-per-cpu=50G
#SBATCH -t 48:00:00

module load nextflow
mkdir fastq_split

while read p; do

  echo "$p"
  cat $PWD/fastq_pass/$p/* > $PWD/fastq_split/$p.gz

done < barcodes.txt

nextflow run rnaseq.nf
