#!/bin/bash
#SBATCH -n 10
#SBATCH --mem-per-cpu=50G
#SBATCH -t 48:00:00

module load nextflow
mkdir fastq_split

files=($PWD/fastq_pass/*)                                 # list of files
idx=0                                       # start index of actual package
len=400                                     # files per package
pcount=0                                    # package counter

while [ $idx -le ${#files[@]} ] ; do
  cat "${files[@]:idx:len}" > ./fastq_split/fastq${pcount}.gz  # process subarray
  ((idx+=len))                              # start of next package
  ((pcount++))                              # number of next package
done

nextflow run wgs.nf
