#!/bin/bash
#SBATCH -n 10
#SBATCH --mem-per-cpu=100G
#SBATCH -t 48:00:00
#SBATCH -o 5_scaffold

directory=`cat directory.txt`
file=$1

longstitch run draft=${directory}/assembly/dbg.cns \
reads=${directory}/QC/${file} \
G=3e9 \
> ${directory}/stitch/longstitch
