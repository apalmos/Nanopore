#!/bin/bash
#SBATCH -n 10
#SBATCH --mem-per-cpu=50G
#SBATCH -t 48:00:00

nextflow run wgs.nf -resume
