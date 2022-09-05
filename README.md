# Nanopore
Useful scripts for Nanopore + makings of a pipeline

## WGS
Under nextflow_WGS you will find the Nextflow script for running the WGS pipeline on Nanopore data. You will need to get all the files from the 'nextflow_WGS' folder & copy them into your project directory (containing 'fastq_pass', 'fastq_fail' etc. folders).

- wgs.nf is the main Nextflow file. Here you will find all the processes. Feel free to fork and add more, if needed! 
- start_nextflow is a bash script to launch Nextflow. It first concatenates all the files into one fastq.gz file and then start the pipeline
- resume_nextflow can be used to simply re-launch Nextflow with all cached analyses (from where you left off)
- nextflow.config is a file used to input configurations into the pipeline. It's currently empty so that the pipeline can remain generic. 

### Note: the pipeline outputs files into alignment_output folder. The VCF (and/or other files) can then be used for downstream applications, as required. 

## RNAseq
Uunder nextflow_RNAseq you wil find similar scripts, but for RNAseq on Nanopore data. You will need to copy them into your project directory (containing 'fastq_pass', 'fastq_fail' etc. folders).

- wgs.nf is the main Nextflow file. Here you will find all the processes. Feel free to fork and add more, if needed! 
- barcodes.txt is an array with barcodes. It tells 'start_nextflow.sh' how to concatenate the barcoded read files. 
- start_nextflow is a bash script to launch Nextflow. It first concatenates all the files into one fastq.gz file and then start the pipeline
- resume_nextflow can be used to simply re-launch Nextflow with all cached analyses (from where you left off)
- nextflow.config is a file used to input configurations into the pipeline. It's currently empty so that the pipeline can remain generic.
- samples.tsv is a file needed for downstream anlayses (see below). It is used to differentiate between groups/conditions. 

### Note: RNAseq pipeline needs to be finished with Ballgown / DEseq2 / EDGER 
