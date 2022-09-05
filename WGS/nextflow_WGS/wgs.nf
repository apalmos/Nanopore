#! /usr/bin/env nextflow

split_fastq = Channel.fromPath("$PWD/fastq_split/fast*")
reference = file("/scratch/prj/sgdp_nanopore/Resources/hg38.fa")

process fastqc {

  publishDir 'alignment_output/fastqc/', mode: 'copy'

  input:
  path query_file from split_fastq.collect()

  output:
  path "*"

  script:
  """
  cat ${query_file} > whole_run.fastq.gz
  module load fastqc
  fastqc whole_run.fastq.gz
  """
}

split_fastq = Channel.fromPath("$PWD/fastq_split/fast*")

process aligAndconvert {

  scratch "$PWD/work"

  input:
  path query_file from split_fastq

  output:
  path "${query_file}.bam" into genomes

  script:
  """
  minimap2 -a --MD -x map-ont -u b -t 10 ${reference} ${query_file} | samtools view -b -o ${query_file}.bam
  """
}

process mergeBam {

  input:
  file sorted from genomes.collect()

  output:
  file "mergedbam.bam" into merged

  script:
  """
  samtools merge mergedbam.bam ${sorted}
  """
}

process sortBam {
  publishDir 'alignment_output/'

  input:
  file merged_file from merged

  output:
  path "sortedbam.bam" into sorted

  script:
  """
  samtools sort ${merged_file} -o sortedbam.bam
  samtools index sortedbam.bam
  """
}

process nanoStat {
  publishDir 'alignment_output/'

  input:
  file sortedbam from sorted

  when:
  sorted.name =~ /sortedbam.*/

  output:
  file "readreport.txt" into stats

  script:
  """
  NanoStat --bam ${sortedbam} > readreport.txt
  """
}

process depthSam {
  publishDir 'alignment_output/'

  input:
  file sortedbam from sorted

  when:
  sorted.name =~ /sortedbam.*/

  output:
  file "depth_all.txt" into depth_all
  file "depth_reads.txt" into depth_reads

  script:
  """
  samtools depth -a ${sortedbam} | awk '{sum+=\$3} END { print "Average = ",sum/NR}' > depth_all.txt
  samtools depth ${sortedbam} | awk '{sum+=\$3} END { print "Average = ",sum/NR}' > depth_reads.txt
  """
}

process sniffles {
  publishDir 'alignment_output/'

  input:
  tuple val(all_files) from sorted

  output:
  file "variants.vcf" into variants

  script:
  """
  sniffles --input ${all_files} --vcf variants.vcf --reference ${reference}
  """
}

process nanPlot {
  publishDir 'alignment_output/'

  input:
  file sortedbam from sorted

  when:
  sorted.name =~ /sortedbam.*/

  output:
  path 'bamplots_*' into nanoplot

  script:
  """
  NanoPlot -t 12 --color yellow --bam ${sortedbam} --downsample 10000 -o bamplots_downsampled
  """
}
