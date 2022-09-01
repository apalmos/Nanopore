split_fastq = Channel.fromPath("$PWD/fastq_split/barcode*")
reference = file("/scratch/prj/sgdp_nanopore/Resources/hg38.fa")

process fastqc {

  publishDir 'alignment_output/fastqc/', mode: 'copy'

  input:
  path query_file from split_fastq

  output:
  path "*"

  script:
  """
  module load fastqc
  fastqc ${query_file}
  """
}

split_fastq = Channel.fromPath("$PWD/fastq_split/barcode*")

process nanoplot {

  publishDir 'alignment_output/nanoplot/', mode: 'copy'

  input:
  path query_file from split_fastq

  output:
  path "${query_file}_*"

  script:
  """
  NanoPlot --fastq ${query_file} -o ${query_file}_*
  """
}

split_fastq = Channel.fromPath("$PWD/fastq_split/barcode*")

process nanofilt {

  scratch "$PWD/work"

  input:
  path query_file from split_fastq

  output:
  file "${query_file}.filtered"  into filter

  script:
  """
  gunzip -c ${query_file} | NanoFilt -l 500 --headcrop 10 | gzip - > ${query_file}.filtered
  """
}

process aligAndconvert {

  publishDir 'alignment_output/'

  input:
  path query_file from filter

  output:
  file "${query_file}.bam" into genomes

  script:
  """
  minimap2 -ax splice ${reference} ${query_file} | samtools sort > ${query_file}.bam
  """
}

genomes.into {
  genomes1
  genomes2
  genomes3
  genomes4
  genomes5
  genomes6
  genomes7
}

process bedtools {

  publishDir 'alignment_output/bigBedWig'

  input:
  path query_file from genomes1

  output:
  file "${query_file}.bed" into bed
  file "${query_file}.big_bed" into bigbed
  file "${query_file}.bedGraph" into bedgraph
  file "${query_file}.bigWig" into bigwig

  script:
  """
  bedtools bamtobed -i ${query_file} > ${query_file}.bed
  /scratch/prj/sgdp_nanopore/software/bedSort ${query_file}.bed ${query_file}.bed
  /scratch/prj/sgdp_nanopore/software/bedToBigBed ${query_file}.bed /scratch/prj/sgdp_nanopore/Resources/hg38.chrom.sizes ${query_file}.big_bed
  bedtools genomecov -ibam ${query_file} -bg | bedtools sort > ${query_file}.bedGraph
  /scratch/prj/sgdp_nanopore/software/bedGraphToBigWig ${query_file}.bedGraph /scratch/prj/sgdp_nanopore/Resources/hg38.chrom.sizes ${query_file}.bigWig
  """
}

process stats {

  publishDir 'alignment_output/stats'

  input:
  file query_file from genomes2

  output:
  file "${query_file}.stats" into stats
  file "${query_file}.flagstat" into flagstat

  script:
  """
  samtools stats ${query_file} |grep ^SN | cut -f 2- > ${query_file}.stats
  samtools flagstat ${query_file} > ${query_file}.flagstat
  """
}

process assemble_transcriptome {

  publishDir 'alignment_output/assembly'

  input:
  file query_file from genomes3

  output:
  file "${query_file}.gtf" into assembly

  script:
  """
  /scratch/prj/sgdp_nanopore/software/stringtie/stringtie -L -G /scratch/prj/sgdp_nanopore/Resources/hg38.knownGene.gtf -o ${query_file}.gtf ${query_file}
  """
}

process merge_transcriptome {

  publishDir 'alignment_output/merged_transcriptome'

  input:
  file gtfs from assembly.collect()

  output:
  file "merged.gtf" into merged

  script:
  """
  /scratch/prj/sgdp_nanopore/software/stringtie/stringtie --merge ${gtfs} -G /scratch/prj/sgdp_nanopore/Resources/hg38.knownGene.gtf -o merged.gtf
  """
}

process assemble_final {

  publishDir 'alignment_output/final_assembly'

  input:
  file merged from merged
  file query_file from genomes4

  output:
  file "${query_file}.gtf" into final_assembly

  script:
  """
  /scratch/prj/sgdp_nanopore/software/stringtie/stringtie -L -G ${merged} -o ${query_file}.gtf ${query_file}
  """
}

process gene_abundance {

  publishDir 'alignment_output/gene_abundance'

  input:
  file merged from merged
  file query_file from genomes5

  output:
  file "${query_file}_gene_abundance.gtf" into gene_abundance

  script:
  """
  /scratch/prj/sgdp_nanopore/software/stringtie/stringtie -L -A -G ${merged} -o ${query_file}_gene_abundance.gtf ${query_file}
  """
}

process count_genes {

  publishDir 'alignment_output/count_genes'

  input:
  file merged from merged
  file query_file from genomes6.collect()

  output:
  file "barcode_counts_gene.txt" into count_genes

  script:
  """
  /scratch/prj/sgdp_nanopore/software/subread-2.0.3-source/bin/featureCounts -L -O -g gene_id -t exon -a ${merged} -o barcode_counts_gene.txt ${query_file}
  """
}

process count_transcripts {

  publishDir 'alignment_output/count_transcripts'

  input:
  file merged from merged
  file query_file from genomes7.collect()

  output:
  file "barcode_counts_transcript.txt" into count_transcripts

  script:
  """
  /scratch/prj/sgdp_nanopore/software/subread-2.0.3-source/bin/featureCounts -L -O -f --primary --fraction -F GTF -g transcript_id -t transcript --extraAttributes gene_id -a ${merged} -o barcode_counts_transcript.txt ${query_file}
  """
}
