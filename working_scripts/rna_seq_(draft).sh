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






















module load fastqc

#fastqc
mkdir fastqc
fastqc -o fastqc barcode1.fastq.gz

#nanoplot
NanoPlot --fastq barcode1.fastq.gz -o nanoplot

#nanofilt
gunzip -c barcode1.fastq.gz | NanoFilt -l 500 --headcrop 10 | gzip - > barcode1_filtered.fastq.gz

#fastqc again to compare
mkdir fastqc_filtered
fastqc -o fastqc_filtered barcode1_filtered.fastq.gz

#align
minimap2 \
-ax splice \
-R '@RG\tID:barcode1\tSM:RNAseq_test' \
/scratch/prj/sgdp_nanopore/Resources/hg38.fa \
barcode1_filtered.fastq.gz | samtools sort > sorted.barcode1.bam

####get statistics
samtools stats sorted.barcode1.bam |grep ^SN | cut -f 2-
samtools flagstat sorted.barcode1.bam

#get chrom sizes
ln -s /scratch/prj/sgdp_nanopore/Resources/hg38.fa .
bn=$(basename hg38.fa)
samtools faidx $bn
cut -f1,2 ${bn}.fai > chrom.sizes

####bedtools

#stringtie to assmelbe transcriptome
 ../../../software/stringtie/stringtie -L -G hg38.knownGene.gtf -o barcode1.gtf sorted.barcode1.bam

#stringtie merge transcriptome
 ../../../software/stringtie/stringtie --merge barcode1.gtf barcode2.gtf -G hg38.knownGene.gtf -o stringtie.merged.gtf

#stringtie assembly final
 ../../../software/stringtie/stringtie -L -G stringtie.merged.gtf -o barcode1.stringtie_final.gtf sorted.barcode1.bam

  ../../../software/stringtie/stringtie -L -G stringtie.merged.gtf -o barcode2.stringtie_final.gtf sorted.barcode2.bam

#Assemble transcriptome and compute RNA-seq expression
 ../../../software/stringtie/stringtie -L -A -G stringtie.merged.gtf -o barcode1.stringtie_transcript_final.gtf sorted.barcode1.bam

  ../../../software/stringtie/stringtie -L -A -G stringtie.merged.gtf -o barcode2.stringtie_transcript_final.gtf sorted.barcode2.bam

#feature count of genes
 ../../../software/subread-2.0.3-source/bin/featureCounts -L -O -g gene_id -t exon -a stringtie.merged.gtf -o barcode_counts_gene.txt sorted.barcode1.bam sorted.barcode2.bam

#feature count of transcripts
../../../software/subread-2.0.3-source/bin/featureCounts -L -O -f --primary --fraction -F GTF -g transcript_id -t transcript --extraAttributes gene_id -a stringtie.merged.gtf -o barcode_counts_transcript.txt sorted.barcode1.bam sorted.barcode2.bam


#Ballgown is used to calculate differential transcript and gene expression levels and test them for significant differences.

##first make a table that can be read into R
../../../software/tablemaker-2.1.1.Linux_x86_64/tablemaker -p 4 -q -W --library-type fr-firststrand -G stringtie.merged.gtf -o barcode1_output sorted.barcode1.bam

../../../software/tablemaker-2.1.1.Linux_x86_64/tablemaker -p 4 -q -W --library-type fr-firststrand -G stringtie.merged.gtf -o barcode2_output sorted.barcode2.bam






tail -n +2 barcode_counts_transcript.txt > gene1
awk '{print $1, $7}' gene1 > barcode1_gene
awk '{SUM+=$2}END{print SUM}' barcode1_gene

tail -n +2 barcode2_counts_gene.txt > gene1
awk '{print $1, $7}' gene1 > barcode2_gene
awk '{SUM+=$2}END{print SUM}' barcode2_gene


tail -n +2 barcode1_counts_transcript.txt > temp
awk '{print $1, $8}' temp > temp2

awk '{SUM+=$2}END{print SUM}' test2
awk '{SUM+=$2}END{print SUM}' temp2



#Differential expression analysis using DESEQ and EDGER. Merge the results of the analysis in a single csv file.

to do









#params.sampleName="barcode01"
#params.fastq_path = "$PWD/fastq_pass/${params.sampleName}/*.gz"
#split_fastq = Channel.fromPath(params.fastq_path)










#!/bin/bash
#SBATCH -n 10
#SBATCH --mem-per-cpu=100G
#SBATCH -t 48:00:00
#SBATCH -o alignment_1KG

files=(/scratch/prj/sgdp_nanopore/Projects/04_Single_cell_RNAseq_AML/sample2/20220427_1056_2E_PAK00128_d5620da9/fastq_pass/*)
idx=0
len=400
pcount=0

while [ $idx -le ${#files[@]} ] ; do
  cat "${files[@]:idx:len}" > ./QC/fastq${pcount}.gz
  ((idx+=len))
  ((pcount++))
done

#!/bin/bash
#SBATCH -n 10
#SBATCH --mem-per-cpu=100G
#SBATCH -t 48:00:00
#SBATCH -o alignment_1KG

directory=$(pwd)
file=$1
reference=(/scratch/prj/sgdp_nanopore/Resources/human_g1k_v37.fasta)


conda create --name rnaseq_nanopore

minimap_options:
    directRNA: '-ax splice -uf -k14'
    cDNA: '-ax splice'
    cDNAStranded: '-ax splice -uf'
    directcDNA: '-ax splice -uf -k14'

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz

#Nanoplot
NanoPlot -o {params.OUT_DIR} --fastq {input}

fastqc:
fastqc -o {params.OUT_DIR} {input}


####align for each barcode
minimap2 \
-ax splice -uf -k14 \
-R '@RG\tID:barcode1\tSM:RNAseq_test' \
/scratch/prj/sgdp_nanopore/Resources/hg38.fa \
barcode1.fastq.gz | samtools sort > sorted.barcode1.bam

minimap2 \
-ax splice -uf -k14 \
-R '@RG\tID:barcode2\tSM:RNAseq_test' \
/scratch/prj/sgdp_nanopore/Resources/hg38.fa \
barcode2.fastq.gz | samtools sort > sorted.barcode2.bam



#check RG tag
samtools view -H barcode1.bam | grep '^@RG'
samtools view -H barcode2.bam | grep '^@RG'


####get statistics
samtools stats sorted.barcode1.bam |grep ^SN | cut -f 2-
samtools stats sorted.barcode2.bam |grep ^SN | cut -f 2-

samtools flagstat sorted.barcode1.bam
samtools flagstat sorted.barcode2.bam

#get chrom sizes
ln -s /scratch/prj/sgdp_nanopore/Resources/hg38.fa .
bn=$(basename hg38.fa)
samtools faidx $bn
cut -f1,2 ${bn}.fai > chrom.sizes

####bedtools
bedtools bamtobed -i sorted.barcode1.bam > sorted.barcode1.bed
./bedSort sorted.barcode1.bed sorted.barcode1.bed
./bedToBigBed sorted.barcode1.bed chrom.sizes barcode1.big_bed
bedtools genomecov -ibam sorted.barcode1.bam -bg | bedtools sort > barcode1.bedGraph
./bedGraphToBigWig barcode1.bedGraph chrom.sizes barcode1.bigWig

bedtools bamtobed -i sorted.barcode2.bam > sorted.barcode2.bed
./bedSort sorted.barcode2.bed sorted.barcode2.bed
./bedToBigBed sorted.barcode2.bed chrom.sizes barcode2.big_bed
bedtools genomecov -ibam sorted.barcode2.bam -bg | bedtools sort > barcode2.bedGraph
./bedGraphToBigWig barcode2.bedGraph chrom.sizes barcode2.bigWig

#stringtie to assmelbe transcriptome
./stringtie/stringtie -L -G hg38.knownGene.gtf -o barcode1.gtf sorted.barcode1.bam
./stringtie/stringtie -L -G hg38.knownGene.gtf -o barcode2.gtf sorted.barcode2.bam


./stringtie/stringtie --merge barcode1.gtf barcode2.gtf -G hg38.knownGene.gtf -o stringtie.merged.gtf


./stringtie/stringtie -L -G stringtie.merged.gtf -o barcode1.stringtie_final.gtf sorted.barcode1.bam
./stringtie/stringtie -L -G stringtie.merged.gtf -o barcode2.stringtie_final.gtf sorted.barcode2.bam


./subread-2.0.3-source/bin/featureCounts -L -O -g gene_id -t exon -a stringtie.merged.gtf -o barcode1_counts_gene.txt sorted.barcode1.bam

./subread-2.0.3-source/bin/featureCounts -L -O -g gene_id -t exon -a stringtie.merged.gtf -o barcode2_counts_gene.txt sorted.barcode2.bam

./subread-2.0.3-source/bin/featureCounts -L -O -f --primary --fraction -F GTF -g transcript_id -t transcript --extraAttributes gene_id -a stringtie.merged.gtf -o barcode1_counts_transcript.txt sorted.barcode1.bam

./subread-2.0.3-source/bin/featureCounts -L -O -f --primary --fraction -F GTF -g transcript_id -t transcript --extraAttributes gene_id -a stringtie.merged.gtf -o barcode2_counts_transcript.txt sorted.barcode2.bam






#merge files
picard MergeSamFiles \
-I barcode1.bam \
-I barcode2.bam \
-O merged_files.bam

#sort file
picard SortSam \
-I merged_files.bam \
-O sorted.bam \
-SORT_ORDER coordinate

NanoStat --bam sorted.bam > readreport.txt

#quality control metrics using RNA-SeQC
picard CollectRnaSeqMetrics \
-I sorted.bam \
-O output.RNA_Metrics \
-REF_FLAT refFlat.txt \
-STRAND SECOND_READ_TRANSCRIPTION_STRAND

#Create wiggle tracks
bedtools genomecov -bg -ibam sorted.bam > sorted.bedGraph
./bedSort sorted.bedGraph resorted.bedGraph

./bedGraphToBigWig resorted.bedGraph hg38.chrom.sizes out.bw


bedtools bamtobed -i sorted.bam > sorted.bed
./bedSort sorted.bed resorted.bed


./bedToBigBed resorted.bed hg38.chrom.sizes out.bb

#cut -f1-3,5 resorted.bed > resorted.bedgraph
#bedtools merge -d 1000 -c 4 -o sum -i resorted.bedgraph > resorted.merged.bedgraph



###Count reads using HTseq
htseq-count sorted.bam gencode.v41.annotation.gtf

python -m HTSeq.scripts.count [options] <alignment_files> <gtf_file>

pip install numpy==1.20.0



./stringtie/stringtie -L -G gencode.v41.annotation.gtf -o test.stringtie.gtf sorted.bam




./subread-2.0.3-source/bin/featureCounts -O -L \
-a test.stringtie.gtf \
-o example_featureCounts_output.txt \
sorted.bam


./subread-2.0.3-source/bin/featureCounts -L -O -g gene_id -t exon -a test.stringtie.gtf -o example_featureCounts_output.txt sorted.bam

#raw count metrics



#stringtie to assmelbe transcriptome



#strintie to merge assemblies to a master transcriptome reference



#stringtie to assemble transcriptome and compute rnaseq expression



#ballgown to calculate differential transcript and gene expresion levels and test for significant differences



#DESEQ2 and EDGER for gene expression analysis (differential gene expression)

















pyBigWig

java -jar RNA-SeQC_v1.1.8.jar \
-s "Sample_ID|sorted.bam|Notes" \
-t test.gtf \
-r /scratch/prj/sgdp_nanopore/Resources/hg38.fa \
-strat gc -gc gencode.v7.gc.txt \
-BWArRNA human_all_rRNA.fasta \
-o ./rna_qc_out



java -jar RNASeQC.jar -n 1000 -s "TestId|ThousandReads.bam|TestDesc" -t gencode.v7.annotation_goodContig.gtf -r Homo_sapiens_assembly19.fasta -o ./testReport/ -strat gc -gc gencode.v7.gc.txt



#create bigWig
#bamCoverage -b all_reads_filtered.bam -o coverage.bw
#create bigBED
conda install -c aleg -c anaconda -c bioconda -c conda-forge pycoqc=2.5.2



###sort by position

#Stringtie
