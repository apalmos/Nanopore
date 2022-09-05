srun -p cpu --pty /bin/bash

pkg <- "pinfsc50"
vcf_file <- system.file("extdata", "pinf_sc50.vcf.gz", package = pkg)
dna_file <- system.file("extdata", "pinf_sc50.fasta", package = pkg)
gff_file <- system.file("extdata", "pinf_sc50.gff", package = pkg)

library(vcfR)

vcf <- read.vcfR("alignment/variants.vcf", verbose=FALSE)

dna <- ape::read.dna("alignment/output.fa", format = "fasta")

gff <- read.table("/scratch/prj/sgdp_nanopore/Resources/Homo_sapiens.GRCh38.107.gff3", sep="\t", quote="")

chrom <- create.chromR(name='Supercontig', vcf=vcf, seq=dna, ann=gff)

chrom <- masker(chrom, min_QUAL = 1, min_DP = 300, max_DP = 700, min_MQ = 59.9,  max_MQ = 60.1)

chrom <- proc.chromR(chrom, verbose=TRUE)

chromoqc(chrom, dp.alpha=20)

chromoqc(chrom, xlim=c(5e+05, 6e+05))
