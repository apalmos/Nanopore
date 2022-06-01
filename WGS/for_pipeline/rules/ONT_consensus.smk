# create environment
rule setup_environment:
  output:
    directory("./assembly")
  shell:
    "mkdir QC; mkdir assembly; mkdir alignment; mkdir error_correction; mkdir consensus; mkdir stitch"

# download dependencies environment
rule install_dependencies
    output:
        touch("resources/software/install_minimap2.done")
    conda:
        "../envs/ONT_consensus.yaml"
    shell:
        "conda install -c bioconda minimap2; conda install -c bioconda miniasm; conda install -c bioconda minipolish; conda install -c bioconda NanoStat; conda install -c bioconda sniffles; conda install -c bioconda longstitch;


# de-novo alignment
rule run_minimap2:
    input:
        "QC/{sample}.fastq.gz"
    output:
        "./assembly/{sample}.paf"
    params:
        dir="result"
    shell:
        "minimap2 -X -t 10 -x ava-ont ./QC/{sample}.fastq.gz ./QC/{sample}.fastq.gz > ./assembly/{sample}.paf"

# file wrangling
rule wrangle_paf:
    input:
        "./assembly/{sample}.paf"
    output:
        "./assembly/{sample}self.paf"
    params:
        dir="result"
    shell:
        "awk '{if($1==$6){print}}' ./assembly/{sample}.paf > ./assembly/{sample}self.paf"

# assembly
rule assemble:
    input:
        "./assembly/{sample}.paf"
    output:
        "./assembly/{sample}self.paf"
    params:
        dir="result"
    shell:
        "miniasm -f ./QC/nanofilt.fa ./assembly/AONT.paf > ${directory}/assembly/AONT.gfa
