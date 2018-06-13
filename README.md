# Pitagora Workflows in CWL

This repository is to share the [Common Workflow Language](https://www.commonwl.org) definition files of tools and workflows, maintained by Pitagora Network, as known as Galaxy Community Japan.

## Tools

- Get Data
  - [download-sra](tools/download-sra)
  - [fastq-dump](tools/fastq-dump)
  - [pfastq-dump](tools/pfastq-dump)
- Data Formatting
  - [samtools](tools/samtools)
- Read Alignment
  - [bwa](tools/bwa)
  - [bowtie](tools/bowtie)
  - [bowtie2](tools/bowtie2)
  - [last](tools/last)
  - [star](tools/star)
  - [HISAT2](tools/hisat2)
- Variant Analysis
  - [gatk](tools/gatk)
  - [annovar](tools/annovar)
  - [picard](tools/picard)
- Transcriptome Analysis
  - splice junction mapper
    - [tophat2](tools/tophat2)
  - transcript assemble / quantification
    - [cufflinks](tools/cufflinks)
    - [stringtie](tools/stringtie)
    - [eXpress](tools/eXpress)
    - [rsem](tools/rsem)
    - [rsubread](tools/rsubread)
    - [kallisto](tools/kallisto)
    - [sailfish](tools/sailfish)
    - [salmon](tools/salmon)
  - de novo transcriptome assemble
    - [trinity](tools/trinity)

## Workflows

- RNA-Seq
  - HISAT2-*
    - [hisat2-cufflinks](workflows/hisat2-cufflinks)
    - [hisat2-stringtie](workflows/hisat2-stringtie)
    - [hisat2-eXpress](workflows/hisat2-eXpress)
    - [hisat2-rsem](workflows/hisat2-rsem)
  - STAR-*
    - [star-cufflinks](workflows/star-cufflinks)
    - [star-stringtie](workflows/star-stringtie)
    - [star-eXpress](workflows/star-eXpress)
    - [star-rsem](workflows/star-rsem)
  - [kallisto](workflows/kallisto)
  - [sailfish](workflows/sailfish)
  - [salmon](workflows/salmon)
