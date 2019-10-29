# CONCUR

CONCUR (**Co**do**n** **c**o**u**nts from **R**ibo-seq) is a tool for calculating codon usage from Ribo-seq data.

## Getting Started

These instructions will allow you to run CONCUR on your computer. CONCUR is run on the command line on Unix, Linux and Apple OS X.

### Prerequisites

You will need [Perl](https://www.perl.org) and [bedtools](http://bedtools.readthedocs.io/en/latest/content/installation.html) to be installed in your system.

The following organisms are pre-installed from the Gencode [https://www.gencodegenes.org] project:

* Human - *hg38*, *hg19*
* Mouse - *mm10*, *mm9*

The following organism is pre-installed from Ensembl [http://ensemblgenomes.org/]:

* Rat - *rn9*
* Yeast - *sc3*

If you need to analyze another organism, you can easily do so provided you have a gtf file with annotated protein-coding genes for that organism. Please see additional instructions below.

### Installation

Download the latest version ([v1.0](https://github.com/susbo/concur/releases)) from the releases tab.

The following commands will install CONCUR in your current directory:
```
tar xvfc concur-1.0.tar.gz
cd concur-1.0
```
Verify that the tool is working with the example in the demo directory.

## Running CONCUR

You need a bam file with read alignments, `alignments.bam`, to run CONCUR. You also need to specify the genome, and an output directory.

```
perl concur.pl -i alignment.bam -g hg38 -o project_name
```

### Parameters

The following parameters are available

| Parameter | Description
| :---: | ---
| -i \-\-input | Input bam file [*mandatory*]
| -g \-\-genome | Genome version (e.g., *hg38*, *hg19*, *mm10*, *mm9*, *rn9* or *sc3*) [*mandatory*]
| -o \-\-out | Output folder name [*mandatory*]
| -n \-\-name | Output file name [*optional*, input file name is used by default]

### Installing Additional Genomes

These instructions will help you to install additional genomes. CONCUR can be run for any organism provided that you have a gtf file containing genes and their coding sequence. If available, it will use the read frame information in column 8, otherwise frame can be calculated manually.

#### Using coding sequences fasta
This is the best option if there is a fasta file with the coding sequences available (e.g., Ensembl annotations).

First, download the coding sequence and annotation files (gtf or gff):
```
wget -O Saccer3.cds.fa.gz ftp://ftp.ensemblgenomes.org/pub/fungi/release-40/fasta/saccharomyces_cerevisiae/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz
wget -O Saccer3.gff.gz ftp://ftp.ensemblgenomes.org/pub/fungi/release-40/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.40.gff3.gz

wget -O Ratnor9.cds.fa.gz ftp://ftp.ensembl.org/pub/release-95/fasta/rattus_norvegicus/cds/Rattus_norvegicus.Rnor_6.0.cds.all.fa.gz
wget -O Ratnor9.gtf.gz ftp://ftp.ensembl.org/pub/release-95/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.95.gtf.gz
```
Next, run the installation tool. Use `--recalculate` if you wish to disregard the read frame information in column 8 of the gtf/gff file.
```
perl concur_install_genome.pl --gtf Saccer3.gff.gz --fasta Saccer3.cds.fa.gz --short sc3
perl concur_install_genome.pl --gtf Ratnor9.gtf.gz --fasta Ratnor9.cds.fa.gz --short rn9
```
This will create two files for yeast: data/sc3.bg.txt and data/sc3.bed.gz, and two files for rat: data/rn9.bg.txt and data/rn9.bed.gz.

#### Using protein coding gene fasta
This is the best option if there is not a separate fasta file with only the coding sequences available, but there is a fasta file with transcript sequences and information about the CDS position (e.g., Gencode annotations).

First, download the coding sequence and annotation files (gtf or gff):
```
wget -O Musmus10.pcg.fa.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M20/gencode.vM20.pc_transcripts.fa.gz
wget -O Musmus10.gtf.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M20/gencode.vM20.primary_assembly.annotation.gtf.gz

wget -O Homsap38.pcg.fa.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.pc_transcripts.fa.gz
wget -O Homsap38.gtf.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.primary_assembly.annotation.gtf.gz
```
Next, run the installation tool. The `--pcg` flag is used to extract the coding sequences from the full transcript sequences. Use `--recalculate` if you wish to disregard the read frame information in column 8 of the gtf/gff file.
```
perl concur_install_genome.pl --gtf Musmus10.gtf.gz --fasta Musmus10.pcg.fa.gz --short mm10 --pcg
perl concur_install_genome.pl --gtf Homsap38.gtf.gz --fasta Homsap38.pcg.fa.gz --short hg38 --pcg
```
The `--pcg` tag will assume that there is a string in the "CDS:61-1041" format in each fasta header line. The start and end position of the coding sequence is retrieved from this string and is used to extract the coding sequence from the full transcript sequence. Coding sequences where the length is not a multiple of three nucleotides will not be used.

This will create two files for mouse: data/mm10.bg.txt and data/mm10.bed.gz, and two files for human: data/hg38.bg.txt and data/hg38.bed.gz.

### Version

The current version is 1.0. For other the versions, see the [tags on this repository](https://github.com/your/project/tags). 

### Authors

* **Susanne Bornelöv** - [susbo](https://github.com/susbo)

### License

This project is licensed under the MIT License - see the [LICENSE.txt](LICENSE.txt) file for details.

### Citation

At the moment you can refer to the github repository [https://github.com/susbo/concur]. We are in the process of submitting a paper describing the tool.
