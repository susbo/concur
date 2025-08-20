# CONCUR

CONCUR (**Co**do**n** **c**o**u**nts from **R**ibo-seq) is a tool for calculating codon usage from Ribo-seq data.

## Getting Started

These instructions will allow you to run CONCUR on your computer. CONCUR is a command line tool developed for Linux and macOS.

### Prerequisites

You will need [Perl](https://www.perl.org), [bedtools](http://bedtools.readthedocs.io/en/latest/content/installation.html) and [samtools](https://github.com/samtools/samtools) to be installed in your system.

In addition, you need to have [R](https://cran.r-project.org) and the two R packages pheatmap and RColorBrewer installed to generate some of the figures. If you chose not to install them, you need to run CONCUR with the `--withoutR` parameter.

The following organisms are pre-installed from the Gencode [https://www.gencodegenes.org] project:

* Human - *hg38*, *hg19*
* Mouse - *mm10*, *mm9*

The following organism is pre-installed from Ensembl [http://ensemblgenomes.org/]:

* Rat - *rn9*
* Yeast - *sc3*

If you need to analyze another organism, you can easily do so provided you have a gtf file with annotated protein-coding genes for that organism. Please see additional instructions below.

### Installation

Download the latest release ([v1.1.0](https://github.com/susbo/concur/releases)) from the release tab.

The following commands will unpack CONCUR in your current directory:
```
tar xvfc concur-1.1.0.tar.gz
cd concur-1.1.0
```
Please verify that the tool is working using the example in the `demo` directory. This will only take a minute or so.

Alternatively, clone this repository for the latest version (v.1.1.1-beta).

## Running CONCUR

You need a bam file with read alignments, `alignments.bam`, to run CONCUR. You also need to specify the genome (`-g`) and an output directory (`-o`).

```
perl concur.pl -i alignment.bam -g hg38 -o project_name
```

### Parameters

The following general parameters are available

| Parameter | Description
| :---: | ---
| -i / \-\-input BAM_FILE | Input bam file [*mandatory*]
| -g / \-\-genome GENOME | Genome version (e.g., *hg38*, *hg19*, *mm10*, *mm9*, *rn9* or *sc3*) [*mandatory*]
| -o / \-\-out FOLDER | Output folder name [*mandatory*]
| -n / \-\-name FILENAME| Output file name [*optional*]. Input file name is used by default.
| -w / \-\-withoutR | Run without creating figures using R. This is useful if R is not installed.
| -h / \-\-help | Print help message and quit
| -m / \-\-man | Print help message and quit
| -v / \-\-version | Print version and quit

The following parameters can be used to change some of the default behavour

| Parameter | Description
| :---: | ---
| -s / \-\-size FROM-TO | This will alter the fragment size range included in the analysis (described in section 2.1 of the manuscript). The default range is **20-50**. Non-informative lengths are automatically detected and excluded and the default range should be suitable for most datasets. [*optional*]
| -r / \-\-reads_min READS | This parameter sets the minimum number of reads near the TIS required to include a read set in the analysis (described in section 2.1 of the manuscript). The default threshold is **1000** reads. Increasing this threshold may improve the analysis of deeply sequenced libraries by excluding low-quality read sets that may affect the read set validation steps. [*optional*]
| -f / \-\-filter_outliers THRESHOLD | This option will change the final filtering of the selected read sets. By default, a read set is used in the final codon usage calculations if S_r >= **0.5**\*S_r^max at the P and A site (described in section 2.2.3 of the manuscript). In a dataset where many read sets have passed the validation filters, this threshold would exclude read sets that are nevertheless outliers compared with the best ones. We believe it is generally useful to apply this filter to focus on the most informative reads. However, the threshold can be lowered if keeping as many read sets as possible is of higher importance (use a threshold <0.5), or increased if stricter filtering is desired (use a threshold >0.5). [*optional*]

### Installing Additional Genomes

These instructions will help you to install additional genomes. CONCUR can be run for any organism provided that you have a gtf file containing genes and their coding sequence. If available, CONCUR will use the reading frame information in column 8, otherwise frame can be calculated manually.

#### Using coding sequences fasta
This is the best option if there is a fasta file with the coding sequences available (e.g., Ensembl annotations).

First, download the coding sequence and annotation files (gtf or gff):
```
wget -O Saccer3.cds.fa.gz ftp://ftp.ensemblgenomes.org/pub/fungi/release-40/fasta/saccharomyces_cerevisiae/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz
wget -O Saccer3.gff.gz ftp://ftp.ensemblgenomes.org/pub/fungi/release-40/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.40.gff3.gz

wget -O Ratnor9.cds.fa.gz ftp://ftp.ensembl.org/pub/release-95/fasta/rattus_norvegicus/cds/Rattus_norvegicus.Rnor_6.0.cds.all.fa.gz
wget -O Ratnor9.gtf.gz ftp://ftp.ensembl.org/pub/release-95/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.1.05.gtf.gz
```
Next, run the installation tool. Use `--recalculate` if you wish to disregard the reading frame information in column 8 of the gtf/gff file.
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
Next, run the installation tool. The `--pcg` flag is used to extract the coding sequences from the full transcript sequences. Use `--recalculate` if you wish to disregard the reading frame information in column 8 of the gtf/gff file.
```
perl concur_install_genome.pl --gtf Musmus10.gtf.gz --fasta Musmus10.pcg.fa.gz --short mm10 --pcg
perl concur_install_genome.pl --gtf Homsap38.gtf.gz --fasta Homsap38.pcg.fa.gz --short hg38 --pcg
```
The `--pcg` tag will assume that there is a string in the "CDS:61-1041" format in each fasta header line. The start and end position of the coding sequence is retrieved from this string and is used to extract the coding sequence from the full transcript sequence. Coding sequences where the length is not a multiple of three nucleotides will not be used.

This will create two files for mouse: `data/mm10.bg.txt` and `data/mm10.bed.gz`, and two files for human: `data/hg38.bg.txt` and `data/hg38.bed.gz`.

### Version

The current version is 1.1.0. For other the versions, see the [releases on this repository](https://github.com/susbo/concur/releases). 

### Authors

* **Susanne Bornelöv** - [susbo](https://github.com/susbo)

### License

This project is licensed under the GNU AGPLv3 License - see the [LICENSE.txt](LICENSE.txt) file for details.

### Citation

If you use CONCUR for your work, please cite:

Michaela Frye, Susanne Bornelöv (2020) CONCUR: quick and robust calculation of codon usage from ribosome profiling data, *Bioinformatics*, bta733, [https://doi.org/10.1093/bioinformatics/btaa733](doi:10.1093/bioinformatics/btaa733)
