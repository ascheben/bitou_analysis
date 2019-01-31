# bitou_analysis
This pipeline was used for the population genetic analysis of bitou bush, an invasive plant 
# dependencies
## core tools
* snakemake 5.4.0
* conda 4.6.1
* python 3.5
## sequence analysis
* trimmomatic 0.36
* fastqc 0.11.4
* multiqc 1.0
* stacks 2.1
## variant analysis
* vcftools 0.1.16
* snprelate 1.14.0
* gdsfmt 1.16.0
* faststructure 1.0
* plink 1.90
* raxml 8.2.12
* vcf2phylip.py 1.5 (https://github.com/edgardomortiz/vcf2phylip)
## data visualisation
* ggtree 1.14.6
* pophelper 2.2.7
* phangorn 2.4.0
# clone repository
```
git clone https://github.com/ascheben/bitou_analysis.git
```
# run
To run the sequence analysis, gzipped fastq files and a stacks popmap file should be in the working directory. User arguments: 1) A file listing all input fastq files with one file per line, 2) a trimmomatic adapter fasta file, 3) number of cores to use for stacks, and 4) the enzyme recognition site to check at the start of each read.
```
./run_stacks.sh fastq_list.txt adapters.fa 24 CGT
```
The snakemake pipeline can be executed after user configuration has been added to the `config.yaml` file.
```
snakemake
```
# notes
Although an effort was made to make most of this pipeline generic, parts of the pipeline are specific to the analysis it was written for. For instance, the enzyme recognition step in `run_stacks.sh` is only required when an enzyme is not available in stacks. It should also be noted that while the snakemake pipeline can be executed within reasonable time (~1h) on a personal notebook, the stacks pipeline for read processing and de novo assembly should be run with the suggested 24 cores or more on a cluster with substantial memory and walltime available. Finally, the snakemake pipeline was executed within a conda environment with all dependencies installed using conda (except fastStructure and ggtree, which were installed manually and via Bioconductor respectively).
