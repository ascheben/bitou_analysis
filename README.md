# bitou_analysis
Pipeline used for population genetic analysis of bitou bush
# dependencies
## core tools
* Snakemake 5.4.0
* conda 4.6.1
* Python 3.5
## sequence analysis
* trimmomatic 0.36
* fastqc 0.11.4
* multiqc 1.0
* stacks 2.1
## variant analysis
* vcftools 0.1.16
* SNPrelate 1.14.0
* gdsfmt 1.16.0
* fastStructure 1.0
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
```
./run_stacks.sh
snakemake
```
# notes
Although an effort was made to make most of this pipeline generic, parts of the pipeline are specific to the analysis it was written for. It should also be noted that while the snakemake pipeline can be executed within reasonable time (~1h) on a personal notebook, the stacks pipeline for read processing and de novo assembly should be run with the suggested 24 cores or more on a cluster with substantial memory and walltime available. Finally, the snakemake pipeline was executed within a conda environment with all dependencies installed using conda (except fastStructure and ggtree, which were installed manually and via Bioconductor respectively).
