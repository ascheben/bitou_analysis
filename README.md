# Bitou analysis
This repository documents scripts used in the analysis of bitou bush for the manuscript 'Genomics reveals the history of a complex plant invasion and improves the management of a biological invasion from the South African-Australian biotic exchange'.

## Outline
Three main analyses are documented here. Aims and results of these analyses are discussed in more detail in the manuscript.
* SNP structure - snakemake pipeline for conducting population structure analyses
* ABC modeling - notebook documentation of code used for modelling population demography
* SNP outlier analysis - notebook documentation for detecting candidate selected SNPs

## Notes
Although an effort was made to make most of this pipeline generic, parts of the pipeline are specific to the analysis it was written for. It should also be noted that while the snakemake pipeline can be executed within reasonable time (~1h) on a personal notebook, the `runs_stacks.sh` script should be run with 24 cores or more on a cluster with substantial memory and walltime, or split into multiple jobs. Finally, the snakemake pipeline was executed within a conda environment with all dependencies installed using conda (except fastStructure and ggtree, which were installed manually and via Bioconductor respectively).

## Citation
Byrne D, Scheben A, Scott JK, Webber BL, Batchelor KL, Severn-Ellis AA, Gooden B, Bell KL. Genomics reveals the history of a complex plant invasion and improves the management of a biological invasion from the South African-Australian biotic exchange. Ecol Evol. 2022 Aug 23;12(8):e9179. doi: 10.1002/ece3.9179. PMID: 36016815; PMCID: PMC9396708.
