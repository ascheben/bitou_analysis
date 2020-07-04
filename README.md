# Bitou analysis
This repository documents scripts used in the analysis of bitou bush for the unpublished manuscript 'Successful invasion into coastal ecosystems despite strong genetic bottlenecks in successive introductions of an alien shrub'.

## Outline
Three main analyses are documented here. Aims and results of these analyses are discussed in more detail in the manuscript.
* SNP structure - snakemake pipeline for conducting population structure analyses
* ABC modeling - notebook documentation of code used for modelling population demography
* SNP outlier analysis - notebook documentation for detecting candidate selected SNPs

## Notes
Although an effort was made to make most of this pipeline generic, parts of the pipeline are specific to the analysis it was written for. It should also be noted that while the snakemake pipeline can be executed within reasonable time (~1h) on a personal notebook, the `runs_stacks.sh` script should be run with 24 cores or more on a cluster with substantial memory and walltime, or split into multiple jobs. Finally, the snakemake pipeline was executed within a conda environment with all dependencies installed using conda (except fastStructure and ggtree, which were installed manually and via Bioconductor respectively).
