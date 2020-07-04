# bitou_analysis
This pipeline was used for the population genetic analysis of bitou bush, an invasive plant in Australia.

## run
The `run_stacks.sh` script is primarily for documentation purposes and contains no exception handling. User arguments: 1) A file listing all input fastq files with one file per line, 2) a trimmomatic adapter fasta file, 3) number of cores to use for stacks, 4) the enzyme recognition site to check at the start of each read, and 4) maximum stacks distance (-n, -M).
```
./run_stacks.sh fastq_list.txt adapters.fa 24 CGT 3
```
## clone repository
```
git clone https://github.com/ascheben/bitou_analysis.git
```
## notes
Although an effort was made to make most of this pipeline generic, parts of the pipeline are specific to the analysis it was written for. It should also be noted that while the snakemake pipeline can be executed within reasonable time (~1h) on a personal notebook, the `runs_stacks.sh` script should be run with 24 cores or more on a cluster with substantial memory and walltime, or split into multiple jobs. Finally, the snakemake pipeline was executed within a conda environment with all dependencies installed using conda (except fastStructure and ggtree, which were installed manually and via Bioconductor respectively).
