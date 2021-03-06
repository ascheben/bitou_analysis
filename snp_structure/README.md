# SNP structure analysis
This pipeline was used for the exploratory population structure analysis of bitou bush, an invasive plant in Australia.

# Input files

The required inputs include a VCF, a tab-separated population map file (format: "sample\tpopulation") and a nexus output file from SplitsTree. Sample names must be consistent between these files. Inputs are provided by modifying the `config.yaml` file. The below lines should be modified.

```
invcf: bitou_core_raw
popmap: popmap_core.txt
nnet: bitou_m3_ibs.dst.nex
```

The snakemake pipeline can be executed after user configuration has been added to the `config.yaml` file.
```
snakemake
```

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
