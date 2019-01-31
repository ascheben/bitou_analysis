# bitou_analysis
Pipeline used for population genetic analysis of bitou bush
# dependencies
## core tools
* Snakemake 5.4.0
* conda 4.6.1
* Python 3.5
## sequence analysis
* trimmomatic
* fastqc
* multiqc 1.0
* stacks 2.1
## variant analysis
* vcftools 0.1.16
* SNPrelate 1.14.0
* gdsfmt 1.16.0
* fastStructure
* plink 1.90
* raxml 8.2.12
* vcf2phylip.py 1.5 (https://github.com/edgardomortiz/vcf2phylip)
## data visualisation
* ggtree 1.14.6
* pophelper 2.2.7
* phangorn 2.4.0

# run
```
./run_stacks.sh
snakemake
```
