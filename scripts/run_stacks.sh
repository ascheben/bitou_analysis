#!/usr/bin/env bash

# This script will:
# 1) clean and QC single end Illumina reads from ddRAD libraries
# 2) Run the manual stacks pipeline 
# Note that fastq.gz files, an adapters file and a popmap.txt file must be in the working directory

# Use example
# ./run_stacks.sh fastq_list.txt adapters.fa 24 CGT 3

# Working directory
WD=$(pwd)
# List of gzipped fastq files
FASTQ_LIST="$1"
# Fasta format trimmomatic adapters file
ADAPTERS="$2"
# Threads to use for stacks analyses
THREADS="$3"
# Enzyme recognition site
SITE="$4"
# Distance between and within loci for stack building (-n,-M)
DIST="$5"


# Create output directories
mkdir clean stacks raw trim

while read file
do
	# Drop all reads with adapter contamination
    trimmomatic SE ${file} ${file%%.fq.gz}.trim.fq.gz ILLUMINACLIP:${ADAPTERS}:2:30:10 MINLEN:95
	#exclude reads not beginning with the restriction enzyme recognition site
	zcat ${file%%.fq.gz}.trim.fq.gz | \
	grep -P "^${SITE}[ACGTNacgtn]*$" -B 1 -A 2 | \
	grep -v -P "^--$" | \
	gzip -c > ${file%%.fq.gz}.${SITE}.trim.fq.gz
   mv ${file%%.fq.gz}.trim.fq.gz ${WD}/trim/
	# Run quality control program on reads
	fastqc ${file%%.fq.gz}.${SITE}.trim.fq.gz
	mv ${file%%.fq.gz}.${SITE}.trim.fq.gz ${WD}/clean/${file}
   mv ${file} ${WD}/raw/
done<${FASTQ_LIST}
# Summarize individual quality control statistics
multiqc .
rm *_fastqc.html *_fastqc.zip
# A fastqc v0.1.11 and MultiQC v1.0 analysis were used for quality control of the final clean reads. An awk script was used to sum rows in the multiqc output, obtaining average base phred score per sample. The mean per base quality across all samples was 36.424.

./scripts/mean_rows.sh ${WD}/multiqc_data/mqc_fastqc_per_base_sequence_quality_plot_1.txt > ${WD}/multiqc_data/mean_per_base_sequence_quality.txt

# Manually run stacks denovo_map pipeline
# Commands modified from http://catchenlab.life.illinois.edu/stacks
id=1

INFASTQ=$(find ${WD}/clean -name '*.fq.gz' | sed 's#.*/##')
echo $INFASTQ
echo $INFASTQ | tr ' ' '\n' | while read infile
do
    ustacks -f ${WD}/clean/$infile -o ${WD}/stacks -i $id --name ${infile%%.fq.gz} -M $DIST -p $THREADS
    let "id+=1"
done
#
# Build catalogue
cstacks -n $DIST -P ${WD}/stacks -M ${WD}/popmap.txt -p $THREADS
sstacks -P ${WD}/stacks -M ${WD}/popmap.txt -p $THREADS
tsv2bam -P ${WD}/stacks -M ${WD}/popmap.txt -t $THREADS
## Genotype individuals
gstacks -P ${WD}/stacks -M ${WD}/popmap.txt -t $THREADS
## Output SNPs in VCF format
populations -P ${WD}/stacks -M ${WD}/popmap.txt -t $THREADS -p 1 --vcf
