Bitou outlier analysis
================

## Outlier analysis

Convert VCF to BAYENV/BAYESCAN using the java program PGDSpider2.

``` bash
java -Xmx1024m -Xms512m -jar PGDSpider_2.1.1.5/PGDSpider2-cli.jar -inputfile bitou_core_snps_nopisi.vcf -inputformat VCF -outputfile bitou_core_snps_nopisi.bayescan -outputformat GESTE_BAYE_SCAN -spid template_VCF_GESTE_BAYE_SCAN.spid
```

When you first attempt to execute the above command with a spid file, it
will generate a template like the one below. I filled it out in the
following way. To successfully convert run, the above command with a
template like below as input.

``` bash
 spid-file generated: Fri Aug 16 18:49:47 AWST 2019

# VCF Parser questions
PARSER_FORMAT=VCF

# Do you want to include a file with population definitions?
VCF_PARSER_POP_QUESTION=true
# Only input following regions (refSeqName:start:end, multiple regions: whitespace separated):
VCF_PARSER_REGION_QUESTION=
# What is the ploidy of the data?
VCF_PARSER_PLOIDY_QUESTION=2
# Only output following individuals (ind1, ind2, ind4, ...):
VCF_PARSER_IND_QUESTION=
# Output genotypes as missing if the read depth of a position for the sample is below:
VCF_PARSER_READ_QUESTION=1
# Take most likely genotype if "PL" or "GL" is given in the genotype field?
VCF_PARSER_PL_QUESTION=false
# Do you want to exclude loci with only missing data?
VCF_PARSER_EXC_MISSING_LOCI_QUESTION=true
# Select population definition file:
VCF_PARSER_POP_FILE_QUESTION=final_pop_all.txt
# Only output SNPs with a phred-scaled quality of at least:
VCF_PARSER_QUAL_QUESTION=1
# Do you want to include non-polymorphic SNPs?
VCF_PARSER_MONOMORPHIC_QUESTION=true
# Output genotypes as missing if the phred-scale genotype quality is below:
VCF_PARSER_GTQUAL_QUESTION=1

# BAYENV Writer questions
WRITER_FORMAT=BAYENV

# Do you want to save an additional sample file with sample sizes?
BAYENV_WRITER_WRITE_SAMPLE_FILE_QUESTION=true
# Save sample file
BAYENV_WRITER_SAMPLE_FILE_QUESTION=true
# Save sample/loci names file
BAYENV_WRITER_INFO_FILE_QUESTION=true
# Do you want to save two additional files with used sample and loci names?
BAYENV_WRITER_WRITE_INFO_FILE_QUESTION=true
# Assign half missing genotypes (one allele missing) as complete missing?
BAYENV_WRITER_HALF_MISSING_QUESTION=false
```

Generate coverance matrices, e.g.:

``` bash
./bayenv2 -i bitou_core_snps_nopisi.bayenv2 -k 200000 -r 5432 -p 5 > matrix1.out
```

After generation of five replicate matrices for bayenv using the bitou
SNP set of 16, we extract the final matrices from each of the replicate
analyses.

``` bash
cd /mnt/d/Documents/UWA/Applied\ Bioinformatics/Bitou/DIYABC_2019/bayenv
for i in {0..4}; do 
  cat matrix$i.out | \
  perl -lne 'if(/ITER = 200000/){$ok=1}elsif($ok){ print }' \
  | sed '/^$/d' > final_matrix_$i.txt
done
```

Calculate the mean of the five replicate matrices. We then scale the
covariance matrix to visualise a correlation matrix. The results can be
assessed to ensure that they are supported by the structure analysis.

``` r
library(corrplot)
setwd("D://Documents//UWA//Applied Bioinformatics//Bitou//DIYABC_2019//bayenv")

# read all final matrices
m1 = as.matrix( read.table(file="final_matrix_0.txt", header=F) )
m2 = as.matrix( read.table(file="final_matrix_1.txt", header=F) )
m3 = as.matrix( read.table(file="final_matrix_2.txt", header=F) )
m4 = as.matrix( read.table(file="final_matrix_3.txt", header=F) )
m5 = as.matrix( read.table(file="final_matrix_4.txt", header=F) )


# make a list of matrices and get mean as explained in:
mat_list = list( m1, m2, m3, m4, m5 )
mean_mat = apply(simplify2array(mat_list), c(1,2), mean)

covmat <- cov(mean_mat)
cormat <- cov2cor(covmat)
#V1 : EAU
#V2: EBEACH
#V3: DURBAN
#V4: STJOHN
#V5: WAU
corrplot(cormat, method="number",number.digits = 3)

# write resulting mean cov matrix
write.table(mean_mat,file="final_matrix_mean.txt",
            sep="\t",row.names=F,col.names=F,quote=F)
```

The above shows that the covariance matrix is in line with the structure
analysis, PCA and phylogeny. We now compare the bayescan and bayenv
outcomes, retaining only the significant loci and finding shared loci
between the methods.

After conversion from VCF, Bayescan can be run quite easily using the
binary provided (or compiled yourself from source). Here the pr\_odds
flag is important and indicates how likely it is for loci to deviate
from neutral. A pr\_odds of 100 means that the chances are 100 to 1,
which is a conservative estimate and sensible when \>1k loci are being
tested.

``` bash
bayescan/BayeScan2.1_linux64bits bitou_core_snps_nopisi.bayescan -od prob100_rep1 -pr_odds 100 -threads 12"
```

The above was run for three replicates. The bayescan “fst.txt” output
file is key for finding the loci under selection. Below I merge the
bayenv and bayescan results to get a list of shared loci considered to
be under selection. I also test the three bayescan replicate runs for
convergence.

``` r
library(coda)
```

    ## Warning: package 'coda' was built under R version 3.5.3

``` r
library(dplyr)
```

    ## Warning: package 'dplyr' was built under R version 3.5.3

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
myxtx <- read.csv("D://Documents//UWA//Applied Bioinformatics//Bitou//DIYABC_2019//bayenv//all.xtx", sep = "\t", header = FALSE)
# Dedupe as some values seem to be erroneously calculated multiple times
myxtxu <- distinct(myxtx,V1, .keep_all= TRUE)
# Calculate 99% quantile to use as conservative cut-off
quantile(myxtxu$V2, probs = c(0.01, 0.99))
```

    ##        1%       99% 
    ##  2.547904 17.676220

``` r
# Retain loci in this quantile
canloci <- myxtxu[which(myxtxu[,2]>17.676220),]
# Modify the locus names to allow merging with bayescan loci
canloci$V1 <-  gsub("batch/snp_batch","",canloci$V1) 
canloci$V1 <- gsub("(?<![0-9])0+", "", canloci$V1, perl = TRUE)
canloci$V1 <- as.numeric(canloci$V1)
canloci$V1 <-  canloci$V1 + 1
names(canloci) <- c("Locus_num", "xtx")
```

``` r
bsfst <- read.table("D://Documents//UWA//Applied Bioinformatics//Bitou//DIYABC_2019//bayescan//bitou_core_snps_nopisi.baye_fst_rep1.txt",header=TRUE)
# Use 5% FDR cut-off to retain significant loci
bscanloci <- bsfst[which(bsfst[,3]<0.05),]
# Add column with row numbers 
bscanloci <- tibble::rownames_to_column(bscanloci, "Locus_num")
# Make locus numbers numeric
bscanloci$Locus_num <- as.numeric(bscanloci$Locus_num)
# Find shared loci
shared <- inner_join(bscanloci, canloci, by = "Locus_num")

### BAYESCAN CONVERGENCE TESTING ###

chain1<-read.table("D://Documents//UWA//Applied Bioinformatics//Bitou//DIYABC_2019//bayescan//bitou_core_snps_nopisi.baye_rep1.sel",header=TRUE)
chain2<-read.table("D://Documents//UWA//Applied Bioinformatics//Bitou//DIYABC_2019//bayescan//bitou_core_snps_nopisi.baye_rep2.sel",header=TRUE)
chain3<-read.table("D://Documents//UWA//Applied Bioinformatics//Bitou//DIYABC_2019//bayescan//bitou_core_snps_nopisi.baye_rep3.sel",header=TRUE)
# Remove first column
chain1<-chain1[-c(1)]
chain2<-chain2[-c(1)]
chain3<-chain3[-c(1)]
#Create a MCMC objects with the correct thinning interval
chain1<-mcmc(chain1,thin=10)
chain2<-mcmc(chain2,thin=10)
chain3<-mcmc(chain3,thin=10)
#Carry out Geweke diagnostic test for each chain
geweke.diag(chain1, frac1=0.1, frac2=0.5)
```

    ## 
    ## Fraction in 1st window = 0.1
    ## Fraction in 2nd window = 0.5 
    ## 
    ##     Fst1     Fst2     Fst3     Fst4     Fst5 
    ##  1.05436  0.04334 -1.47040  0.16835 -0.33703

``` r
geweke.diag(chain2, frac1=0.1, frac2=0.5)
```

    ## 
    ## Fraction in 1st window = 0.1
    ## Fraction in 2nd window = 0.5 
    ## 
    ##    Fst1    Fst2    Fst3    Fst4    Fst5 
    ##  1.4697 -0.6552  1.2129 -1.7263 -1.0293

``` r
geweke.diag(chain3, frac1=0.1, frac2=0.5)
```

    ## 
    ## Fraction in 1st window = 0.1
    ## Fraction in 2nd window = 0.5 
    ## 
    ##    Fst1    Fst2    Fst3    Fst4    Fst5 
    ## -1.2517 -0.4891  0.4423 -0.3204  0.4592

``` r
# A plot can also be reported
#geweke.plot(chain1, frac1 = 0.1, frac2 = 0.5, nbins = 20,
#            pvalue = 0.05, auto.layout = TRUE)
# When two chains are available a Gelman plot can be used to test for convergence
#combined = mcmc.list(chain1,chain2)
#plot(combined)
# A factor of 1 means that between- and within-chain variance are the same
# The rule of thumb is that values below 1.1 or so are OK
#gelman.diag(combined)
```

The Geweke tests for the three chains shown above suggest that the runs
have converged.

The 21 SNPs were gathered from the VCF.

``` bash
#The list of outlier loci numbers was pasted into a text file
#The numbers were then extracted from the original input vcf file
grep '#' bitou_core_snps_nopisi.vcf > outlier_loci.vcf
while read l; do grep -v '#' bitou_core_snps_nopisi.vcf | head -${l} | tail -1 >> outlier_loci.vcf.tail;done<outlier_loci.txt
```

Inspecting these 21 SNPs that were detected as outliers in TASSEL, we
can see that these SNPs are all show and AUS vs RSA divide, with the
EBEACH group revealing affinities to the AUS samples.
