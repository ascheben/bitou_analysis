Bitou ABC analysis
================

This document outlines the different steps in the ABC random forest and
outlier analyses carried out on bitou bush ddRAD SNPs in 2019.

## Convert VCF to DIYABC input file

Starting with the full bitou SNP set, the first step in the ABC analysis
is to prepare a SNP matrix as input for DIYABC. It is worth noting that
the order of populations in the input SNP matrix, is the order that is
used in the ABC analysis\! So “sample 1” in the ABC parameters
corresponds to the first population in the SNP input file.

``` bash
# Enter working directory containing VCF file
cd "/mnt/d/Documents/UWA/Applied Bioinformatics/Bitou/DIYABC_2019/final_input_snps"
# To reduce the size of the SNP set, we can remove all SNPs with any missing data
vcftools --vcf bitou_core_snps_nopisi.vcf --max-missing 1 --recode --stdout > bitou_core_snps_nopisi_nomiss.vcf
# We can then check how many SNPs are left
#grep -v '#' bitou_core_snps_nopisi_nomiss.vcf | wc -l
# Alternatively we can retain SNPs with up to 95% missing genotypes
##vcftools --vcf bitou_core_snps_nopisi.vcf --max-missing 0.95 --recode --stdout > bitou_core_snps_nopisi_mm95.vcf
```

Next we should ensure that each population of the five defined in our
popmap file is represented. Stacks populations is installed in a conda
environment, so we need to load this first. We can skip this step since
we have no missing data in our dataset.

``` bash
# Enter working directory containing VCF file
cd "/mnt/d/Documents/UWA/Applied Bioinformatics/Bitou/DIYABC_2019/final_input_snps"
conda activate runstacks
populations -V bitou_core_snps_nopisi_nomiss.vcf -O . -M final_popmap.txt -p 5 --vcf --no_hap_exports
conda deactivate
```

Now we need to remove monomorphic SNPs. These should have been removed
with the MAF filter earlier, however, this does not exclude monomorphic
heterozygous genotypes. A python script can filter out any monomorphic
SNPs.

``` python
import sys
import csv

vcf = sys.argv[1] #input vcf file

with open(vcf) as tsv:
    for line in csv.reader(tsv, dialect="excel-tab"):
        # skip comment lines
        if not line[0].startswith("#"):
            # number of non-genotype fields in standard VCF format
            num_info_fields = 9
        # get number of samples
            num_samples = len(line) - num_info_fields
            het_count = 0
            alt_count = 0
            hom_count = 0
            null_count = 0
        # loop through alleles for each SNP and count genotypes
            for allele in line[num_info_fields:]:
                g = allele.split(":")[0]
                if g == "0/1":
                    het_count = het_count + 1
                elif g == "1/1":
                    alt_count = alt_count + 1
                elif g == "0/0":
                    hom_count = hom_count + 1
                elif g == "./.":
                    null_count = null_count + 1

            # print SNP if user-defined thresholds are met
            tot = null_count + het_count
            tot2 = null_count + hom_count
            tot3 = null_count  + alt_count
            if tot == num_samples or tot2 == num_samples or tot3 == num_samples:
                pass
            else:
                print(*line, sep = '\t')
    # print all comment lines
        else:
            print(*line, sep = '\t')
```

We can use the script like so:

``` bash
# Enter working directory containing VCF file
cd "/mnt/d/Documents/UWA/Applied Bioinformatics/Bitou/DIYABC_2019/final_input_snps"
python3 filtmono.py bitou_core_snps_nopisi_nomiss.p.vcf > bitou_core_snps_nopisi_nomiss_nomono.p.vcf
```

If there are now too many SNPs, we can downsample to a chosen number of
SNPs, e.g. 4000

``` bash
# Enter working directory containing VCF file
cd "/mnt/d/Documents/UWA/Applied\ Bioinformatics/Bitou/DIYABC_2019/final_input_snps"
 (grep '#' bitou_core_snps_nopisi_nomiss_nomono.p.vcf && grep -v '#'bitou_core_snps_nopisi_nomiss_nomono.p.vcf | shuf -n 4000) >bitou_core_snps_nopisi_nomiss_nomono_4k.p.vcf
```

Finally we can convert to DIYABC SNP format. I have modified the
vcf2diyabc.py (<https://github.com/loire/vcf2DIYABC.py>) script to
directly take user arguments.

``` bash
# Enter working directory containing VCF file
cd "/mnt/d/Documents/UWA/Applied Bioinformatics/Bitou/DIYABC_2019/final_input_snps"
# To convert to DIYABC SNP format we need a popmap file with a dummy sex column, for example:
echo "DIYABC popmap file looks like this:"
head final_pop_diyabc_all.txt
printf "\n\n"
# We can now conver to SNP format using vcf2diyabc.py
python2 vcf2DIYABC_clean.py bitou_core_snps_nopisi_nomiss_nomono_4k.p.vcf final_pop_diyabc_all.txt
# The resulting SNP file can then be sorted by population
 (head -n 3 bitou_core_snps_nopisi_nomiss_nomono_4k.p.DIYABC.snp && tail -n +4 bitou_core_snps_nopisi_nomiss_nomono_4k.p.DIYABC.snp| sort) > bitou_core_snps_nopisi_nomiss_nomono_4k.p.DIYABC.sort.snp
# We shouldn't forget to add the sex ratio to the SNP, as this leads to a loading error in DIYABC
sed -i '1s/insert title here <NM=$insert_sex_ratio$NF>/bitou <NM=1.0NF>/' bitou_core_snps_nopisi_nomiss_nomono_4k.p.DIYABC.sort.snp
# The sex ratio is used for non-autosomal variants in animals and will not be used with our SNPs
# We therefore use a dummy value of 1.0
```

    ## DIYABC popmap file looks like this:
    ## EAU_NSW_WOL_BB_02    9   EAU
    ## EAU_NSW_WOL_BB_03    9   EAU
    ## EAU_NSW_DUN_BB_02    9   EAU
    ## EAU_NSW_DUN_BB_03    9   EAU
    ## EAU_NSW_ILUKA_BB_02  9   EAU
    ## EAU_NSW_ILUKA_BB_03  9   EAU
    ## EAU_NSW_MINNIE_BB_01 9   EAU
    ## EAU_NSW_MINNIE_BB_02 9   EAU
    ## EAU_NSW_WOL_BB_04    9   EAU
    ## EAU_NSW_DUN_BB_01    9   EAU
    ## 
    ## 
    ## Parsing bitou_core_snps_nopisi_nomiss_nomono_4k.p.vcf...
    ## Writing outputfile as: bitou_core_snps_nopisi_nomiss_nomono_4k.p.DIYABC.snp

## DIYABC parameters - using header files

The header file specifies all model parameters and scenarios to be
tested. Together with the SNP input file, it makes up the core input for
DIYABC analysis. There are two data dryad repositories for published
data that provide excellent examples of a range of SNP input data and
ABC paramter header files: [Extraordinarily rapid speciation in a marine
fish](https://datadryad.org/resource/doi:10.5061/dryad.f6154) and
[Molecular data reveal recent genetic diversification in a
hummingbird](https://datadryad.org/resource/doi:10.5061/dryad.88589).
These header files can be generated via text editor, or automatically
using the DIYABC user interface.

Below is an example of a scenario that was tested in bitou. Note that in
the first line we have parameters for the population size of 7
populations, although we only sampled 5 populations. The populations 6
(RSA\_NA) and 7 (RSA) represent an unsampled South African population
and an unsampled South African ancestral population respectively. The
ancestral population is required because we do not know how the four
other South African populations are related.

``` bash
EAU DUR EBEACH STJOHN WAU RSA_NA RSA
0 sample 1 # population 1 (first occuring in SNP input file) sampled at time 0
0 sample 2 # population 2 (second occuring in SNP input file) sampled at time 0
0 sample 3 # population 3 (third occuring in SNP input file) sampled at time 0
0 sample 4 # population 4 (fourth occuring in SNP input file) sampled at time 0
0 sample 5 # population 5 (fifth occuring in SNP input file) sampled at time 0
t_invW-db2 VarNe 5 WAU0 # bottleneck after colonization of WAU
t_adm split 5 6 1 ra # admixture between EAU population and South African ghost population
t_invE-db1 VarNe 1 EAU0 # bottleneck after colonization of EAU
t_invE merge 6 1 # South African ghost population and EAU population diverge at the time of invasion t_invE
t_anc merge 7 2
t_anc merge 7 3
t_anc merge 7 4
t_anc merge 7 6

# three priors for effective population size of four populations
0 sample 1 
0 sample 2
0 sample 3 
```

The priors set for the above scenario are as follows.

``` bash
# Below are the population sizes uniform distributions
EAU N UN[10.0,10000.0,0.0,0.0] # Uniform distribution of 10 diploid individuals to 10k diploid individuals
DUR N UN[10.0,10000.0,0.0,0.0]
EBEACH N UN[10.0,10000.0,0.0,0.0]
STJOHN N UN[10.0,10000.0,0.0,0.0]
WAU N UN[10.0,10000.0,0.0,0.0]
RSA_NA N UN[10.0,10000.0,0.0,0.0]
RSA N UN[10.0,10000.0,0.0,0.0]

# Below here are the time parameters (in generations) and additional bottleneck population sizes
t_invW T UN[8.0,10.0,0.0,0.0] # time of invasion of Western Australia in generations
db2 T UN[1.0,5.0,0.0,0.0] # WAU bottleneck duration is fixed to 1-5 generations
WAU0 N UN[5.0,1000.0,0.0,0.0] # WAU bottleneck population size
t_adm T UN[9.0,44.0,0.0,0.0] # time of admixture between EAU and RSA_NA pop
ra A UN[0.001,0.5,0.0,0.0] # admixture rate is bounded to maximum 0.5, with 1.0 being maximum admixture
t_invE T UN[38.0,45.0,0.0,0.0] # time of invasion of Eastern Australia in generations
db1 T UN[1.0,5.0,0.0,0.0] # EAU bottleneck duration is fixed to 1-5 generations
EAU0 N UN[5.0,1000.0,0.0,0.0] # EAU bottleneck population size
t_anc T UN[10.0,10000.0,0.0,0.0] # The time all the populations split from the ancestral population

# Time constraints
t_anc>t_invW
t_anc>t_invE
t_invE>t_invW
t_anc>t_adm
t_invE>t_adm
t_adm>t_invW
```

## Running DIYABC on a cluster

Based on my experience, I do not recommend running DIYABC on a Linux
cluster. Many analyses resulted in segmentation faults, which are hard
to trace, e.g. for model checking and parameter estimation.
Additionally, I was not able to import the reftables generated via the
command line into the DIYABC GUI program to carry out these analyses.

To run the analysis on a cluster, the DIYABC linux binary is required,
which is easy to compile. In the GUI, under ‘settings’ it is possible to
select cluster, which when selected will generate an archive of useful
files when the user hits the main ‘analyse’ button to start the
simulations. The key files are the node.sh and the launch.sh files. The
launch.sh can be used to generate RNG files required for the
parallelized analysis. The node.sh script can launch the analysis. These
scripts did not work out of the box on my server, so I modified them.
The launch script should be modified to generate the desired number of
RNG files (eg 200 for 200 parallel jobs). These files should then be
moved into separate directores named “mtymp\_$RNGNUM”. I modified the
node.sh script to use this name for temporary working directories in
which reftables are then generated. The node.sh script is then launched
in parallel for each RNG file.

``` bash
DIYABCEXE="/scratch/pawsey0149/ascheben/bitou/diyabc3/diyabc_core-2.1.0-linux-x64"
CPU=1
# Random number generator file
RNGNUM=1
# Number of generations to run simulations
NGEN=1000000
# Working dir
WD="/scratch/pawsey0149/ascheben/bitou/diyabc3"
# Input SNP file
INSNP="/scratch/pawsey0149/ascheben/bitou/diyabc3/bitou_core_snps_nopisi_nomiss.p.snps.DIYABC.snp"

bash node.sh '$DIYABCEXE' $CPU $NGEN '$WD' $RNGNUM '$INSNP'
```

This is the modified node.sh script.

``` bash
set -o nounset
echo "host =  `hostname`"

# ===============================================
# check script parameters
function usage {
    echo -e "\tThis script runs DIYABC on a node cluster"
    echo -e "\tUsage : $0 <diyabcPath> <numberOfCores> <numToGenerate> <workingDirectory> <numIdOfTheCurrentJob> <dataFile>"
    echo -e "\texample : "
    echo -e "\t\t $0 \"/usr/local/bin/diyabc\" 1 5000 \"$PWD\" 13 \"dyabcData.mss\""
}

# Test  parameters number
if [ "$#" -lt 5 ]
    then
    echo -e "\nERROR :"
    echo -e "      Require 6  parameters, $# given\n\n"
    usage
    exit 2
fi

DIYABCPATH=$1
NBCORES=$2
NBTOGEN=$3
USERDIR=$4
MYNUMBER=$5
DATAFILE=$6
MYTMPDIR="/scratch/pawsey0149/ascheben/bitou/diyabc2/mytmp_$MYNUMBER"
#mkdir $MYTMPDIR
# Get full path for diyabc, dir and file
DIYABCPATH=`readlink -f "$DIYABCPATH"`
USERDIR=`readlink -f "$USERDIR"`
DATAFILE=`readlink -f "$DATAFILE"`

if [ ! -r "$DIYABCPATH" ]; then
    echo -e "\nERROR :" >&2
    echo -e "      diyabc not found in $DIYABCPATH or you can not read it\n" >&2
    usage
    exit 3
fi

if [ ! -r "$DATAFILE" ]; then
    echo -e "\nERROR :" >&2
    echo -e "      You do not have rights to read $DATAFILE  but you should !\n" >&2
    usage
    exit 3
fi

if [ ! -w "$USERDIR" ]; then
    echo -e "\nERROR :"  >&2
    echo -e "      You do not have the right to write in $USERDIR but you should ! (or dir does)\n" >&2
    usage
    exit 3
fi

echo -e "All parameters accepted"
# ===============================================
# start working

jobId="`hostname`-n${MYNUMBER}-pid${BASHPID}-${RANDOM}"

# initialisation values for trapActions
errStatus=""
rngFile=""
rngFileName=""
rngLockFile=""
rngFlagFile=""
rngBackupFile=""
scriptMadeTMPDIR="false"

#copy files to tmp dir
cp "$DIYABCPATH" "$MYTMPDIR/"
chmod +x "$MYTMPDIR/`basename $DIYABCPATH`"
cp header.txt "$MYTMPDIR/"
cp "$DATAFILE" "$MYTMPDIR/"

# Computations
##rngNum=`echo ${rngFileName:10} | cut -d'.' -f1`

rngNum=`expr $5 - 1`
cmd="$MYTMPDIR/`basename $DIYABCPATH` -t $NBCORES -p $MYTMPDIR/ -w $rngNum -r $NBTOGEN "
echo -e "Running computations :"
echo -e "$cmd"
"$MYTMPDIR/`basename $DIYABCPATH`" -t $NBCORES -p "$MYTMPDIR/" -w $rngNum -r $NBTOGEN &> "$MYTMPDIR/diyabc.out"
#copy results to USERDIR and clean
echo -e "computation finished with succes, bringing back results"
mv "$MYTMPDIR"/reftable.bin "$USERDIR"/reftable_$MYNUMBER.bin
mv "$MYTMPDIR"/reftable.log "$USERDIR"/reftable_$MYNUMBER.log
##mv "$MYTMPDIR/$rngFileName" "$rngFile"
echo -e "cleaning temp files"
#rm -f "$rngBackupFile" "$rngLockFile" "$rngFlagFile"
##if [ $scriptMadeTMPDIR = 'true' ]; then rm -rf "$MYTMPDIR"; fi
echo -e "End with succes"
echo -e " "
#END
exit 0
```

All the completed reftable\_${RNG}.bin files will be moved to the
working directory. They can be concatenated using the DIYABC binary as
follows. This will generate a final “reftable.bin” file.

``` bash
$diyabcPath -p "$PWD/" -q  &> concat.out
```

## ABC Random forest analysis

Based on the literature, DIYABC simulations are commonly run for 1
million or more generations per scenario. However, with a dataset of
\>1k SNPs and a large number of scenarios to test, this is not feasible.
ABC random forest analysis enables scenario selection with accuracies as
high as DIYABC but using only \~10k generations of simulations. The R
package abcrf can be used to implement this method. Three input files
generated by DIYABC are required:

1.  reftable.bin
2.  statobs.txt
3.  header.txt

These files include the simulations, the scenarios and the observed
summary statistics. To understand the below it is important to note that
DIYABC will not generate absolutely exactly equal numbers of simulations
per scenario. For example when requesting 30k simulations for 3
scenatios, DIYABC may generate 9998 for one scenario, 1002 for another
and 1000 for the third. Therefore one should also generate slightly more
scenarios and then extract an exact number during the data preparation
stage.

``` r
# Load libraries
library(abcrf)
library(data.table)
# Set working directory
setwd("D://Documents//UWA//Applied Bioinformatics//Bitou//DIYABC_2019//randomforest")

# Set total number of generations (usually 10k per scenario in the reftable)
i <- 60000

# Set names of output files
outerr <- paste(i,"_bagerr_1K_WAU_sort_best_cor_2db.txt")
outpred <- paste(i,"_prediction_1K_WAU_sort_best_cor_2db.txt")
outobject <- paste(i,"_objects_1K_WAU_sort_best_cor_2db.RData")
outplot <- paste (i, "_errplot_1K_WAU_sort_best_cor_2db.pdf")

# Load DIYABC reftable binary and corresponding header file
myref <- readRefTable(filename = "reftable_1K_WAU_sort_best_cor_2db.bin", header = "header_1K_WAU_sort_best_cor_2db.txt", N = 62000)
# Extract scenario indices
modindex <- myref$scenarios
# Extract summary statistics
sumsta <- myref$stats
# get n simulations for each modindex (=scenario)
data1 <- data.table(modindex, sumsta)
data1 <- data1[, head(.SD, 10000), by=modindex]
# Generate dataframe
data1 <- data.frame(data1)
# Generate model with 1000 trees
# The below uses five CPU cores and reads in i total simulations from the reftable
model.rf1 <- abcrf(modindex~., data = data1, ntree=1000, paral=TRUE, sampsize = i, ncores = 5)
#estimate error rate of model with ntrees and save plot
pdf(outplot) 
bagerr1 <- err.abcrf(model.rf1, data1)
dev.off()
sink(file = outerr)
print(bagerr1)
sink()
# Reading the observed dataset
obs.poi <- read.table("statobs_1K_WAU_sort_best_cor_2db.txt", header = TRUE)
# Prediction based on observed statistics
pred.obsPoi1 <- predict(object = model.rf1, obs = obs.poi,
                        training = data1, ntree=1000)
# Save prediction in file
sink(file = outpred)
print(do.call(rbind.data.frame, pred.obsPoi1))
sink()
# Save objects from analysis in R object for later use
save(model.rf1,pred.obsPoi1,bagerr1, file = outobject)
```

We can observe that with 1000 trees in the random forest the prior error
rate stablises.

``` r
setwd("D://Documents//UWA//Applied Bioinformatics//Bitou//DIYABC_2019//randomforest")
load("70000 _objects_1k_WAU_sort_best.RData")
plot(bagerr1)
```

![alt text](https://github.com/ascheben/bitou_analysis/src/images/bag_errors_plot.png "Decrease in out-of-bag error rate with increasing tree number")

We can check the outcome of the prediction. Here we can see as the
header row the scenario numbers. The votes are how many of the trees
supported each scenario. The posterior probability is only show here for
the winning scenario (1), which is shown in the allocation row.

``` r
print(do.call(rbind.data.frame, pred.obsPoi1))
```

    ##                   1        2        3       4        5       6
    ## allocation   1.0000   1.0000   1.0000  1.0000   1.0000  1.0000
    ## vote       347.0000 205.0000 186.0000 75.0000 126.0000 61.0000
    ## post.prob    0.5169   0.5169   0.5169  0.5169   0.5169  0.5169

It is sometimes recommended to ensure that the scenario outcome remains
the same with a different number of input simulations, so the analysis
can be redone with 5,000 and 7,000 simulations per scenario. If the
scenario selected remains the same, then 10k simulations are sufficient.

Finally, model checking and parameter estimation can be carried out for
the best model using DIYABC. For example, the best model is rerun, this
time with 100k-1m simulations in DIYABC. The DIYABC functions under
‘Analysis’ can then be used to see whether the posterior distribution
clusters around the observed values as it should. It is sometimes
suggested that model checking should make use of different summary
statistics from model selection, and this can be done by not giving the
whole table of summary statistics to abcrf and specifying the other
statistics to DIYABC in the model checking analysis. The number of
summary statistics significantly different from those in the posterior
distribution of the best model can also be check by selecting ‘numerical
results’ after the model check analysis. Many of these approaches are
outlined in this [Molecular Ecology
review](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1365-294X.2010.04773.x)
and in the [DIYABC
manual](http://www1.montpellier.inra.fr/CBGP/diyabc/). The parameter
estimation for the best scenario can be used to identify (with
confidence intervals) the size of populations and the time of divergence
or admixture events in the best model.
