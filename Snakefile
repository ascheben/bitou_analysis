# vim: syntax=snakemake
##execute with: snakemake --cores 8 --use-conda
import os

configfile: 'config.yaml'

K = ["1","2"]
# main rule, specify exact output file
#Files which should be generated in the course of the analysis
rule all:
    input:
        "output/" + config["invcf"] + ".imiss",
        "output/" + config['invcf'] + "_mim" + str(config['im_cut1']) + "_biallelic_minDP" + str(config['mindp']) + "_mm" + str(config['mm']) + "_maf" + str(config['maf']) +".vcf",
        "output/" + config['invcf']  + "_mim" + str(config['im_cut1']) + "_biallelic_minDP" + str(config['mindp']) + "_mm" + str(config['mm']) + "_maf" + str(config['maf']) + ".imiss",
        "output/" + config['invcf'] + "_mim" + str(config['im_cut1']) + "_biallelic_minDP" + str(config['mindp']) + "_mm" + str(config['mm']) + "_maf" + str(config['maf']) + "_thin" + str(config['thin']) + "_mim" + str(config['im_cut2']) + ".vcf",
        "output/" + config['invcf'] + ".bed",
        "output/" + config['invcf'] + ".bim",
        "output/" + config['invcf'] + ".fam",
        expand("output/{sample}.{K}.{ext}",
            sample = config['invcf'],
            K = config['K'],
            ext = ['meanQ','meanP']),
        "output/" + config['invcf'] + "_hetfilt.vcf",
        "output/" + config['invcf'] + "_hetfilt.phy",
        "output/" + "RAxML_bipartitions." + config['invcf']

# sub-rules
rule filter_vcf_1:
    input:
        config["invcf"] + ".vcf"
    output:
        "output/" + config["invcf"] + ".imiss"
    shell:
        "vcftools --vcf {input} --missing-indv --stdout | tail -n +2 | awk '$5>{config[im_cut1]}' | cut -f1 > {output}"

rule filter_vcf_2:
    input:
        vcf = config["invcf"] + ".vcf",
        imiss = "output/" + config["invcf"] + ".imiss"
    output:
        "output/" + config['invcf'] + "_mim" + str(config['im_cut1']) + "_biallelic_minDP" + str(config['mindp']) + "_mm" + str(config['mm']) + "_maf" + str(config['maf']) +".vcf"
    shell:
        "vcftools --vcf {input.vcf} --remove {input.imiss} --maf {config[maf]} --max-missing {config[mm]} --remove-indels --max-alleles 2 --min-alleles 2 --minDP {config[mindp]} --recode --stdout > {output}"

rule filter_vcf_3:
    input:
        "output/" + config['invcf'] + "_mim" + str(config['im_cut1']) + "_biallelic_minDP" + str(config['mindp']) + "_mm" + str(config['mm']) + "_maf" + str(config['maf']) +".vcf"
    output:
        "output/" + config['invcf']  + "_mim" + str(config['im_cut1']) + "_biallelic_minDP" + str(config['mindp']) + "_mm" + str(config['mm']) + "_maf" + str(config['maf']) + ".imiss"
    shell:
        "vcftools --vcf {input} --missing-indv --stdout | tail -n +2 | awk '$5>{config[im_cut2]}' | cut -f1 > {output}"

rule filter_vcf_4:
    input:
        vcf = "output/" + config['invcf'] + "_mim" + str(config['im_cut1']) + "_biallelic_minDP" + str(config['mindp']) + "_mm" + str(config['mm']) + "_maf" + str(config['maf']) +".vcf",
        imiss = "output/" + config['invcf'] + "_mim" + str(config['im_cut1']) + "_biallelic_minDP" + str(config['mindp']) + "_mm" + str(config['mm']) + "_maf" + str(config['maf']) +".imiss"
    output:
        "output/" + config['invcf'] + "_mim" + str(config['im_cut1']) + "_biallelic_minDP" + str(config['mindp']) + "_mm" + str(config['mm']) + "_maf" + str(config['maf']) + "_thin" + str(config['thin']) + "_mim" + str(config['im_cut2']) + ".vcf"
    shell:
        "vcftools --vcf {input.vcf} --thin 500 --recode --stdout | vcftools --vcf - --remove {input.imiss} --recode --stdout  > {output}"


#vim command required to non-greedily replace chr names
rule vcf2plink:
    input: "output/" + config['invcf'] + "_mim" + str(config['im_cut1']) + "_biallelic_minDP" + str(config['mindp']) + "_mm" + str(config['mm']) + "_maf" + str(config['maf']) + "_thin" + str(config['thin']) + "_mim" + str(config['im_cut2']) + ".vcf"
    output:
        protected("output/" + config['invcf'] + ".bed"),
        protected("output/" + config['invcf'] + ".bim"),
        protected("output/" + config['invcf'] + ".fam")
    params:
        prefix = config['invcf']

    shell:
        """
        vim -c 'g/^[^#]/s/.\{{-}}\t/chrUn\t/' -c 'wq' {input}
        plink --vcf {input} --double-id --allow-extra-chr --recode --out output/{params.prefix} --make-bed
        """
#trying to execute all values of K in parallel
#snakemake keeps putting all the values in a single command so it fails
rule fastStructure:
    output:
        "output/{sample}.{K}.meanQ",
        "output/{sample}.{K}.meanP"
    params:
        prefix = config['invcf'],
        maxK = config['maxK']
    shell:
        "for l in {{1..{params.maxK}}};do structure.py -K $l --input=output/{params.prefix} --output=output/{params.prefix};done"

# Infer ML phylogeny from SNPs
rule remove_het_snps:
    input: "output/" + config['invcf'] + "_mim" + str(config['im_cut1']) + "_biallelic_minDP" + str(config['mindp']) +     "_mm" + str(config['mm']) + "_maf" + str(config['maf']) + "_thin" + str(config['thin']) + "_mim" + str(config['im_cut2']) + ".vcf"
    output:
       "output/" + config['invcf'] + "_hetfilt.vcf"
    params:
        minalt = config['minalt'],
        maxhet = config['maxhetprop']
    shell:
        "./scripts/filter_hets.py {input} {params.maxhet} {params.minalt} > {output}"

rule vcf2phy:
    input:
        "output/" + config['invcf'] + "_hetfilt.vcf"
    output:
        "output/" + config['invcf'] + "_hetfilt.phy"
    shell:
        "./scripts/vcf2phylip/vcf2phylip.py -i {input}"

rule raxml:
    input:
        "output/" + config['invcf'] + "_hetfilt.phy"
    output:
        "output/" + "RAxML_bipartitions." + config['invcf']
    params:
        bs = config['bs'],
        outgroup = config['outgroup'],
        prefix = config['invcf'],
        outdir = os.getcwd() + "/output"
    threads: 4
    shell:
        "raxmlHPC-PTHREADS-SSE3 -f a -V -T {threads} -m ASC_GTRCAT --asc-corr lewis -p 12345 -x 12345 -# {params.bs} -s {input} -n {params.prefix} -o {params.outgroup} -w {params.outdir}"

rule ggtree:
    input:
        "output/" + "RAxML_bipartitions." + config['invcf'],
        config['popmap']
    output:
        "output/" + config['invcf'] + "_tree.svg",
        "output/" + config['invcf'] + "_tree.pdf"


# Snakemake notification
onerror:
  print("Error: Snakemake aborted!")
#  shell("mail -s 'Snakemake Job Error: See log inside!' {config[email]} < {log}")


onsuccess:
  print("Success: Snakemake completed!")
#  shell("mail -s 'Snakemake Job Completed: Have a Beer!' {config[email]} < {log}")

Rscript ../code/plot_network.R /home/arminps/ws/ddrad/bitou/snake/bitou_m3_ibs.dst.nex bitou_core_raw /home/arminps/ws/ddrad/bitou/snake/output
Rscript ../code/plot_structure.R /home/arminps/ws/ddrad/bitou/snake/output/ /home/arminps/ws/ddrad/bitou/snake/popmap_core.names.txt /home/arminps/ws/ddrad/bitou/snake/popmap_core.groups.txt 2
Rscript /home/arminps/ws/ddrad/bitou/code/ggtree_bitou.R /home/arminps/ws/ddrad/bitou/snake/output/RAxML_bipartitions.bitou_core_raw /home/arminps/ws/ddrad/bitou/snake/popmap_core.txt
