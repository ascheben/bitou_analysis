# vim: syntax=snakemake
##execute with: snakemake --cores 8 --use-conda
import os

configfile: 'config.yaml'
length_k =  str(len(config['K']))
k_list = list(map(int, config['K']))
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
        "output/" + "RAxML_bipartitions." + config['invcf'],
        "output/RAxML_bipartitions." + config['invcf'] + "_nobrln_tree.svg",
        "output/RAxML_bipartitions." + config['invcf'] + "_brln_tree.svg",
        "output/" + config['invcf'] + "_label_pca.pdf",
        "output/" + config['invcf'] + "_pca.pdf",
        "output/structure_plot_lenK" + length_k + ".pdf",
        expand("output/structure_plot_K{K}.pdf", K = config['K']),
        "output/" + config['invcf'] + "_annot_net.svg",
        "output/" + config['invcf'] + "_net.svg"

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
        maxK = max(k_list),
        minK = min(k_list)
    shell:
        "for l in {{{params.minK}..{params.maxK}}};do structure.py -K $l --input=output/{params.prefix} --output=output/{params.prefix};done"

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
        tree = "output/" + "RAxML_bipartitions." + config['invcf'],
        map = config['popmap']
    output:
        "output/RAxML_bipartitions." + config['invcf'] + "_nobrln_tree.svg",
        "output/RAxML_bipartitions." + config['invcf'] + "_brln_tree.svg"
    shell:
        "Rscript scripts/ggtree.R {input.tree} {input.map}"

rule pca:
    input:
        vcf = "output/" + config['invcf'] + "_mim" + str(config['im_cut1']) + "_biallelic_minDP" + str(config['mindp']) + "_mm" + str(config['mm']) + "_maf" + str(config['maf']) + "_thin" + str(config['thin']) + "_mim" + str(config['im_cut2']) + ".vcf",
        map = config['popmap']
    output:
        "output/" + config['invcf'] + "_label_pca.pdf",
        "output/" + config['invcf'] + "_pca.pdf"
    params:
         os.getcwd() + "/output/" + config['invcf']
    shell:
        "Rscript scripts/pca.R {input.vcf} {input.map} {params}"

rule plot_structure:
    output:
        allk = "output/structure_plot_lenK" + str(len(config['K'])) + ".pdf",
        singlek = expand("output/structure_plot_K{K}.pdf",
                        K = config['K'])
    params:
        lenK = str(len(config['K'])),
        path = os.getcwd() + "/output",
        names = config['names'],
        groups = config['groups']
    shell:
        "Rscript scripts/plotStructure.R {params.path} {params.names} {params.groups} {params.lenK}"

rule plot_network:
    input:
        nnet =  config['nnet']
    output:
        annot = "output/" + config['invcf'] + "_annot_net.svg",
        plain= "output/" + config['invcf'] + "_net.svg"

    params:
        path = os.getcwd() + "/output",
        base = config['invcf']
    shell:
        "Rscript scripts/network.R {input.nnet} {params.base} {params.path}"

# Snakemake notification
onerror:
  print("Error: Snakemake aborted!")
#  shell("mail -s 'Snakemake Job Error: See log inside!' {config[email]} < {log}")


onsuccess:
  print("Success: Snakemake completed!")
#  shell("mail -s 'Snakemake Job Completed: Have a Beer!' {config[email]} < {log}")

