#!/usr/bin/env python3

import sys
import csv

vcf = sys.argv[1] #input vcf file
cutoff = float(sys.argv[2]) #max proportion heterozygous alleles
minalt = int(sys.argv[3]) #min number of alternative alleles

with open(vcf) as tsv:
    for line in csv.reader(tsv, dialect="excel-tab"):
        if not line[0].startswith("#"):
            num_info_fields = 9
            num_samples = len(line) - num_info_fields
            het_count = 0
            alt_count = 0
            hom_count = 0
            for allele in line[num_info_fields:]:
                g = allele.split(":")[0]
                if g == "0/1":
                    het_count = het_count + 1
                elif g == "1/1":
                    alt_count = alt_count + 1
                elif g == "0/0":
                    hom_count = hom_count + 1
            if (het_count / num_samples) <= cutoff and alt_count >= minalt :
                print(*line, sep = '\t')
        else:
            print(*line, sep = '\t')
