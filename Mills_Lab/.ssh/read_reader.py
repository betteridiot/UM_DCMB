#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import fnmatch


def catcher(path, keyword):
    returns = []
    for d, s, f in os.walk(path):
        lister = fnmatch.filter(f, keyword)
        returns.extend(os.path.join(d, f) for f in lister)
    return returns

def main():
    snp_chrom, snp_pos = sys.argv[1].split(":")
    region = "{}-{}".format(int(snp_pos)-1000, int(snp_pos)+1000)

    reference = '/home/mdsherm/Project/Reference/hg19/Sequence/Chromosomes'

    rna_path = '/home/mdsherm/Rotation/RNA_fq/tophat_hg19'
    rna_bams = catcher(rna_path, "accepted_hits.bam")
    rna_cross = {rna_bam:rna_bam.split('_')[2].split('/')[1] for rna_bam in rna_bams}

    ribo_path = '/home/mdsherm/Rotation/ribosomal/bwa_alignment_hg19'
    ribo_bams = catcher(ribo_path, "*_sort.bam")
    ribo_cross = {}
    for row in open('/home/mdsherm/Project/ribosamplescross', 'rb'):
        for ribo_bam in ribo_bams:
            if ribo_bam.split('/')[6].split('_')[0] in row.split()[0]:
                ribo_cross.update({ribo_bam:row.split()[1]})

    rna = [(rna_cross.get(f) , (os.popen("samtools view -q 10 {} {}:{} | wc -l".format(f, snp_chrom, region)))) for f in rna_bams]
    ribo = [(ribo_cross.get(f) , (os.popen("samtools view -q 10 {} {}:{} | wc -l".format(f, snp_chrom, region)))) for f in ribo_cross.]


if __name__ == "__main__":
    main()
