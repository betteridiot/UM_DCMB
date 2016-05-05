#!/usr/bin/env python
from __future__ import print_function
import os
import csv
import sys
import fnmatch
import numpy as np
# correction

def catcher(path, keyword):
    returns = []
    for d, s, f in os.walk(path):
        lister = fnmatch.filter(f, keyword)
        returns.extend(os.path.join(d, f) for f in lister)
    return returns

def main():
    if "-h" in sys.argv[1]:
        print('Usage: ./read_reader.py chr#_SNPpos:ID')
    else:
        snp_pos, snp_id = sys.argv[1].split(':')
        snp_chrom, snp_pos = snp_pos.split('_')
        region = "{}-{}".format(int(snp_pos)-1000, int(snp_pos)+1000)

        reference = '/home/mdsherm/Project/Reference/hg19/Sequence/Chromosomes'
        project_root = '/home/mdsherm/Project/SNuPer_results'
        snp_file = os.popen('locate -d {}/chr{}/mlocate.db {}.snp'.format(
            project_root, snp_chrom, sys.argv[1])).read().rstrip()
        snp_geno = {row[0]:row[1] for row in csv.reader(open(snp_file, "rb"), delimiter='\t')}

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

        rna = [(rna_cross.get(f) , int(os.popen(
            "samtools view -q 10 {} {}:{} | wc -l".format(f, snp_chrom, region)).rstrip())
                ) for f in rna_bams]
        rna_noreps = [(sample[0],) for sample in rna]
        ribo = [(ribo_cross.get(f) , int(os.popen(
            "samtools view -q 10 {} {}:{} | wc -l".format(f, snp_chrom, region)))
                 ) for f in ribo_cross]

        mix_cross = [(rnsample[0], snp_geno.get(rnsample[0]), rnsample[1], rbsample[1]) for rnsample in rna
                     for rbsample in ribo if rbsample[0] in rnsample[0]]
        mix_ref = [row for row in mix_cross if "0|0" in row[1]]
        mix_ref = [(row[2], row[3]) for row in mix_ref]
        mix_het = [row for row in mix_cross if "1|0" in row[1] or "0|1" in row[1]]
        mix_het = [(row[2], row[3]) for row in mix_het]
        mix_alt = [row for row in mix_cross if "1|1" in row[1]]
        mix_alt = [(row[2], row[3]) for row in mix_alt]

        print('RNA mean:std / Ribo mean:std')
        print('REF:: {}:{} / {}:{}'.format(
            np.mean(mix_ref, axis=0)[0], np.std(mix_ref, axis=0)[0],
            np.mean(mix_ref, axis=0)[1], np.std(mix_ref, axis=0)[1]))
        print('HET:: {}:{} / {}:{}'.format(
            np.mean(mix_het, axis=0)[0], np.std(mix_het, axis=0)[0],
            np.mean(mix_het, axis=0)[1], np.std(mix_het, axis=0)[1]))
        print('ALT:: {}:{} / {}:{}'.format(
            np.mean(mix_alt, axis=0)[0], np.std(mix_alt, axis=0)[0],
            np.mean(mix_alt, axis=0)[1], np.std(mix_alt, axis=0)[1]))

# TODO command for building locate DB:
# updatedb -U /home/mdsherm/Project/SNuPer_results/chr10 -o /home/mdsherm/Project/SNuPer_results/chr10/mlocate.db --require-visibility 0

# TODO command to use locate:
# locate -d /home/mdsherm/Project/SNuPer_results/chr10/mlocate.db 10_15429889:XY.snp

if __name__ == "__main__":
    main()


a = [('asdf', 7), ('asdf', 9), ('jkl;', 8)]
reduct = []
for b in  range(len(a)): print(a[b])


