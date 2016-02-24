#!/usr/bin/env python

from __future__ import print_function
import os
import gzip
import glob
import time
# import argparse
import math
import csv
import string
import random
from multiprocessing.dummy import Pool as ThreadPool

"""
Created on Thu Jan 28 10:58:33 2016
Last Modified on Sat Feb 15 12:22 2016

Description: ORFSNuPer uses identified SNPs from the 1000 Genomes Project that meet certain MAF criteria related to a
    given population. From that SNP, the hg19 reference genome is indexed with regard to the identified SNP and looks
    to see if a potential ORF was created. IFF an ORF is identified, ORFSNuPer looks upstream and downstream for stop
    codons. IFF a potential novel reading frame is found, using reads from RNA-seq and ribosome profiling data,
    ORFSNuPer will check to see if transcription and/or translation occurs.
@author: Marcus D. Sherman
@email: mdsherm@umich.edu
"""
startTime, startasc = time.time(), time.asctime()

# # ARGPARSE START
# parser = argparse.ArgumentParser(description='Finds novel ORFs dues to SNPs')
# # define where the reference sequence is
# parser.add_argument('-r', action='store', dest='ref', help='Directory of reference chromosomes',
#                     default='/home/mdsherm/Project/Reference/hg19/Sequence/Chromosomes')
# # define what VCF file you will be working from
# parser.add_argument('-v', action='store', dest='vcf', help='Path/to/<vcf.gz>',
#                     default='/home/mdsherm/Project/YRI_vcfsubsets/filteredGenotypeVCF/unannotatedchr22.vcf.gz')
# # how many nucleotides do you want to look upstream and downstream of potential ORFs
# parser.add_argument('-t', action='store', dest='threshold', type=int, help='Up/downstream threshold', default=3000)
# # Define your output filename and directory
# parser.add_argument('-o', action='store', dest='output', help='Set output filename',
#                     default='/home/mdsherm/Project/SNuPer_results/pythonTest')
# # Where are the ribosome profiling BAM files
# parser.add_argument('--ribosome', action='store', dest='ribo', help='Directory of Ribosomal BAM files',
#                     default='/home/mdsherm/Rotation/ribosomal/bwa_alignment')
# # Where are the RNA-seq BAM files
# parser.add_argument('--rna', action='store', dest='rna', help='Directory of RNA BAM files',
#                     default='/home/mdsherm/Rotation/RNA_fq/tophat_hg19')
# parser.add_argument('--alternative', action='store_true', help='Use alternative start codon GTG', default='False')
# args = parser.parse_args()
# vcf = args.vcf
# riboDir = args.ribo
# RNADir = args.rna
# outfile = args.output
# reference = args.ref
# threshold = args.threshold
reference = '/home/mdsherm/Project/Reference/hg19/Sequence/Chromosomes'
vcf = '/home/mdsherm/Project/YRI_vcfsubsets/filteredGenotypeVCF/unannotatedchr22.vcf.gz'
RNADir = '/home/mdsherm/Rotation/RNA_fq/tophat_hg19'
riboDir = '/home/mdsherm/Rotation/ribosomal/bwa_alignment'
outDir = '/home/mdsherm/Project/SNuPer_results/pythonTest/'
threshold = 3000
orfcount = 0  # use when debugging
# ARGPARSE END

# Adding in ribosome sample name to SRA ID file
ribosamples = []
with open('/home/mdsherm/Project/ribosamplescross') as ribo:
    reader = csv.reader(ribo, delimiter='\t')
    ribosamples.extend([row for row in reader])

# # CODONS START
negStops = ['TTA', 'CTA', 'TCA']
plusStops = ['ATT', 'ATC', 'ACT']
# if args.alternative is True:
#     plusStart, negStart = ('TAC', 'CAC'), ('CAT', 'CAC')
# else:
#     plusStart = ("TAC",)
#     negStart = ("CAT",)

plusStart = ("TAC",)
negStart = ("CAT",)
# CODONS END

# makes lists of all RNA-seq and ribosome profiling BAM files
RNAbams, Ribobams = [], []
RNAbams_total, Ribobams_total = [], []
# dir, _, _ implies dir, dirs, files...however, I only want the top object
os.chdir(RNADir)
for dir, _, _ in os.walk(os.getcwd()):
    RNAbams.extend(glob.glob(os.path.join(dir, "*hits.bam")))
os.chdir(riboDir)
for dir, _, _ in os.walk(os.getcwd()):
    Ribobams.extend(glob.glob(os.path.join(dir, "*sort.bam")))
with open("/home/mdsherm/Project/RNAsamples", 'r') as rnasamples:
    samplenames = rnasamples.read().splitlines()
RNAsamp_crossref = [(name, [entry for entry in RNAbams if name in entry]) for name in samplenames]
Ribosamp_crossref = [(name[1], [entry for entry in Ribobams if name[0] in entry]) for name in ribosamples]


# Make a dict of Samples and FPKM
def read_indexer(fileDIR, LIST):
    for item in fileDIR:
        fullcounter = os.popen("samtools idxstats " + item + " |awk '{sum+=$3} END {print sum}'")
        fullcount = float(fullcounter.readline().rstrip())
        LIST.extend([[item, fullcount]])
    return {key: value for (key, value) in LIST}


# Matches sample names to bam files and pulls out total reads of
RNAbams_dict = read_indexer(RNAbams, RNAbams_total)
RNAtotalreads = dict(zip(samplenames, [sum([RNAbams_dict.get(j) for j in i[1]])for i in RNAsamp_crossref]))
Ribobams_dict = read_indexer(Ribobams, Ribobams_total)
Ribototalreads = dict(zip(samplenames, Ribobams_dict.values()))


# Used to iterate potential ORF class instantiation
def portORF(CHROM, START, END, ID, STRAND, GENO):
    portedORF = potORF(CHROM, START, END, ID, STRAND, GENO)
    return portedORF


# iteratively pulls FPKM over a given region across all BAMs
def readCheck(RNAorRIBO, CHROM, START, STOP, LENGTH):
    bamlist, samplereads = [], []
    typeCheck = RNAorRIBO
    if typeCheck:  # True = RNA-seq, False = Ribosome profiling
        bamlist = [i[1] for i in RNAsamp_crossref]
        bamlist_dict = RNAbams_dict
    elif typeCheck is False:
        bamlist = [i[1] for i in Ribosamp_crossref]
        bamlist_dict = Ribobams_dict
    for filelist in bamlist:
        counter, fullcount = [], []
        for filename in filelist:
            readcount = os.popen('samtools view -q 10 ' + filename + ' chr%d:%d-%d | wc -l'
                                 % (int(CHROM), int(START), int(STOP)))
            counter.append(float(readcount.readline().rstrip()))
            fullcount.append(bamlist_dict[filename])
        try:
            samplereads.append((sum(counter)/(sum(fullcount)*LENGTH))*math.pow(10, 9))
        except ZeroDivisionError:
            samplereads.append(float(0))
    if sum(samplereads) == (float(0) or None):
        return "NA"
    else:
        return samplereads


# Pull out the reference sequence at a given position
def SNP_search(CHROM, START, STOP):
    seq = os.popen('samtools faidx %s/chr%s.fa chr%s:%d-%d' % (reference, CHROM, CHROM, START, STOP))
    seq.readline()
    seq = seq.read().rstrip()
    return seq


def strand_checker(SEQ, ROW, LIST):
    global orfcount
    posCheckplus = sum([SEQ.find(codon) for codon in plusStart if codon in SEQ])+1
    posCheckneg = sum([SEQ.find(codon) for codon in negStart if codon in SEQ])+1
    if posCheckplus > 0:
        seqPos = int(ROW[1]) - posCheckplus
        LIST.append(portORF(ROW[0], seqPos, seqPos + 2, ROW[2], True, ROW[9:]))
        # orfcount += 1

    # Check to see if (-) strand ORF is found
    elif posCheckneg > 0:
        seqPos = int(ROW[1]) + posCheckneg
        LIST.append(portORF(ROW[0], seqPos, seqPos - 2, ROW[2], False, ROW[9:]))
        # orfcount += 1
    else:
        pass


class potORF(object):
    """Create an potential ORF object with the following attributes:

    chrom: identifies which chromosome the potential ORF is on
    start: position of the first nt in start codon
    end: postion of the last nt in the start codon
    strand: (+) or (-) strand
    up: upstream UTR of potential ORF
    upCheck: is there an upstream stop
    down: downstream sequence of potential ORF
    downCheck: is there a downstream stop
    up/downPos: how many codons away is a given stop codon (relative to strand)
    RNAcount/ribocount = how many reads were present in the region of the ORF
    Homo/heterocount = Sum of the samples that are homozygous or heterozygous for SNP
    Length = length of the ORF
    """

    # instantiate potential ORF object and get upstream and downstream sequences
    def __init__(self, CHROM, START, END, ID, STRAND, GENO):
        self.chrom = CHROM
        self.start = START
        self.end = END
        self.SNPid = ID
        self.strand = STRAND
        self.genotypes = GENO
        self.length = 0
        self.up = os.popen('samtools faidx %s/chr%s.fa chr%s:%d-%d' % (
            reference, self.chrom, self.chrom, int(self.start) - threshold, int(self.start) - 1))
        self.up.readline()
        self.up = ((self.up.read()).rstrip()).upper().replace('\n', '')
        self.down = os.popen('samtools faidx %s/chr%s.fa chr%s:%d-%d' % (
            reference, self.chrom, self.chrom, int(self.end) + 1, int(self.end) + threshold))
        self.down.readline()
        self.down = ((self.down.read()).rstrip()).upper().replace('\n', '')
        self.upcheck, self.downcheck = False, False
        self.upPos, self.downPos = [], []
        self.RNAcount, self.ribocount = None, None

    # check to see if a stop codon is within upstream sequence (downstream if (-) strand)
    def lookUp(self):
        codonCount = 0
        if self.strand:
            stops = plusStops
        else:
            stops = negStops
        for w in range(0, threshold, 3):
            codon = self.up[w:w + 3]
            codonCount += 1
            if any(string in codon for string in stops):
                self.upcheck = True
                self.upPos.append(codonCount)
        return self

    # check to see if a stop codon is within downsteam sequence (upstream if (-) strand)
    def lookDown(self):
        codonCount = 0
        if self.strand:
            stops = plusStops
        else:
            stops = negStops
        for z in range(0, threshold, 3):
            codon = self.down[z:z + 3]
            codonCount += 1
            if any(string in codon for string in stops):
                self.downcheck = True
                self.downPos.append(codonCount)
        return self

    # Looks at all specified BAM files and determine the number of reads found in a given region and corrects for length
    def WordCount(self):
        if self.upcheck and self.downcheck:
            if self.strand:  # is it a (+) strand?
                begin, end = int(self.start), int(self.start) + (int(self.downPos[0]) * 3) + 2
                self.length = end-begin
                self.RNAcount = readCheck(True, int(self.chrom), begin, end, self.length)
                if self.RNAcount == ("NA" or float(0)):
                    pass
                else:
                    self.ribocount = readCheck(False, int(self.chrom), begin, end, self.length)

            else:  # is it a (-) strand?
                begin, end = int(self.start - int(self.upPos[-1]) * 3), int(self.start)
                self.length = end-begin
                self.RNAcount = readCheck(True, int(self.chrom), begin, end, self.length)
                if self.RNAcount == ("NA" or float(0)):
                    pass
                else:
                    self.ribocount = readCheck(False, int(self.chrom), begin, end, self.length)


# Removes bams from the list to be checked if they are not in the samples of the vcf
def popper(LIST):
    index, counter = 0, []
    for i in range(len(LIST)):
        if LIST[i][0] in samples:
            index += 1
            pass
        else:
            counter.append(index)
            index += 1
    for i in counter[::-1]:
        LIST.pop(i)


# function identifies SNPs, extracts sequence from reference +/-2 nt and looks for start codon within sequence
def ORFSNuper():
    global potORFs
    # global orfcount
    potORFs = []
    global samples
    with gzip.open(vcf, 'rt')as VCF:
        # while orfcount < 50:   # use when debugging
        for line in VCF:
            if "##" in line:
                continue
            elif "#CHROM" in line:
                header = line.split()
                samples = header[9:]
                popper(RNAsamp_crossref)
                popper(Ribosamp_crossref)
            else:
                columns = line.split()

                # Check to see if it is a SNP
                if len(columns[3]) and len(columns[4]) == 1:

                    # look for the reference sequence around SNP
                    seq = SNP_search(columns[0], int(columns[1]) - 2, int(columns[1]) + 2)
                    seq_step = (seq[:1] + columns[4] + seq[3:]).upper()

                    # iff sequence creates start codon does it make a class instance
                    strand_checker(seq_step, columns, potORFs)
                else:
                    continue
                    # For debugging
                    # if orfcount >= 50:
                    #     print("orfcount met!")
                    #     break


ORFSNuper()


# Make a unique filename for each SNP written
def id_generator(size=7, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))+".snp"


def file_writer(POTORF):
    try:
        filename = outDir + id_generator()
        if POTORF.strand:
            info = [str(POTORF.chrom), "+", str(POTORF.start), str((POTORF.start + int(POTORF.downPos[0]) * 3)+2)]
        else:
            info = [str(POTORF.chrom), "-", str(POTORF.start), str(POTORF.start - int(POTORF.upPos[-1]) * 3)]
        with open(filename, 'w') as writer:
            header1 = ["#ID", "CHROM", "STRAND", "START", "Nearest_STOP"]
            print('\t'.join(header1), file=writer)
            print('\t'.join(info), file=writer)
            header2 = ["Sample", "Genotype", "RNA_FPKM", "Ribo_FPKM"]
            print('\t'.join(header2), file=writer)
            SNuPed = [(a, b, c, d) for i, (a, b, c, d) in enumerate(zip(samples, POTORF.genotypes,
                                                                        POTORF.RNAcount, POTORF.ribocount))]
            writer.write('\n'.join('%s\t%s\t%s\t%s' % x for x in SNuPed))
    except TypeError:
        pass


# Look upstream and downstream for stop codons and read count of ribosome/RNA bam files by class instance multithreading
pool = ThreadPool()
pool.map(lambda obj: obj.lookUp().lookDown().WordCount(), potORFs)
pool.close()
pool.join()

pool = ThreadPool()
pool.map(lambda obj: file_writer(obj), potORFs)
pool.close()
pool.join()

# # Write each SNP to file
# for SNP in potORFs:
#     file_writer(SNP)

# find out how long the process took
endTime, endasc = time.time(), time.asctime()
m, s = divmod(endTime - startTime, 60)
h, m = divmod(m, 60)
d, h = divmod(h, 24)

# Write a small report for start time, end time, and elapsed time
with open(outDir + "out.log", 'w') as f:
    print("Program started:", file=f)
    print(startasc, file=f)
    print('', file=f)
    print("Program completed:", file=f)
    print(endasc)
    print('', file=f)
    print("Elapsed time:", file=f)
    print("%s days, %s hours, %s seconds" % (str(h), str(m), str(round(s, 2))), file=f)
    print('', file=f)