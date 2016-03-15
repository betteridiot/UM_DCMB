#! /usr/bin/env python

from __future__ import print_function
import os
import gzip
# import glob
import time
import argparse
import csv
import dill as pickle
import numpy as np
from random import randint
from pathos.multiprocessing import ProcessingPool as Pool
# from multiprocessing import Pool, cpu_count
from multiprocessing.dummy import Pool as ThreadPool

"""
Created on Thu Jan 28 10:58:33 2016
Last Modified on Thur March 10 16:30 2016

Description: ORFSNuPer uses identified SNPs from the 1000 Genomes Project that meet certain MAF criteria related to a
    given population. From that SNP, the hg19 reference genome is indexed with regard to the identified SNP and looks
    to see if a potential ORF was created. IFF an ORF is identified, ORFSNuPer looks upstream and downstream for stop
    codons. IFF a potential novel reading frame is found, using reads from RNA-seq and ribosome profiling data,
    ORFSNuPer will check to see if transcription and/or translation occurs.
@author: Marcus D. Sherman
@email: mdsherm@umich.edu
"""
startTime, startasc = time.time(), time.asctime()

# ARGPARSE START
parser = argparse.ArgumentParser(description='Finds novel ORFs dues to SNPs')
# define where the reference sequence is
parser.add_argument('-r', action='store', dest='ref', help='Directory of reference chromosomes',
                    default='/home/mdsherm/Project/Reference/hg19/Sequence/Chromosomes')
# define what VCF file you will be working from
parser.add_argument('-v', action='store', dest='vcf', help='Path/to/<vcf.gz>',
                    default='/home/mdsherm/Project/YRI_vcfsubsets/filteredGenotypeVCF/unannotatedchr22.vcf.gz')
# how many nucleotides do you want to look upstream and downstream of potential ORFs
parser.add_argument('-t', action='store', dest='threshold', type=int, help='Up/downstream threshold', default=3000)
# Define your output filename and directory
parser.add_argument('-o', action='store', dest='output', help='Set output path',
                    default='/home/mdsherm/Project/SNuPer_results/pythonTest/')
# Where are the ribosome profiling BAM files
parser.add_argument('--ribosome', action='store', dest='ribo', help='Directory of Ribosomal BAM files',
                    default='/home/mdsherm/Rotation/ribosomal/bwa_alignment')
# Where are the RNA-seq BAM files
parser.add_argument('--rna', action='store', dest='rna', help='Directory of RNA BAM files',
                    default='/home/mdsherm/Rotation/RNA_fq/tophat_hg19')
parser.add_argument('--alternative', action='store_true', help='Use alternative start codon GTG', default='False')
args = parser.parse_args()
vcf = args.vcf
riboDir = args.ribo
RNADir = args.rna
outDir = args.output
reference = args.ref
threshold = args.threshold
# reference = '/home/mdsherm/Project/Reference/hg19/Sequence/Chromosomes'
# vcf = '/home/mdsherm/Project/YRI_vcfsubsets/filteredGenotypeVCF/testing_smll.vcf.gz'
# RNADir = '/home/mdsherm/Rotation/RNA_fq/tophat_hg19'
# riboDir = '/home/mdsherm/Rotation/ribosomal/bwa_alignment'
# outDir = '/home/mdsherm/Project/SNuPer_results/pythonTest/'
# threshold = 3000
# orfcount = 0  # use when debugging
# ARGPARSE END

# CODONS START
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


list_files = ["RNAsamp_crossref", "Ribosamp_crossref", "RNAbams_dict", "Ribobams_dict"]
bamDir = "/home/mdsherm/Project/UM_DCMB/Mills_Lab/pkl/"
for name in list_files:
    with open('%s%s.pkl' % (bamDir, name), 'rb') as f:
        n = name
        globals()[n] = pickle.load(f)


# Used to iterate potential ORF class instantiation
def portORF(CHROM, START, END, ID, STRAND, SNP_POS, GENO):
    portedORF = potORF(CHROM, START, END, ID, STRAND, SNP_POS, GENO)
    return portedORF


# iteratively pulls FPKM over a given region across all BAMs
def readCheck(RNAorRIBO, CHROM, START, STOP, LENGTH):
    """Iteratively pulls FPKM over a given region across all BAMs.

    Args:
        RNAorRIBO (bool): True if RNA, False for ribosome
        CHROM (int): What chromosome the SNP is on
        START (int): The position of the first nt of start codon
        STOP (int): The position of the last nt of stop codon
        LENGTH (int): How long the ORF is
    """

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
            samplereads.append((sum(counter)/(sum(fullcount)*LENGTH))*10**9)
        except ZeroDivisionError:
            samplereads.append(0)
    if RNAorRIBO and max(samplereads) < 1:
        return None
    else:
        # print("Good read count ORF")
        return samplereads


# Pull out the reference sequence at a given position
def SNP_search(CHROM, START, STOP):
    """Pulls out the sequence +/- 2 nt around SNP from reference genome.

    Args:
        CHROM (int): The chromosome the SNP is on
        START (int): 2 nt upstream of SNP
        STOP (int): 2 nt downstream of SNP
    """

    seq = os.popen('samtools faidx %s/chr%s.fa chr%s:%d-%d' % (reference, CHROM, CHROM, START, STOP))
    seq.readline()
    seq = seq.read().rstrip()
    return seq


def strand_checker(SEQ, ROW, LIST):
    """Looks to see if the SNP ORF has a start codon for both the plus and minus strands.
    If a start codon is found, it pulls upstream and downstream sequences of the same size
    as set in the threshold parameter.

    Args:
        SEQ (str): The sequence +/- 2 nt of SNP
        ROW (list): The row for the SNP from the VCF
        LIST (list): the list that will hold the generated potORFs
    """

    # global orfcount
    posCheckplus = sum([SEQ.find(codon) for codon in plusStart if codon in SEQ])+1
    posCheckneg = sum([SEQ.find(codon) for codon in negStart if codon in SEQ])+1
    if not 1 <= ROW.count("1|1") < len(ROW[9:]):
        pass
    else:
        if posCheckplus > 0:
            seqPos = int(ROW[1]) - posCheckplus
            LIST.append(portORF(ROW[0], seqPos, seqPos + 2, ROW[2], True, ROW[1], ROW[9:]))
            # orfcount += 1

        # Check to see if (-) strand ORF is found
        elif posCheckneg > 0:
            seqPos = int(ROW[1]) + posCheckneg
            LIST.append(portORF(ROW[0], seqPos, seqPos - 2, ROW[2], False, ROW[1], ROW[9:]))
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
    def __init__(self, CHROM, START, END, ID, STRAND, SNP_POS, GENO):
        self.chrom = CHROM
        self.start = START
        self.end = END
        self.SNPid = ID
        self.strand = STRAND
        self.SNPpos = SNP_POS
        self.genotypes = GENO
        self.length = None
        up = os.popen('samtools faidx %s/chr%s.fa chr%s:%d-%d' % (
            reference, self.chrom, self.chrom, int(self.start) - threshold, int(self.start) - 1))
        up.readline()
        self.up = ((up.read()).rstrip()).upper().replace('\n', '')
        down = os.popen('samtools faidx %s/chr%s.fa chr%s:%d-%d' % (
            reference, self.chrom, self.chrom, int(self.end) + 1, int(self.end) + threshold))
        down.readline()
        self.down = ((down.read()).rstrip()).upper().replace('\n', '')
        self.upcheck, self.downcheck = False, False
        self.upPos, self.downPos = [], []
        self.RNAcount, self.ribocount = [], []
        self.RNApercent_0, self.ribopercent_0 = 0, 0
        self.RNApercent_1, self.ribopercent_1 = 0, 0
        self.RNApercent_5, self.ribopercent_5 = 0, 0
        self.RNApercent_10, self.ribopercent_10 = 0, 0

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
        # global orfcount
        """Looks at all specified BAM files and determine the number of reads
         found in a given region and corrects for length."""
        if self.upcheck and self.downcheck:
            if self.strand:  # is it a (+) strand?
                begin, end = int(self.start), int(self.start) + (int(self.downPos[0]) * 3) + 2
                self.length = end-begin
                self.RNAcount = readCheck(True, int(self.chrom), begin, end, self.length)
                if self.RNAcount is None:
                    pass
                else:
                    self.ribocount = readCheck(False, int(self.chrom), begin, end, self.length)
                    # orfcount += 1

            else:  # is it a (-) strand?
                begin, end = int(self.start - int(self.upPos[-1]) * 3), int(self.start)
                self.length = end-begin
                self.RNAcount = readCheck(True, int(self.chrom), begin, end, self.length)
                if self.RNAcount is None:
                    pass
                else:
                    self.ribocount = readCheck(False, int(self.chrom), begin, end, self.length)
                    # orfcount += 1
        return self

    def metadata(self):
        if self.RNAcount is not None:
            self.RNApercent_0 = sum(value > 0 for value in self.RNAcount)/float(len(self.RNAcount))
            self.ribopercent_0 = sum(value > 0 for value in self.ribocount)/float(len(self.ribocount))
            self.RNApercent_1 = sum(value > 1 for value in self.RNAcount)/float(len(self.RNAcount))
            self.ribopercent_1 = sum(value > 1 for value in self.ribocount)/float(len(self.ribocount))
            self.RNApercent_5 = sum(value > 5 for value in self.RNAcount)/float(len(self.RNAcount))
            self.ribopercent_5 = sum(value > 5 for value in self.ribocount)/float(len(self.ribocount))
            self.RNApercent_10 = sum(value > 10 for value in self.RNAcount)/float(len(self.RNAcount))
            self.ribopercent_10 = sum(value > 10 for value in self.ribocount)/float(len(self.ribocount))
        else:
            pass


# Removes bams from the list to be checked if they are not in the samples of the vcf
def popper(LIST):
    """Iteratively removes entries from LIST that don't match samples in VCF

    Args:
        LIST (list): list of bam files in the given bam directory

    Returns:
        Gives an abbreviated list based on the samples present in the given VCF
        """
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


def modulo_check(LIST):
    mod_check = []
    for i in range(1, 8):
        if len(LIST) % i is 0:
            mod_check.append(i)
        else:
            pass
    return max(mod_check)


# function identifies SNPs, extracts sequence from reference +/-2 nt and looks for start codon within sequence
def ORFSNuper():
    """Identifies SNPs, extracts their sequences from reference +/- 2 nt and looks for start codons."""
    global pre_potORFs
    # global orfcount
    pre_potORFs = []
    global linecount
    linecount = 0
    firstLine = True
    global samples
    with gzip.open(vcf, 'rt')as VCF:
        sampleCheck = False
        # while orfcount <= 49:   # use when debugging
        for line in VCF:
            linecount += 1
            if "##" in line:
                continue
            elif "#CHROM" in line:
                sampleCheck = True
                header = line.split()
                samples = header[9:]
                name = outDir.split('part', 1)[0] + "samples.pkl"
                pickle.dump(samples, open(name, 'wb'))
                popper(RNAsamp_crossref)
                popper(Ribosamp_crossref)
                firstLine = False
            else:
                if not sampleCheck:
                    name = outDir.split('part', 1)[0] + 'samples.pkl'
                    with open(name, 'rb') as f:
                        samples = pickle.load(f)
                if firstLine:
                    popper(RNAsamp_crossref)
                    popper(Ribosamp_crossref)
                    firstLine = False
                columns = line.split()

                # Check to see if it is a SNP
                if len(columns[3]) and len(columns[4]) == 1:

                    # look for the reference sequence around SNP
                    seq = SNP_search(columns[0], int(columns[1]) - 2, int(columns[1]) + 2)
                    seq_step = (seq[:1] + columns[4] + seq[3:]).upper()

                    # iff sequence creates start codon does it make a class instance
                    strand_checker(seq_step, columns, pre_potORFs)
                else:
                    continue
                # For debugging
                # if orfcount >= 50:
                #     print("orfcount met!")
                #     break
    # return pre_potORFs

ORFSNuper()


def dumper(LIST, lst_NUM, OUT):
    dump = LIST
    pickle.dump(dump, open('%s%s%s' % (outDir + "DILL/", str(lst_NUM), OUT), 'wb'))


dumper(pre_potORFs, "", "pre_potORFs.pkl")
# pickle.dump(pre_potORFs, open(outDir + "DILL/" + "pre_potORFs.pkl", 'wb'))
potORFs = np.split(np.asarray(pre_potORFs), modulo_check(pre_potORFs))
potORFs = [sublist.tolist() for sublist in potORFs]


def file_writer(LIST):
    """Outputs a file to outputDir based on the object that is passed to it

    Args:
        LIST (list): one of the potential ORF SNPs located in potORFs generated by ORFSNuPer
    """
    for snp in LIST:
        try:
            filename = '%s%s_%s' % (outDir + "SNPs/", str(snp.chrom), str(snp.SNPpos))
            if snp.RNAcount is None:
                pass
            else:
                if snp.strand:
                    info = ["##" + snp.SNPid, str(snp.chrom), "+", str(snp.start),
                            str((snp.start + int(snp.downPos[0]) * 3)+2)]
                else:
                    info = ["##" + snp.SNPid, str(snp.chrom), "-", str(snp.start),
                            str(snp.start - int(snp.upPos[-1]) * 3)]
                with open(filename, 'w') as writer:
                    header1 = ["##ID", "CHROM", "STRAND", "START", "Nearest_STOP"]
                    print('\t'.join(header1), file=writer)
                    print('\t'.join(info), file=writer)
                    header2 = ["#Sample", "Genotype", "RNA_FPKM", "Ribo_FPKM"]
                    print('\t'.join(header2), file=writer)
                    SNuPed = [(a, b, c, d) for i, (a, b, c, d) in enumerate(zip(samples, snp.genotypes,
                                                                                snp.RNAcount, snp.ribocount))]
                    writer.write('\n'.join('%s\t%s\t%s\t%s' % x for x in SNuPed))
        except TypeError:
            pass


def threader(lst, STEP):
    cmd_lst = ['obj.lookUp()', 'obj.lookDown()', 'obj.WordCount()', 'obj.metadata()']
    pool = ThreadPool()
    exec('pool.map(lambda obj: %s , lst)' % (cmd_lst[STEP]))
    pool.close()
    pool.join()


def threadpool(lst):
    tmp1 = lst
    i = randint(1, 10)
    threader(tmp1, 0)
    tmp1 = [snp for snp in tmp1 if snp.upcheck]
    dumper(tmp1, i, "UPpotORFs.pkl")
    threader(tmp1, 1)
    tmp1 = [snp for snp in tmp1 if snp.downcheck]
    dumper(tmp1, i, "DWNpotORFs.pkl")
    threader(tmp1, 2)
    tmp1 = [snp for snp in tmp1 if snp.RNAcount is not None]
    dumper(tmp1, i, "COUNTpotORFs.pkl")
    threader(tmp1, 3)
    return tmp1


pooler = Pool()
results = pooler.map(threadpool, potORFs)
potORFs = [snp for sublist in results if len(sublist) > 0 for snp in sublist]
dumper(potORFs, "", "potORFs_full.pkl")
file_writer(potORFs)


# Writes a master file containing metadata for each SNP
with open(outDir+"metadata", 'w') as meta:
    writer = csv.writer(meta, delimiter='\t', lineterminator='\n')
    header = ["#CHROM_SNPPOSITION", "%RNA-FPKM>0", "%ribo-FPKM>0", "%RNA-FPKM>1",
              "%ribo-FPKM>1", "%RNA-FPKM>5", "%ribo-FPKM>5", "%RNA-FPKM>10", "%ribo-FPKM>10"]
    writer.writerow(header)
    metadata = [[str(snp.chrom) + "_" + str(snp.SNPpos),
                 str(snp.RNApercent_0), str(snp.ribopercent_0),
                 str(snp.RNApercent_1), str(snp.ribopercent_1),
                 str(snp.RNApercent_5), str(snp.ribopercent_5),
                 str(snp.RNApercent_10), str(snp.ribopercent_10)] for snp in potORFs if snp.RNAcount is not None]
    writer.writerows(metadata)

# find out how long the process took
endTime, endasc = time.time(), time.asctime()
m, s = divmod(endTime - startTime, 60)
h, m = divmod(m, 60)
d, h = divmod(h, 24)

# Write a small report for start time, end time, and elapsed time
with open(outDir + "out.log", 'w') as f:
    print(str(len(potORFs)) + " potential ORFs found", file=f)
    print("Program started:", file=f)
    print(startasc, file=f)
    print('', file=f)
    print("Program completed:", file=f)
    print(endasc)
    print('', file=f)
    print("Elapsed time:", file=f)
    print("%s days, %s hours, %s minutes, %s seconds" % (str(d), str(h), str(m), str(round(s, 2))), file=f)
    print('', file=f)

pickle.dump(potORFs, open(outDir + "DILL/" + "potORFs.pkl", 'wb'))
