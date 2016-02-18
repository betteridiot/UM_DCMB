# !/usr/bin/env python

from __future__ import print_function
import os
import gzip
import glob
import time
import argparse
import math
import csv
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
parser.add_argument('-o', action='store', dest='output', help='Set output filename',
                    default='/home/mdsherm/Project/SNuPer_results/pythonTest')
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
outfile = args.output
reference = args.ref
threshold = args.threshold
# ARGPARSE END

# Adding in ribosome sample name to SRA ID file
ribosamples = []
with open('/home/mdsherm/Project/ribosamplescross') as ribo:
    reader = csv.reader(ribo, delimiter='\t')
    ribosamples.extend([row for row in reader])


# CODONS START
negStops = ['TTA', 'CTA', 'TCA']
plusStops = ['ATT', 'ATC', 'ACT']
if args.alternative is True:
    plusStart, negStart = ['TAC', 'CAC'], ['CAT', 'CAC']
else:
    plusStart = "TAC"
    negStart = "CAT"
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


# Finds total reads per bam file for FPKM
def read_indexer(fileDIR, LIST):
    for item in fileDIR:
        fullcounter = os.popen("samtools idxstats " + item + " |awk '{sum+=$3} END {print sum}'")
        fullcount = float(fullcounter.readline().rstrip())
        LIST.extend([[item, fullcount]])


read_indexer(RNAbams, RNAbams_total)
read_indexer(Ribobams, Ribobams_total)

# TODO may not need
# Gets a list of all the bam files that have the given genotype
def sampleFinder(LIST, RNAorRibo):  # True for RNA, False for Ribosome
    templist1, templist2 = [], []
    if RNAorRibo:
        bamlist = RNAbams
        templist2.extend([bamlist[j] for j in range(len(bamlist)) for line in LIST if line in bamlist[j]])
    elif not RNAorRibo:
        bamlist = Ribobams
        templist1.extend([line[0] for line in ribosamples for element in LIST if element in line[1]])
        templist2.extend([bamlist[n] for n in range(len(bamlist)) for line in templist1 if line in bamlist[n]])
    return templist2


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
    def __init__(self, CHROM, START, END, ID, STRAND):
        self.chrom = CHROM
        self.start = START
        self.end = END
        self.SNPid = ID
        self.strand = STRAND
        self.length = 0
        self.up = os.popen('samtools faidx %s/chr%s.fa chr%s:%d-%d' % (
            reference, self.chrom, self.chrom, int(self.start) - threshold, int(self.start) - 1))
        self.up.readline()
        self.up = ((self.up.read()).rstrip()).upper()
        self.down = os.popen('samtools faidx %s/chr%s.fa chr%s:%d-%d' % (
            reference, self.chrom, self.chrom, int(self.end) + 1, int(self.end) + threshold))
        self.down.readline()
        self.down = ((self.down.read()).rstrip()).upper()
        self.upcheck, self.downcheck = False, False
        self.upPos, self.downPos = [], []
        # self.RNAcount, self.ribocount = None, None

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
                self.upPos.extend([codonCount, ])
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
                self.downPos.extend([codonCount, ])
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


# Used to iterate potential ORF class instantiation
def portORF(CHROM, START, END, ID, STRAND):
    portedORF = potORF(CHROM, START, END, ID, STRAND)
    return portedORF


# iteratively pulls FPKM over a given region across all BAMs
def readCheck(RNAorRIBO, CHROM, START, STOP, LENGTH):
    bamlist, WC = [], []
    typeCheck = RNAorRIBO
    if typeCheck:  # True = RNA-seq, False = Ribosome profiling
        bamlist = RNAbams
    elif typeCheck is False:
        bamlist = Ribobams
    for step in bamlist:
        readcount = os.popen('samtools view -q 10 ' + step + ' chr%d:%d-%d | wc -l'
                             % (int(CHROM), int(START), int(STOP)))
        count = float(readcount.readline().rstrip())

        # TODO tie FPKM of sample to total read count of sample
        fullcounter = os.popen("samtools idxstats " + step + " |awk '{sum+=$3} END {print sum}'")
        fullcount = float(fullcounter.readline().rstrip())
        WC.extend([round((count/(LENGTH*fullcount))*math.pow(10, 9), 4)])
    if sum(WC) == float(0):
        return "NA"
    else:
        return round(float(sum(WC)) / int(len(WC)), 3)


def genoCheck(DIR, CHROM, START, STOP, LENGTH):
    bamlist, WC = [], []
    bamlist = DIR
    for level in bamlist:
        readcount = os.popen('samtools view -q 10 ' + level[0] + ' chr%d:%d-%d | wc -l'
                             % (int(CHROM), int(START), int(STOP)))
        count = float(readcount.readline().rstrip())
        WC.extend([round((count/(LENGTH*level[1]))*math.pow(10, 9), 4)])
    return WC


# Pull out the reference sequence at a given position
def SNP_search(CHROM, START, STOP):
    seq = os.popen('samtools faidx %s/chr%s.fa chr%s:%d-%d' % (reference, CHROM, CHROM, START, STOP))
    seq.readline()
    seq = seq.read().rstrip()
    return seq


# function identifies SNPs, extracts sequence from reference +/-2 nt and looks for start codon within sequence
def ORFSNuper():
    global potORFs
    global samplenames
    potORFs = []
    # orfcount = 0  # use when debugging
    with gzip.open(vcf, 'rt')as VCF:
        # while orfcount < 15:   # use when debugging
        for line in VCF:
            if "##" in line:
                continue
            elif "#CHROM" in line:
                header = line.split()
                samplenames = header[9:]
            else:
                columns = line.split()

                # Check to see if it is a SNP
                if len(columns[3]) and len(columns[4]) == 1:

                    # look for the reference sequence around SNP
                    seq = SNP_search(columns[0], int(columns[1]) - 2, int(columns[1]) + 2)
                    seq_step = (seq[:1] + columns[4] + seq[3:]).upper()
                    # Check to see if (+) strand ORF is found
                    if plusStart in seq_step:
                        posCheck = seq_step.find(plusStart) + 1
                        if 1 <= posCheck < 3:
                            seqPos = int(columns[1]) - posCheck
                        elif posCheck > 3:
                            seqPos = int(columns[1]) + posCheck
                        else:
                            pass
                        # Create potential ORF class instance
                        potORFs.extend([portORF(columns[0], seqPos, seqPos + 2, columns[2], True)])
                        # orfcount += 1  # use when debugging

                    # Check to see if (-) strand ORF is found
                    if negStart in seq_step:
                        posCheck = seq_step.find(seq_step, negStart) + 1
                        if 1 <= posCheck < 3:
                            seqPos = int(columns[1]) + posCheck
                        elif posCheck > 3:
                            seqPos = int(columns[1]) - posCheck
                        else:
                            pass
                        # Create potential ORF class instance
                        potORFs.extend([portORF(columns[0], seqPos, seqPos - 2, columns[2], False)])
                        # orfcount += 1  # use when debugging
                else:
                    continue
            # For debugging
            # if orfcount >= 15:
            #     print("orfcount met!")
            #     break


# Find the potential ORFs
ORFSNuper()

# Look upstream and downstream for stop codons and read count of ribosome/RNA bam files by class instance multithreading
pool = ThreadPool()
pool.map(lambda obj: obj.lookUp().lookDown().WordCount(), potORFs)
pool.close()
pool.join()


# TODO think of writing file iteratively
# Using the joined instances of potential ORFs, cleans & coalesces the data for output
SNuPed = []
for i in range(len(potORFs)):
    if potORFs[i].upcheck and potORFs[i].downcheck:  # does it have an up/downstream stop?
        if potORFs[i].RNAcount == (0 or None):
            continue
        else:
            if potORFs[i].strand:  # is it a (+) strand?
                # if there were RNA-seq reads, check for Ribosome profiling reads (translation)
                SNuPed.extend(['\t'.join([str(potORFs[i].chrom), "+", str(potORFs[i].start),
                                          str((potORFs[i].start + int(potORFs[i].downPos[0]) * 3)+2),
                                          str(potORFs[i].RNAcount), str(potORFs[i].ribocount),
                                          str(potORFs[i].length)])])
            else:
                SNuPed.extend(['\t'.join([str(potORFs[i].chrom), "-", str(potORFs[i].start),
                                          str(potORFs[i].start - int(potORFs[i].upPos[-1]) * 3),
                                          str(potORFs[i].RNAcount), str(potORFs[i].ribocount),
                                          str(potORFs[i].length)])])
    else:
        continue

# find out how long the process took
endTime, endasc = time.time(), time.asctime()
m, s = divmod(endTime - startTime, 60)
h, m = divmod(m, 60)
d, h = divmod(h, 24)

# Print the list of potential ORFs in a tab-delimited file
with open(outfile, 'w') as f:
    f.write("Sequencing read counts normalized by FPKM")
    print("Genotype counts are sums across all samples in VCF", file=f)
    print("\t".join(["CHROM", "STRAND", "START", "Nearest_STOP", "RNA_ReadCount",
                     "Ribo_ReadCount", "ORF_Length", "0|0", "0|1", "1|1"]), file=f)
    # print >> f, "\n".join(SNuPed)

# Writes a tab-delimited file of the SNP versus genotype FPKM list
# outfilelist = [SNuPed_RNAhoref, SNuPed_RNAhosnp, SNuPed_RNAhet, SNuPed_RIBOhoref, SNuPed_RIBOhosnp, SNuPed_RIBOhet]
# outfileext = ["_RNAhoref", "_RNAhosnp", "_RNAhet", "_RIBOhoref", "_RIBOhosnp", "_RIBOhet"]
# for entry, ext in outfilelist, outfileext:
#     with open(outfile+outfileext[ext], 'w') as g:
#         print >> g, "Sequencing read counts normalized by FPKM"
#         print >> g, '\t'.join(["SNP_ID", "CHROM", "STRAND", "START", "STOP"])+"\t"+'\t'.join(samplenames)
#         print >> g, '\n'.join(entry)
#

# Write a small report for start time, end time, and elapsed time
with open(outfile + ".log", 'w') as f:
    print("Program started:", file=f)
    print(startasc, file=f)
    print('', file=f)
    print("Program completed:", file=f)
    print(endasc)
    print('', file=f)
    print("Elapsed time:", file=f)
    print("%s days, %s hours, %s seconds" % (str(h), str(m), str(round(s, 2))), file=f)
    print('', file=f)
