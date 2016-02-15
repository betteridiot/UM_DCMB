# -*- coding: utf-8 -*-
# !/usr/bin/env python

import os
import gzip
import glob
import time
import argparse
import math
import csv
import numpy as np
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

# TODO add for loop implementation of alternative start codons
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
os.chdir(RNADir)
for dir, _, _ in os.walk(os.getcwd()):
    RNAbams.extend(glob.glob(os.path.join(dir, "*hits.bam")))
os.chdir(riboDir)
for dir, _, _ in os.walk(os.getcwd()):
    Ribobams.extend(glob.glob(os.path.join(dir, "*sort.bam")))

# TODO link sampleFinder to read count


# Gets a list of all the bam files that have the given genotype
def sampleFinder(LIST, RNAorRibo):  # True for RNA, False for Ribosome
    templist1, templist2 = [], []
    if RNAorRibo:
        bamlist = RNAbams
        templist2.extend([bamlist[i] for i in range(len(bamlist)) for line in LIST if line in bamlist[i]])
    elif not RNAorRibo:
        bamlist = Ribobams
        templist1.extend([line[0] for line in ribosamples for element in LIST if element in line[1]])
        templist2.extend([bamlist[i] for i in range(len(bamlist)) for line in templist1 if line in bamlist[i]])
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
    def __init__(self, CHROM, START, END, ID, STRAND,
                 HOMO_ref, HOMO_SNP, HETERO, REFSAMP, SNPSAMP, HETSAMP):
        self.chrom = CHROM
        self.start = START
        self.end = END
        self.SNPid = ID
        self.strand = STRAND
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
        self.RNAcount, self.ribocount = None, None
        self.homo_ref = HOMO_ref
        self.homo_SNP = HOMO_SNP
        self.heterozygous = HETERO
        self.length = None
        self.homorefsamp, self.homosnpsamp, self.hetsamp = REFSAMP, SNPSAMP, HETSAMP
        # Produce a list of bams that present with a given genotype for a given SNP
        self.RNArefsamp = sampleFinder(self.homorefsamp, True)
        self.RNAsnpsamp = sampleFinder(self.homosnpsamp, True)
        self.RNAhetsamp = sampleFinder(self.hetsamp, True)
        self.RIBOrefsamp = sampleFinder(self.homorefsamp, False)
        self.RIBOsnpsamp = sampleFinder(self.homosnpsamp, False)
        self.RIBOhetsamp = sampleFinder(self.hetsamp, False)
        self.RNArefcount, self.RNAsnpcount, self.RNAhetcount = [], [], []
        self.RIBOrefcount, self.RIBOsnpcount, self.RIBOhetcount = [], [], []

    # check to see if a stop codon is within upstream sequence (downstream if (-) strand)
    def lookUp(self):
        codonCount = 0
        if self.strand:
            stops = plusStops
        else:
            stops = negStops
        for y in range(0, threshold, 3):
            codon = self.up[y:y + 3]
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
        for x in range(0, threshold, 3):
            codon = self.down[x:x + 3]
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
                if self.RNAcount == (None or float(0)):
                    pass
                else:
                    self.ribocount = readCheck(False, int(self.chrom), begin, end, self.length)

            else:  # is it a (-) strand?
                begin, end = int(self.start - int(self.upPos[-1]) * 3), int(self.start)
                self.length = end-begin
                self.RNAcount = readCheck(True, int(self.chrom), begin, end, self.length)
                self.RNArefcount = genoCheck(self.RNArefsamp, int(self.chrom), begin, end, self.length)
                self.RNAsnpcount = genoCheck(self.RNAsnpsamp, int(self.chrom), begin, end, self.length)
                self.RNAhetcount = genoCheck(self.RNAhetsamp, int(self.chrom), begin, end, self.length)
                if self.RNAcount == (None or float(0)):
                    pass
                else:
                    self.ribocount = readCheck(False, int(self.chrom), begin, end, self.length)
                    self.RIBOrefcount = genoCheck(self.RIBOrefsamp, int(self.chrom), begin, end, self.length)
                    self.RIBOsnpcount = genoCheck(self.RIBOsnpsamp, int(self.chrom), begin, end, self.length)
                    self.RIBOhetcount = genoCheck(self.RIBOhetsamp, int(self.chrom), begin, end, self.length)




# Used to iterate potential ORF class instantiation
def portORF(CHROM, START, END, ID, STRAND, HOREF, HOSNP, HETERO, REFSAMP, SNPSAMP, HETSAMP):
    portedORF = potORF(CHROM, START, END, ID, STRAND, HOREF, HOSNP, HETERO,
                       REFSAMP, SNPSAMP, HETSAMP)
    return portedORF


# iteratively pulls FPKM over a given region across all BAMs
def readCheck(RNAorRIBO, CHROM, START, STOP, LENGTH):
    bamlist, WC = [], []
    typeCheck = RNAorRIBO
    if typeCheck:  # True = RNA-seq, False = Ribosome profiling
        bamlist = RNAbams
    elif typeCheck is False:
        bamlist = Ribobams
    for entry in bamlist:
        readcount = os.popen('samtools view -q 10 ' + entry + ' chr%d:%d-%d | wc -l'
                             % (int(CHROM), int(START), int(STOP)))
        count = float(readcount.readline().rstrip())
        fullcounter = os.popen("samtools idxstats " + entry + " |awk '{sum+=$3} END {print sum}'")
        fullcount = float(fullcounter.readline().rstrip())
        WC.extend([round((count/(LENGTH*fullcount))*math.pow(10, 9), 4)])
    if sum(WC) == float(0):
        return "NA"
    else:
        return round(float(sum(WC)) / int(len(WC)), 3)

def genoCheck(DIR, CHROM, START, STOP, LENGTH):
    bamlist, WC = [], []
    bamlist = DIR
    for entry in bamlist:
        readcount = os.popen('samtools view -q 10 ' + entry + ' chr%d:%d-%d | wc -l'
                             % (int(CHROM), int(START), int(STOP)))
        count = float(readcount.readline().rstrip())
        fullcounter = os.popen("samtools idxstats " + entry + " |awk '{sum+=$3} END {print sum}'")
        fullcount = float(fullcounter.readline().rstrip())
        WC.extend([round((count/(LENGTH*fullcount))*math.pow(10, 9), 4)])
    return WC


# function identifies SNPs, extracts sequence from reference +/-2 nt and looks for start codon within sequence
def ORFSNuper():
    global potORFs
    global samplenames
    potORFs = []
    # orfcount = 0  # use when debugging
    with gzip.open(vcf, 'rt')as VCF:
        # while orfcount < 15:   # use when debugging
        for line in VCF:
            # skip all of the lines before content
            if "##" in line:
                continue
            elif "#CHROM" in line:
                header = line.split()
                samplenames = header[9:]
            else:
                columns = line.split()

                # Pull out genotype and sample name information for each SNP
                heterlist = np.asarray([z for z, genotype in enumerate(columns) if genotype in ["1|0" or "0|1"]])
                heter = [header[index] for index in heterlist]
                horeflist = np.asarray([z for z, genotype in enumerate(columns) if genotype == "0|0"])
                horef = [header[index] for index in horeflist]
                hosnplist = np.asarray([z for z, genotype in enumerate(columns) if genotype == "1|1"])
                hosnp = [header[index] for index in hosnplist]

                # Pull out genotype counts for a given row in VCF
                homozygous_ref = columns.count("0|0")
                homozygous_SNP = columns.count("1|1")
                heterozygous = columns.count("1|0")+columns.count("0|1")

                # Check to see if it is a SNP
                if len(columns[3]) and len(columns[4]) == 1:

                    # look for the reference sequence around SNP
                    seq = os.popen('samtools faidx %s/chr%s.fa chr%s:%d-%d' % (
                        reference, columns[0], columns[0], int(columns[1]) - 2, int(columns[1]) + 2))
                    seq.readline()
                    seq = seq.read().rstrip()
                    seq_step = (seq[:1] + columns[4] + seq[3:]).upper()

                    # Check to see if (+) strand ORF is found
                    if plusStart in seq_step:
                        posCheck = str.find(seq_step, plusStart) + 1
                        if 1 <= posCheck < 3:
                            seqPos = int(columns[1]) - posCheck
                        elif posCheck > 3:
                            seqPos = int(columns[1]) + posCheck
                        else:
                            pass

                        # Create potential ORF class instance
                        potORFs.extend([portORF(columns[0], seqPos, seqPos + 2, columns[2], True,
                                                homozygous_ref, homozygous_SNP, heterozygous,
                                                horef, hosnp, heter)])
                        # orfcount += 1  # use when debugging

                    # Check to see if (-) strand ORF is found
                    if negStart in seq_step:
                        posCheck = str.find(seq_step, negStart) + 1
                        if 1 <= posCheck < 3:
                            seqPos = int(columns[1]) + posCheck
                        elif posCheck > 3:
                            seqPos = int(columns[1]) - posCheck
                        else:
                            pass

                        # Create potential ORF class instance
                        potORFs.extend([portORF(columns[0], seqPos, seqPos - 2, columns[2], False,
                                                homozygous_ref, homozygous_SNP, heterozygous,
                                                horef, hosnp, heter)])
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
                                          str(potORFs[i].length), str(potORFs[i].homo_ref),
                                          str(potORFs[i].homo_snp), str(potORFs[i].heterozygous)])])
            else:
                SNuPed.extend(['\t'.join([str(potORFs[i].chrom), "-", str(potORFs[i].start),
                                          str(potORFs[i].start - int(potORFs[i].upPos[-1]) * 3),
                                          str(potORFs[i].RNAcount), str(potORFs[i].ribocount),
                                          str(potORFs[i].length), str(potORFs[i].homo_ref),
                                          str(potORFs[i].homo_snp), str(potORFs[i].heterozygous)])])
    else:
        continue

# Creates a list of the SNP and a sample's FPKM given its genotype
SNuPed_RNAhoref, SNuPed_RNAhosnp, SNuPed_RNAhet = [], [], []
SNuPed_RIBOhoref, SNuPed_RIBOhosnp, SNuPed_RIBOhet = [], [], []
for i in range(len(potORFs)):
    if potORFs[i].upcheck and potORFs[i].downcheck:
        if potORFs[i].strand:
            SNuPed_RNAhoref.extend(['\t'.join([str(potORFs[i].SNPid), str(potORFs[i].chrom), str(potORFs[i].start),
                                       str((potORFs[i].start + int(potORFs[i].downPos[0]) * 3)+2)])+"\t"+'\t'.join(potORFs[i].RNArefcount)])
            SNuPed_RNAhosnp.extend(['\t'.join([str(potORFs[i].SNPid), str(potORFs[i].chrom), str(potORFs[i].start),
                                       str((potORFs[i].start + int(potORFs[i].downPos[0]) * 3)+2)])+"\t"+'\t'.join(potORFs[i].RNAsnpcount)])
            SNuPed_RNAhet.extend(['\t'.join([str(potORFs[i].SNPid), str(potORFs[i].chrom), str(potORFs[i].start),
                                       str((potORFs[i].start + int(potORFs[i].downPos[0]) * 3)+2)])+"\t"+'\t'.join(potORFs[i].RNAhetcount)])
            SNuPed_RIBOhoref.extend(['\t'.join([str(potORFs[i].SNPid), str(potORFs[i].chrom), str(potORFs[i].start),
                                       str((potORFs[i].start + int(potORFs[i].downPos[0]) * 3)+2)])+"\t"+'\t'.join(potORFs[i].RIBOrefcount)])
            SNuPed_RIBOhosnp.extend(['\t'.join([str(potORFs[i].SNPid), str(potORFs[i].chrom), str(potORFs[i].start),
                                       str((potORFs[i].start + int(potORFs[i].downPos[0]) * 3)+2)])+"\t"+'\t'.join(potORFs[i].RIBOsnpcount)])
            SNuPed_RIBOhet.extend(['\t'.join([str(potORFs[i].SNPid), str(potORFs[i].chrom), str(potORFs[i].start),
                                       str((potORFs[i].start + int(potORFs[i].downPos[0]) * 3)+2)])+"\t"+'\t'.join(potORFs[i].RIBOhetcount)])
        else:
            SNuPed_RNAhoref.extend(['\t'.join([str(potORFs[i].SNPid), str(potORFs[i].chrom), str(potORFs[i].start),
                                       str(potORFs[i].start - int(potORFs[i].upPos[-1]) * 3)])+"\t"+'\t'.join(potORFs[i].RNArefcount)])
            SNuPed_RNAhosnp.extend(['\t'.join([str(potORFs[i].SNPid), str(potORFs[i].chrom), str(potORFs[i].start),
                                       str(potORFs[i].start - int(potORFs[i].upPos[-1]) * 3)])+"\t"+'\t'.join(potORFs[i].RNAsnpcount)])
            SNuPed_RNAhet.extend(['\t'.join([str(potORFs[i].SNPid), str(potORFs[i].chrom), str(potORFs[i].start),
                                       str(potORFs[i].start - int(potORFs[i].upPos[-1]) * 3)])+"\t"+'\t'.join(potORFs[i].RNAhetcount)])
            SNuPed_RIBOhoref.extend(['\t'.join([str(potORFs[i].SNPid), str(potORFs[i].chrom), str(potORFs[i].start),
                                       str(potORFs[i].start - int(potORFs[i].upPos[-1]) * 3)])+"\t"+'\t'.join(potORFs[i].RIBOrefcount)])
            SNuPed_RIBOhosnp.extend(['\t'.join([str(potORFs[i].SNPid), str(potORFs[i].chrom), str(potORFs[i].start),
                                       str(potORFs[i].start - int(potORFs[i].upPos[-1]) * 3)])+"\t"+'\t'.join(potORFs[i].RIBOsnpcount)])
            SNuPed_RIBOhet.extend(['\t'.join([str(potORFs[i].SNPid), str(potORFs[i].chrom), str(potORFs[i].start),
                                       str(potORFs[i].start - int(potORFs[i].upPos[-1]) * 3)])+"\t"+'\t'.join(potORFs[i].RIBOhetcount)])


# find out how long the process took
endTime, endasc = time.time(), time.asctime()
m, s = divmod(endTime - startTime, 60)
h, m = divmod(m, 60)
d, h = divmod(h, 24)

# Print the list of potential ORFs in a tab-delimited file
with open(outfile, 'w') as f:
    print >> f, "Sequencing read counts normalized by FPKM"
    print >> f, "Genotype counts are sums across all samples in VCF"
    print >> f, "\t".join(["CHROM", "STRAND", "START", "Nearest_STOP", "RNA_ReadCount",
                           "Ribo_ReadCount", "ORF_Length", "0|0", "0|1", "1|1"])
    print >> f, "\n".join(SNuPed)

# Writes a tab-delimited file of the SNP versus genotype FPKM list
outfilelist = [SNuPed_RNAhoref, SNuPed_RNAhosnp, SNuPed_RNAhet, SNuPed_RIBOhoref, SNuPed_RIBOhosnp, SNuPed_RIBOhet]
outfileext = ["_RNAhoref", "_RNAhosnp", "_RNAhet", "_RIBOhoref", "_RIBOhosnp", "_RIBOhet"]
for entry, ext in outfilelist, outfileext:
    with open(outfile+outfileext[ext], 'w') as g:
        print >> g, "Sequencing read counts normalized by FPKM"
        print >> g, '\t'.join(["SNP_ID", "CHROM", "START", "STOP"])+"\t"+'\t'.join(samplenames)
        print >> g, '\n'.join(entry)


# Write a small report for start time, end time, and elapsed time
with open(outfile + ".log", 'w') as f:
    print >> f, "Program started:"
    print >> f, startasc
    print >> f, ""
    print >> f, "Program completed:"
    print >> f, endasc
    print >> f, ""
    print >> f, "Elapsed time:"
    print >> f, str(d) + " days", str(h) + " hours", str(m) + " minutes", str(round(s, 2)) + " seconds"