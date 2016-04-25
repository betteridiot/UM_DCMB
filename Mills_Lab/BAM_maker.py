#!/usr/bin/env python
import os
import glob
import csv
import dill as pickle

RNADir = '/home/mdsherm/Rotation/RNA_fq/tophat_hg19'
riboDir = '/home/mdsherm/Rotation/ribosomal/bwa_alignment'
outDir = '/home/mdsherm/Project/UM_DCMB/Mills_Lab/pkl/'

# Adding in ribosome sample name to SRA ID file
ribosamples = []
with open('/home/mdsherm/Project/ribosamplescross') as ribo:
    reader = csv.reader(ribo, delimiter='\t')
    ribosamples.extend([row for row in reader])

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

pickle.dump(RNAsamp_crossref, open(outDir + "RNAsamp_crossref.pkl", 'wb'))
pickle.dump(Ribosamp_crossref, open(outDir + "Ribosamp_crossref.pkl", 'wb'))
pickle.dump(RNAbams_dict, open(outDir + "RNAbams_dict.pkl", 'wb'))
pickle.dump(RNAtotalreads, open(outDir + "RNAtotalreads.pkl", 'wb'))
pickle.dump(Ribobams_dict, open(outDir + "Ribobams_dict.pkl", 'wb'))
pickle.dump(Ribototalreads, open(outDir + "Ribototalreads.pkl", 'wb'))