# This is a test of the GitHub System
#This version is for vcf file of 1000 genome hg37
import sys
import os
import gzip
import re

infile1, infile2, outfile2, outfile3 = sys.argv[1:]
#/mnt/EXT/Mills-scratch/datasets/1000genomes/ftp/release/20130502/ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz; reference file, /mnt/EXT/Mills-data/yifwang/human_reference/hg37_ucsc; #upstream of existing CAT; Inside an existing CAT; New One (respectively)


outf1 = open(outfile1,'w')
print >> outf1,"\t".join(["New_created_CAT","Existing_Downstream_CAT"])
outf2 = open(outfile2,'w')
print >> outf2, "\t".join(["New_created_CAT","Existing_Upstream_CAT","Stop_Codon_Downstream"])
outf3 = open(outfile3,'w')
print >> outf3, "\t".join(["New_created_CAT","Stop_Codon_Upstream","Stop_Codon_Downstream"])

inf = gzip.open(infile1)



def find_seq(pos,chrm):
    start = pos-2
    end = pos+2
    seq_prev = os.popen('samtools faidx %s/%s.fa %s:%d-%d' %(infile2,chrm,chrm,start,end))
    seq = seq_prev.read().split('\n')[1]
    return seq

#whether there is an CAT to non or non to CAT
def judge_start_change(sequence):
    pos_seq = []
    for i in range(0,3):
        codon=sequence[i:i+3]
        if codon == "CAT":
            pos_seq.append(i)
        else:
            continue
    return pos_seq
#Actually looking for downstream here
def look_for_signal_upstream(start_original,j,chrm):
    start = start_original
    end = start_original + 2999
    seq_read = os.popen('samtools faidx %s/%s.fa %s:%d-%d' %(infile2,chrm,chrm,start,end))
    seq_read.readline()
    seq_read = seq_read.read()
    seq = ""
    for seq_line in seq_read:
        seq += seq_line.rstrip()
    seq = seq.upper()

    for j in range(0,3000,3):
        sub_seq= seq[j:j+3]
        #print sub_seq
        if sub_seq == "CAT":
            return [1,start_original+j+2]

        elif sub_seq == "CTA" or  sub_seq == "TTA" or sub_seq == "TCA":
            ####print start
            return [2,start_original+j+2]
    return None

#It is upstream now for minus strand
def look_for_signal_downstream(start_pos,i,chrm):

    ####print seq
    ####print start
    #There is a stop codon before any existing "CAT"
    #Look for start or stop upstream of the new "CAT"
    start = start_pos - 3002
    end = start_pos - 3
    seq_read = os.popen('samtools faidx %s/%s.fa %s:%d-%d' %(infile2,chrm,chrm,start,end))
    seq_read.readline()
    seq_read = seq_read.read()
    seq = ""
    for seq_line in seq_read:
        seq += seq_line.rstrip()
    #reverse the sequence since it is on minus strand
    seq = seq[::-1]
    seq = seq.upper()
    #print seq
    for i in range(0,3000,3):
        sub_seq= seq[i:i+3]
        #print sub_seq
        if sub_seq == "ATC" or  sub_seq == "ATT" or sub_seq == "ACT":
            start_original = start_pos + 1
            result_li = look_for_signal_upstream(start_original,0,chrm)
            if result_li == None:
                continue
            else:
                if result_li[0] == 1:
                    return 1,result_li[1],start_pos-5-i
                elif result_li[0] == 2:
                    return 2,result_li[1],start_pos-5-i
        #sicne the seuqnce is already reversed, we only need to take the complementary of ATG here as a start codon
        elif sub_seq == "TAC":
            return 0,start_pos-5-i,0
        else:
            continue
    #if there is nothing in +- 3000, return nothing
    return None


for line in inf:
    if '#' in line:
        continue
    else:
        #print line
        pos = int(line.rstrip().split("\t")[1])
        chrm = "chr" + line.rstrip().split("\t")[0]
        sequence = find_seq(pos,chrm)
        alt_base_1 = line.rstrip().split("\t")[3]
        alt_base_2 = line.rstrip().split("\t")[4]

        if len(alt_base_1) == 1 & len(alt_base_2) ==1:
            seq1 = sequence.replace(sequence[2],alt_base_1)
            seq2 = sequence.replace(sequence[2],alt_base_2)

            pos_seq_1 = judge_start_change(seq1)
            pos_seq_2 = judge_start_change(seq2)

            if len(pos_seq_1) != 0:
                for j in pos_seq_1:
                    start_pos = pos + j
                    result_li = look_for_signal_downstream(start_pos,0,chrm)
                    if result_li == None:
                        continue
                    else:
                        #XXXXXXXXX"CAT"(newly created)XXXXXXXXXCATXXXXXXXXXATT
                        if result_li[0] == 0:
                           print >> outf1, "\t".join([chrm,str(start_pos),str(result_li[1]-2)])
                        # CATXXXXXXXXXXXX"CAT"(newly created)XXXXXXXXXXXXATT
                        elif result_li[0] == 1:
                            print >> outf2, "\t".join([chrm,str(start_pos),str(result_li[1]-2),str(result_li[2])])
                        #ATTXXXXXXXXX"CAT"(newly created)XXXXXXXXXXXXXXXATT
                        elif result_li[0] == 2:
                            print >> outf3, "\t".join([chrm,str(start_pos),str(result_li[1]-2),str(result_li[2])])
            else:
                continue

            if len(pos_seq_2) != 0:
                for j in pos_seq_2:
                    start_pos = pos + j
                    result_li = look_for_signal_downstream(start_pos,0,chrm)
                    if result_li == None:
                        continue
                    else:
                        #XXXXXXXXX"CAT"(newly created)XXXXXXXXXCATXXXXXXXXXATT
                        if result_li[0] == 0:
                           print >> outf1, "\t".join([chrm,str(start_pos),str(result_li[1]-2)])
                        # CATXXXXXXXXXXXX"CAT"(newly created)XXXXXXXXXXXXATT
                        elif result_li[0] == 1:
                            print >> outf2, "\t".join([chrm,str(start_pos),str(result_li[1]-2),str(result_li[2])])
                        #ATTXXXXXXXXX"CAT"(newly created)XXXXXXXXXXXXXXXATT
                        elif result_li[0] == 2:
                            print >> outf3, "\t".join([chrm,str(start_pos),str(result_li[1]-2),str(result_li[2])])
            else:
                continue

inf.close()
outf1.close()
outf2.close()
outf3.close()