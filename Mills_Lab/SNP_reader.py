#Psuedocode
#glob all of the SNP files
#iteratively open and read in as tab delimited
#skip lines with pound signs
#find the ones that have both rna-seq and ribosome profiling data

from __future__ import print_function
import os
import glob
import matplotlib as mpl

os.chdir('/home/mdsherm/Project/SNuPer_results/pythonTest/100ktest')
snp_files = glob.glob(os.path.join(os.getcwd(), '*.snp'))