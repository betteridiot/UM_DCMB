#!/usr/bin/env python
from __future__ import print_function
import argparse
import csv
import glob
from operator import itemgetter
import os
from tempfile import TemporaryFile

# ARGPARSE START
parser = argparse.ArgumentParser(description='Concatenate meta files for a chromosome')
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-l', action='store', dest='logs', help='Path/to/logfile/parent/', default=False)
group.add_argument('-m', action='store', dest='meta', help='Path/to/metadata/parent', default=False)
args = parser.parse_args()
if not args.logs:
    parent = args.meta
    search_term = "metadata"
else:
    parent = args.logs
    search_term = "out.log"


def globber():
    """Walks the directory and globs all of the *.log files

    Returns:
        A list of *.log files for a given chromosome.
    """
    cwd = parent
    os.chdir(cwd)
    globs = []
    for directory, _, _ in os.walk(os.getcwd()):
        globs.extend(glob.glob(os.path.join(directory, search_term)))
    return globs


def logger(lst):
    """Returns a sum of potential ORFs found for a given chromosome.

    Args:
        lst (list): list of *.log files from all ORFSNuPer runs

    Returns:
        Sum of potential ORFs for a given chromosome.
    """
    num_orfs = []
    for entry in lst:
        with open(entry, 'rb') as f:
            reader = f.readline()
            reader = (int(string) for string in reader.split() if string.isdigit() and string > 0)
        num_orfs.extend(reader)
    return str(sum(num_orfs))


def meta_join(lst):
    """Writes a 'meta'-metafile to chromosome directory taken from the different
    ORFSNuPer runs.

    Args:
        lst (list): list of metadata files for a given chromosome
    """
    first = True
    line_count = 0
    with TemporaryFile('a+b') as outfile:
        for entry in lst:
            with open(entry, "r+b") as infile:
                if first:
                    line_count + 1
                    header = infile.readline().strip().split()
                    for line in infile:
                        line_count += 1
                        outfile.write(line)
                    first = False
                else:
                    next(infile)
                    for line in infile:
                        line_count += 1
                        outfile.write(line)
        outfile.seek(0)
        table = [line.rstrip().split() for line in outfile]
        dup_check = set()
        unique = {x[0]: x[1] for x in table if x[0] not in dup_check and not dup_check.add(x[0])}
        dupes = [x for x in table if unique.has_key(x[0]) and unique.get(x[0]) is not x[1]]
        print(str(len(table) - len(dup_check)) +
              " duplicate SNP positions (possibly different alleles).")
        with open(parent + "duplicate_SNPs", 'wb') as dup:
            duper = csv.writer(dup, delimiter='\t', lineterminator='\n')
            note = "# Duplicate SNP postion entries (possible different alleles)."
            duper.writerow([note])
            duper.writerow(header)
            duper.writerows(dupes)
        table.sort(key=itemgetter(0, 1))
        with open(parent + "meta_full", 'wb') as final:
            writer = csv.writer(final, delimiter='\t', lineterminator='\n')
            writer.writerow(header)
            writer.writerows(table)


def main():
    """Takes the glob and iteratively pulls out number of potential ORFs
    or writes a concatenates a 'meta'-metafile to chromosome directory

    Returns:
        Sum of potential ORFs found in a given chromosome or writes a
        'meta'-metafile to chromosome directory.
    """
    files = globber()
    if not args.logs:
        meta_join(files)
    else:
        print(logger(files))


if __name__ == '__main__':
    main()


