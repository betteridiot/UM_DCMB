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
    pass
else:
    parent = args.logs
    string = "*.logs"
if not args.meta:
    pass
else:
    parent = args.meta
    string = "metadata"


def globber():
    """Walks the directory and globs all of the *.log files

    Returns:
        A list of *.log files for a given chromosome.
    """
    cwd = parent
    os.chdir(cwd)
    globs = []
    term = string
    for dir, _, _ in os.walk(os.getcwd()):
        globs.extend(glob.glob(os.path.join(dir, term)))
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
            reader = (int(string) for string in reader.split() if string.isdigit())
        num_orfs.extend(reader)
    return str(sum(num_orfs))


def meta_join(lst):
    """Writes a 'meta'-metafile to chromosome directory taken from the different
    ORFSNuPer runs.

    Args:
        lst (list): list of metadata files for a given chromosome
    """
    first = True
    temp = TemporaryFile("a+b")
    with open(temp, "ab+") as outfile:
        for entry in lst:
            with(entry, "rb") as infile:
                if first:
                    header = infile.readline().strip().split()
                    for line in infile:
                        outfile.write(line)
                    first = False
                else:
                    next(infile)
                    for line in infile:
                        outfile.write(line)
        outfile.seek(0)
        table = [line.strip().split() for line in outfile]
        table.sort(key=itemgetter(0))
        with open(parent + "meta_full", 'wb') as final:
            writer = csv.writer(final, delimiter='\t')
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
        print(logger(files))
    else:
        meta_join(files)


if __name__ == '__main__':
    main()
