#!/usr/bin/env python
from __future__ import print_function
import glob
import os
import sys


def globber():
    """Walks the directory and globs all of the *.log files

    Returns:
        A list of *.log files for a given chromosome.
    """
    cwd = sys.argv[1]
    os.chdir(cwd)
    globs = []
    for dir, _, _ in os.walk(os.getcwd()):
        globs.extend(glob.glob(os.path.join(dir, "*.log")))
    return globs


def main():
    """Takes the glob and iteratively pulls out number of potential ORFs

    Returns:
        Sum of potential ORFs found in a given chromosome.
    """
    files = globber()
    num_orfs = []
    for entry in files:
        with open(entry, 'rb') as f:
            reader = f.readline()
            reader = (int(string) for string in reader.split() if string.isdigit())
        num_orfs.extend(reader)
    print(str(sum(num_orfs)))


if __name__ == '__main__':
    main()
#