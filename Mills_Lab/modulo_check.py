#!/usr/bin/env python
from __future__ import print_function
import sys
import os
import re


def modulo_check(NUM):
    """gives greatest whole divisor from 1 to 50

    Args:
        NUM (int): value to be evaluated
    Returns:
        greatest divisor between 1 and 50
    """
    mod_check = []
    for i in range(1, 53):
        if NUM % i is 0:
            mod_check.append(i)
        else:
            pass
    return max(mod_check)


def main():
    pat1 = re.compile(r'unannotatedchr\d{1,2}')
    pat2 = re.compile(r'/home.*/chr\d{1,2}/')
    infile = sys.argv[1]
    output = pat1.search(infile).group()
    path = pat2.search(infile).group()
    line_count = os.popen("gzip -cd %s | wc -l" % infile)
    num = int((line_count.next()).rstrip())
    div = modulo_check(num)
    out = num/div
    os.popen("gzip -cd %s | split -l %d  - %s" % (infile, out, output))
    print("Greatest Divisor: %d" % div)
    print("Lines per file: %d" % out)
    os.chdir(path)
    os.popen("find . -type f -exec mv '{}' '{}'.vcf \;")
    os.popen("rm -f *.gz.vcf")
    os.popen('for f in *.vcf; do   bgzip "$f"; done')

if __name__ == '__main__':
    main()
