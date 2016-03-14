#!/usr/bin/env python
from __future__ import print_function
import sys


def modulo_check(NUM):
    """gives greatest whole divisor from 1 to 50

    Args:
        NUM (int): value to be evaluated
    Returns:
        greatest divisor between 1 and 50
    """
    mod_check = []
    for i in range(1, 50):
        if NUM % i is 0:
            mod_check.append(i)
        else:
            pass
    return max(mod_check)


def main():
    num = int(sys.argv[1])
    div = modulo_check(num)
    out = num/div
    print("Greatest Divisor: %d" %(div))
    print("Lines per file: %d" %(out))


if __name__ == '__main__':
    main()
