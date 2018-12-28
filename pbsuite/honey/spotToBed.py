#!/usr/bin/env python
import sys
import argparse
from pbsuite.utils import VCFIO

USAGE = "Turn .spots results into a .bed"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=USAGE, \
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", metavar="SPOTS", default='-', nargs="?", \
                        help="Results to convert (stdin)")
    parser.add_argument("-b", "--brief", action="store_true",\
                        help="Only output columns 1-4 of .bed")
    
    args = parser.parse_args()
    if args.input == '-':
        fh = sys.stdin
    else:
        fh = open(args.input,'r')
       
    for line in fh.readlines():
        if line.startswith("##"):
            continue
        elif line.startswith("#"):
            h = line.strip()[1:]
            header = {}
            for pos, item in enumerate(h.split('\t')):
                header[item] = pos
            continue
        chrom, s, e, type, size, info = line.strip().split('\t')
        name = "%s.%s" % (type, size)
        print "\t".join([chrom, s, e, name, type, size])
    fh.close()
            
            

