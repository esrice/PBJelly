#!/usr/bin/env python
import sys, argparse
from pbsuite.utils.setupLogging import *

USAGE = "Split input.[fastq|fasta] file into N subfiles"

parser = argparse.ArgumentParser(description=USAGE, \
        formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("input", metavar="INPUT", \
                    help="File to be split")
parser.add_argument("nSplits", metavar="N", type=int, \
                    help="Number of splits to create")
parser.add_argument("-o", "--outPrefix", default=None,\
                    help="Prefix for output files")
#parser.add_argument("-m", "--multiline", action="store_true",\
                    #help="Input file is expected to be multiline per entry")

args = parser.parse_args()
format = 0
if args.input.lower().endswith(".fastq"):
    format = 4
    suffix = ".fastq"
elif args.input.lower().endswith(".fasta"):
    format = 2
    suffix = ".fasta"
if not format:
    logging.error("Input must be .fastq or .fasta")
    exit(0)

if args.outPrefix is None:
    args.outPrefix = args.input[:len(".fasta")] 

outFiles = []
for i in range(args.nSplits):
    outFiles.append(open(args.outPrefix + str(i) + suffix, 'w'))

fh = open(args.input, 'r')

index = 0

while True:
    name = fh.readline()
    seq = fh.readline()
    
    if name == "":
        break
    outFiles[index].write(name + seq)
    if format == 4:
        #plus
        outFiles[index].write(fh.readline())
        #qual
        outFiles[index].write(fh.readline())

    index += 1#I could mod here..
    if index >= args.nSplits:
        index = 0

for f in outFiles:
    f.close()
