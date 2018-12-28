#!/usr/bin/env python
import pysam, sys, argparse
from pbsuite.utils.FileHandlers import revComp

USAGE="Use pysam to extract fastq sequences from a bam"

parser = argparse.ArgumentParser(description=USAGE, \
        formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("BAM", help="input bam file")
parser.add_argument("POS", default=None, nargs="?", \
                    help="`samtools view` format of the location to grab reads from (optional)")
parser.add_argument("-f", "--flag", default=None,
                    help="Only print reads with flags (comma sep)")
parser.add_argument("-F", "--Flag", default=None,
                    help="Only print reads WITHOUT flags (comma sep)")
args = parser.parse_args()
args.flag = [int(x) for x in args.flag.split(',')] if args.flag is not None else []
args.Flag = [int(x) for x in args.Flag.split(',')] if args.Flag is not None else []

s = pysam.Samfile(args.BAM)

if args.POS is not None:
    chrom, extra = args.POS.split(':')
    start, end = [int(x) for x in extra.split('-')]
    get = s.fetch(reference=chrom, start=start, end=end)
else:
    get = s

for read in get:
    if [read.flag & int(x) for x in args.flag].count(False) \
       or  [read.flag & int(x) for x in args.Flag].count(True):
        continue
        
    if read.is_reverse:
        sys.stdout.write("@{0}\n{1}\n+\n{2}\n".format(read.qname, read.seq.translate(revComp)[::-1], read.qual[::-1]))
    else:
        sys.stdout.write("@{0}\n{1}\n+\n{2}\n".format(read.qname, read.seq, read.qual))
