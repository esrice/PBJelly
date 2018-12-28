#!/usr/bin/env python
import sys
from optparse import OptionParser
from FileHandlers import FastaFile

USAGE = """USAGE: %prog <input.fasta> <output.qual> [--options]
Creates fake quality information for fasta files that do not have associated \
quality information"""

if __name__ == '__main__':
    parser = OptionParser(USAGE)
    parser.add_option("-s", "--score", default="40", \
                      help=("Score to give every base in fasta file.\n" \
                            "DEFAULT = 40"))
    
    opts, args = parser.parse_args()
    
    if len(args) != 2:
        parser.error("Error! Expected exactly 2 arguments.")
    fastaName, outName = args
    fasta  = FastaFile(fastaName)
    fout = open(outName,'w')
    
    for entry in fasta:
        fout.write(">" + entry + "\n" + \
                  ((opts.score + " ") * len(fasta[entry])) + "\n")
    
    fout.close()
