#!/usr/bin/env python

import argparse, json
from pbsuite.jelly.Jelly import JellyProtocol
from pbsuite.utils.FileHandlers import FastaFile, FastqFile
from pbsuite.utils.summarizeAssembly import getStats

USAGE = """Get statistics on fasta/fastq sequences recorded in a Protocol.xml"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=USAGE, \
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("xml", metavar="XML", type=str, \
                        help="Protocol.xml with inputs listed")
    args = parser.parse_args()
    protocol = JellyProtocol(args.xml)
    seqLengths = []
    for i in protocol.inputs:
        if i.endswith(".fasta"):
            f = FastaFile(i)
            for j in f.values():
                seqLengths.append(len(j))
        if i.endswith(".fastq"):
            f = FastqFile(i)
            for j in f.values():
                seqLengths.append(len(j.seq))
    print "Read Stats", json.dumps(getStats(seqLengths), indent=4)
