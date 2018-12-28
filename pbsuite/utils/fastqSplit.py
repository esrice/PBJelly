#!/usr/bin/env python

import sys
from optparse import OptionParser
from collections import namedtuple
from FileHandlers import wrap, qwrap, FastqFile
from StringIO import StringIO

USAGE = """Usage: %prog <input.fastq> <baseName>
Splits a fastq into <baseName>.fasta and <baseName>.qual
Assumes Sanger Encoded Phred Scores in fastq
"""

def __parseArgs():
    parser = OptionParser(usage=USAGE)
    
    opts, args = parser.parse_args(sys.argv)
    if len(args) != 3: parser.error('Expected 2 arguments')
    
    return args[1:]

def fastqIter( fn ):
    fh = open(fn, 'r')
    FastQEntry = namedtuple("FastQEntry", "name seq qual")
    while True:
        name = fh.readline().strip()[1:]
        if name == "": break
        #seq grab
        line = fh.readline().strip()
        seq = StringIO()
        
        while not line.startswith('+'):#Assuming no name...
            seq.write(line)
            line = fh.readline().strip()
        seq = seq.getvalue()
        seqLen = len(seq)

        qual = ""
        curLen = 0

        while curLen != len(seq):
            line = fh.readline().strip()
            if line == "":
                sys.stderr.write("Bad Fastq File: Last attempted entry = %s\n" % (name))
                exit(10)
            curLen += len(line)
            qual += line
        

        yield FastQEntry(name, seq, qual)

def phredToQual( qual ):
    """
    Take a qual string that is phred/sanger encoded
    turn it into a list of quals
    """
    return map(lambda x: ord(x)-33, list(qual))
    
if __name__ == '__main__':
    fastq, baseName = __parseArgs()
    
    fout = open(baseName+".fasta", 'w')
    qout = open(baseName+".qual", 'w')
    fastq = FastqFile(fastq)
    for name in fastq:
        entry = fastq[name]
        fout.write(">%s\n%s\n" % (entry.name, wrap(entry.seq)))
        qout.write(">%s\n%s\n" % (entry.name, qwrap(phredToQual(entry.qual))))
    
    fout.close()
    qout.close()
