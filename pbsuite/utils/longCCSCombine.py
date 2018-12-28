#!/usr/bin/env python
import sys, argparse
from collections import defaultdict, namedtuple
from pbsuite.utils.FileHandlers import FastqFile


USAGE = "Replace single pass reads from a ZMW with CCS read when created"
class Sequence():
    def __init__(self, name, seq, qual):
        self.name = name
        self.seq  = seq
        self.qual = qual
    def toString(self):
        return "@%s\n%s\n+\n%s\n" % (self.name, self.seq, self.qual)
    
def fastqIterator(fn):
        """
        Uses yield to allow fastqFile iteration
        Input, a file name
        """
        fh = open(fn,'r')
        while True:
                
                name = fh.readline().strip()[1:]
                seq = fh.readline().strip()
                plus = fh.readline().strip()
                qul = fh.readline().strip()
                
                if qul == "":
                        break
                
                yield Sequence(name, seq, qul)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=USAGE)
    parser.add_argument("filtered_subreads", type=str, \
                        help="Fastq of single pass reads")
    parser.add_argument("ccs_reads", type=str, \
                        help="Fastq of ccs reads")
    
    parser.add_argument("-o", "--output", type=str, default=None, \
                        help="Output fastq file (STDOUT)")
    args = parser.parse_args()

    ccs = FastqFile(args.ccs_reads)
    
    cKeys = ccs.keys()
    if args.output != None:
        output = open(args.output,'w')
    else:
        output = sys.stdout
    
        
    #name: numBases
    ccsReads = defaultdict(int)
    subReads = {}
    ccsPases = defaultdict(int)
    
    #sub = FastqFile(args.filtered_subreads)
    for read in fastqIterator(args.filtered_subreads):
        ccsKey = "/".join(read.name.split('/')[:2])
        try:
            seq = ccs[ccsKey]
            ccsPases[ccsKey] += 1
            subReads[read] = len(read.seq)
            if ccsReads[ccsKey] > 0:
                continue
            ccsReads[ccsKey] = len(read.seq)
        except KeyError:
            seq = read
            
        output.write(seq.toString())
    
    output.close()

    def medianLen(lst):
        lst.sort()
        p = len(lst)/2
        return lst[p]
    
    a = len(ccsReads.keys())
    sys.stderr.write("+CCS reads    : %d\n" % (a))
    b = sum(ccsReads.values())
    sys.stderr.write("+CCS bases    : %d\n" % (b))
    c = len(subReads.keys())
    sys.stderr.write("-SUB reads    : %d\n" % (c))
    d = sum(subReads.values())
    sys.stderr.write("-SUB bases    : %d\n" % (d))
    sys.stderr.write("Avg SUB Len   : %d\n" % (int(sum(subReads.values()))/float(len(subReads))))
    sys.stderr.write("Avg CCS Len   : %d\n" % (int(sum(ccsReads.values()))/float(len(ccsReads))))
    sys.stderr.write("Avg SUB / CCS : %.3f\n" % (sum(ccsPases.values())/float(len(ccsPases))))
    
    #sys.stderr.write("Med CCS RdLen : %.2f\n" % (medianLen(ccsReads.values())))
    #sys.stderr.write("Med Sub RdLen : %.2f\n" % (medianLen(subReads.values())))
