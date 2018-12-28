#!/usr/bin/env python

import sys, glob, re
from optparse import OptionParser
from StringIO import StringIO

from pbsuite.utils.setupLogging import *
USAGE = """%prog <chunk> <stride>
Run this in the directory where you've made a whole bunch of alignment chunks.
For any Chunk A that has been Mapped to B, we'll make the reciprocal of B mapped to A"""

subreadGrab = re.compile(".*/(\d+)_(\d+)$")

class M4File(list):

    def __init__(self, file):
        super(list)
        if type(file) != type(StringIO()):
            file = open(file,'r')
        self.fileHandler = file
        self.__parse()
        self.fileHandler.close()
        
    def __parse(self):
        for line in self.fileHandler.readlines():
            try:
                self.append(M4Line(line.strip()))
            except TypeError, err:
                sys.stderr.write("BadM4Line! \n%s\n%s\n" % (line, str(err)) )
                sys.exit(1)
            

class M4Line():
    
    def __init__(self, line):
        data = re.split("\s+",line)
        
        self.qname          = data[0]
        self.tname          = data[1]
        self.score          = int(data[2])
        self.pctsimilarity  = float(data[3]) 
        self.qstrand        = data[4]
        self.qstart         = int(data[5])
        self.qend           = int(data[6])
        self.qseqlength     = int(data[7])
        self.tstrand        = data[8]
        tstart              = int(data[9])
        tend                = int(data[10])
        self.tseqlength     = int(data[11]) 
        #self.mapqv          = int(data[12])
        #self.clusterScore   = float(data[13])
        #self.probScore      = float(data[14])
        #self.numSigClusters = int(data[15])
        
        if self.tstrand == '1':
            self.tstart, self.tend = self.tseqlength - tend, \
                                     self.tseqlength - tstart 
        else:
            self.tstart, self.tend = (tstart, tend)
            
        #Collect subread Information
        try:
            subStart, subEnd     = subreadGrab.match(self.qname).groups()
        except AttributeError:
            subStart, subEnd = self.qstart, self.qend
        
        self.qsubstart       = int(subStart)
        self.qsubend         = int(subEnd)
        self.qsublength      = self.qsubend - self.qsubstart
        self.queryPctAligned = (self.qend - self.qstart) \
                               / float(self.qsubend - self.qsubstart)
            
        
        self.flag = 0
        self.trim = False
    
    def __str__(self, convert=False):
        #if convert:#Put back into original space
        if self.tstrand == '1':
            tstart, tend = self.tseqlength - self.tend, \
                            self.tseqlength - self.tstart 
        else:
            tstart = self.tstart
            tend = self.tend

        return " ".join(map(str, [self.qname, self.tname, self.score, \
                                  self.pctsimilarity, self.qstrand,   \
                                  self.qstart, self.qend, self.qseqlength, \
                                  self.tstrand, tstart, tend, \
                                  self.tseqlength]))
    
    def toBed(self):
        """
        Returns the M4Line to a bed file
        """
        if self.tstrand == '1':
            strand = "-"
            chromStart = str(self.tstart - (self.qseqlength - self.qend))
            chromEnd = str(self.tend + self.qstart)
        else:
            strand = "+"
            chromStart = str(self.tstart - self.qstart)
            chromEnd = str(self.tend + (self.qseqlength - self.qend))
        #chromStart = str(self.tstart - self.qstart)
        #chromEnd = str(self.tend + (self.qseqlength - self.qend))
        chrom = self.tname
        name = self.qname
        score = str(self.score)
        thickStart = str(self.tstart)
        thickEnd = str(self.tend)
        itemRgb = str(int(self.pctsimilarity))
        return "\t".join([chrom, chromStart, chromEnd, name, score, strand, \
                          thickStart, thickEnd, itemRgb])


def makeReciprocal(chunks, stride):
    for i in range(chunks):
        for s in range(stride):
            for j in range(i,chunks):
                if i != j: #A vs A
                    #Make reciprocal
                    alignments = M4File("chunk_%d.chunk_%d.stride_%d.m4" % (i,j,s))
                    fout = open("chunk_%d.chunk_%d.stride_%d.m4" % (j, i, s),'w')
                    for align in alignments:
                        fout.write(str(switchTargetQuery(align))+"\n")
                    fout.close()
    

def switchTargetQuery(align):
    """
    Turns the query into the target and the target into the query
    """
    if align.tstrand == '1':
        ret = M4Line(" ".join(map(str,[align.tname, align.qname, align.score, align.pctsimilarity, '0', \
                     align.tseqlength - align.tend, align.tseqlength - align.tstart, align.tseqlength, \
                     '1', align.qseqlength - align.qend, align.qseqlength - align.qstart, align.qseqlength])))
    else:
        ret = M4Line(" ".join(map(str,[align.tname, align.qname, align.score, align.pctsimilarity, '0', \
                     align.tstart, align.tend, align.tseqlength, \
                     '0', align.qstart, align.qend, align.qseqlength])))
    return ret


if __name__ == '__main__':
    parser = OptionParser(USAGE)
    opts, args = parser.parse_args()
    try:
        chunk, stride = args
        chunk = int(chunk)
        stride = int(stride)
    except:
        parser.error("Couldn't Parse Arguments")
    
    makeReciprocal(chunk, stride)
    
