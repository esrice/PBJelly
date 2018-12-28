#!/usr/bin/env python

import sys, glob, re, logging, json, os
from optparse import OptionParser
from collections import deque

from pbsuite.utils.setupLogging import *
from pbsuite.jelly.Support import AlignmentConnector, SUPPORTFLAGS

subreadGrab = re.compile(".*/(\d+)_(\d+)$")

def natural_sort(l): 
    """
    how unix -V sorts
    """
    l.sort()
    return l
    #Depricated, but still cool
    #convert = lambda text: int(text) if text.isdigit() else text.lower() 
    #alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    #return sorted(l, key = alphanum_key)

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


class AlignmentFileStack():
    """
    Will treat a file like a stack and help you pop of information
    that you want
    """
    def __init__(self, fn, maxEntries):
        self.fileName = fn
        self.maxEntries = maxEntries
        self.fileHandler = open(fn,'r')
        self.memory = deque()#Entries in memory
        self.stack = []
        self.wait = None
        self.finished = False
        self.__loadMem()
        
    def __loadMem(self):
        """
        Initialize the stack and take maxEntries 
        alignments
        """
        if self.finished:
            return
        
        line = None
        
        while not self.finished and len(self.memory) < self.maxEntries \
                                and line != "":
            line = self.fileHandler.readline()
            if line == "":
                self.finished = True
                self.fileHandler.close()
                return None
            self.memory.append(M4Line(line))
        
        
    def pop(self):
        """
        Just returns a single hit from the file
        """
        if self.finished:
            return None
        
        if len(self.memory) == 0:
            self.__loadMem()
            if len(self.memory) == 0: 
                #still -- this can happen in rare edge cases
                #where the memory fills up on the last entry 
                #in a file
                return None
        
        return self.memory.popleft()
        
    def gleam(self):
        """
        Return all of the reads that have query name
        """
        if len(self.stack) != 0:
            return self.stack[0].qname #Nothing to do on this stack
            
        if self.wait is None:
            self.wait = self.pop()
        
        if self.finished:
            return None
            
        name = self.wait.qname
        self.stack.append(self.wait)
        self.wait = None
        while True:
            nextGuy = self.pop()
            if nextGuy == None:
                break
            if nextGuy.qname == name:
                self.stack.append(nextGuy)
            else:
                self.wait = nextGuy
                break

        return name
                
class OverlapTableCreator():
    
    def __init__(self):
        self.parseArgs()
        self.findFiles()
        self.skimReadGroup()

    def parseArgs(self):
        parser = OptionParser()
        parser.add_option("-i", "--inputDir", default=None,\
                help="Input directory to find chunks.m4 [DEFAULT=pwd]")
        parser.add_option("-j", "--json", action="store_true",\
                help="Output table in JSON format instead of lined [DEFAULT=False]")
        parser.add_option("-o", "--output", default=None,\
                help="Output file name [DEFAULT=stdout]")
        parser.add_option("-t", "--tailMax", type="int", default=-1, \
                help=("Use PBJelly's Support module to remove discordant "
                      "alignments with greater than specified tail length "
                      "[DEFAULT=off]"))
        parser.add_option("-l", "--lengthMin", type="int", default=0, \
                help="Ignore reads (query or target) less than specified length [DEFAULT=off]")
        parser.add_option("-b", "--bestn", type="int", default=sys.maxint, \
                help=("Report only the top bestn alignment scores for a query"
                      " [DEFAULT=all]"))
        parser.add_option("-e", "--extends", action="store_true", \
                help="Only report alignments that extend query [DEFAULT=False]")
        parser.add_option("-m", "--maxEntries", type="int", default=10000, \
                help="Max number of alignments to hold in memory from each file [DEFAULT=10000]")
        parser.add_option("--debug", action="store_true",\
                help="Verbose logging")
        
        opts, args = parser.parse_args()
        
        setupLogging(opts.debug)
        
        if opts.inputDir is not None:
            self.inputDir = opts.inputDir
        else:
            self.inputDir = os.getcwd()
        if not os.path.exists(self.inputDir):
            parser.error("Input directory (%s) does not exist" % inputDir)
       
        self.outputJson = opts.json
        
        if opts.output is not None:
            self.output = open(opts.output,'w')
        else:
            self.output = sys.stdout
        
        #Filter params
        self.tailMax = opts.tailMax
        self.lengthMin = opts.lengthMin
        self.bestn = opts.bestn
        self.extends = opts.extends
        self.maxEntries = opts.maxEntries
        self.debug = opts.debug
    
    def findFiles(self):
        self.fileNames = glob.glob(os.path.join(self.inputDir, "chunk_*.chunk_*.stride_*.m4"))
        logging.info("Found %d Alignment.m4 Files" % len(self.fileNames))
        if len(self.fileNames) == 0:
            logging.error("No Files Found... Exiting.")
            exit(1)
    
    def skimReadGroup(self):
        handlers = []
        for fn in self.fileNames:
            handlers.append(AlignmentFileStack(fn, self.maxEntries))
            
        alignCon = AlignmentConnector()
        ovlcount = 0
        readCount = 0
        while True:
            availNames = []
            delete = []
            for i in handlers:
                success = i.gleam()#grabs the next read's set of hits
                if success is not None:
                    availNames.append(success)
                else:
                    delete.append(i)
            
            #get rid of the finished ones
            for d in delete:
                logging.info("Removing %d fileHandlers" % (len(delete)))
                handlers.remove(d)
            
            if len(availNames) == 0:
                logging.info("No Reads Left")
                return
            
            leastName = natural_sort(availNames)[0]
            
            readGroup = []
            for i in handlers:
                if i.stack[0].qname == leastName:
                    if i.stack[0].qseqlength < self.lengthMin:
                        i.stack = []
                        readCount += 1
                        continue
                    for read in i.stack:
                        #filtering 
                        #The not [option] or [test] structure
                        #says if this option is not on, short circuit and don't goto or.
                        if (read.tseqlength >= self.lengthMin) and \
                           (not self.extends or alignCon.extendsTarget(read) not in \
                            [SUPPORTFLAGS.none, SUPPORTFLAGS.span]) and \
                           (not self.tailMax >= 0 or not \
                            alignCon.isDiscordant(read, self.tailMax)):
                                readGroup.append(read)
                    #reset handler's stack
                    i.stack = []
            
            readCount += 1
            if len(readGroup) > 0:
                self.makeTableEntry(readGroup)
            if readCount % 10000 == 0:
                logging.info("%d reads parsed" % (readCount))
            if len(handlers) == 0:
                logging.info("Out of File Handlers")
                return
    
    def makeTableEntry(self, alignments):
        out = []
        out.append(alignments[0].qname)
        for align in alignments:
            out.extend([align.tname, align.score, align.pctsimilarity, \
                        align.qstrand, align.qstart,align.qend, \
                        align.qseqlength, align.tstrand, align.tstart, \
                        align.tend, align.tseqlength])
                        
            # align.mapqv,
            # align.clusterScore,
            # align.probScore,
            # align.numSigClusters,
        
        if self.outputJson:
            d = {out[0]:out[1:]}
            json.dump(d, self.output)
            self.output.write("\n")
        else:
            self.output.write(" ".join(map(str, out))+"\n")

def iload_json(buff, decoder=None, _w=json.decoder.WHITESPACE.match):
    """
    Internet lifter code. This is how one would load the json output
    of this script
    ...Original Documentation...
    Generate a sequence of top-level JSON values declared in the
    buffer.

    >>> list(iload_json('[1, 2] "a" { "c": 3 }'))
    [[1, 2], u'a', {u'c': 3}]
    """
    decoder = decoder or json._default_decoder
    idx = _w(buff, 0).end()
    end = len(buff)
    ret = {}
    try:
        while idx != end:
            (val, idx) = decoder.raw_decode(buff, idx=idx)
            ret.update(val)
            idx = _w(buff, idx).end()
    except ValueError as exc:
        raise ValueError('%s (%r at position %d).' % (exc, buff[idx:], idx))
    return ret
    
if __name__ == '__main__':
    OverlapTableCreator()
