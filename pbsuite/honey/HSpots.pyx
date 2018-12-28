#!/usr/bin/env python

import sys, time, re, gc, argparse, logging, bisect, math, copy, shutil, multiprocessing
import pysam, numpy, h5py
from collections import defaultdict, Counter
from tempfile import NamedTemporaryFile

from pbsuite.utils.setupLogging import *
from pbsuite.utils.CommandRunner import exe
from pbsuite.utils.FileHandlers import M5File, revComp
from pbsuite.banana.Polish import consensus
from pbsuite.jelly.Assembly import blasr

VERSION = "14.10.14"
AUTHOR = "Adam English"

USAGE = """\
Detect 'smaller' SVs via measuring discordance between sample and reference in long reads.
"""

#########################
## --- Global Vars --- ##
#########################

#Biggest integer I want to deal with
BIGINT  = 2000
BIGINTY = numpy.float32
# NUMPY ARRAY HDF5 COLUMNS AND SIZES
COLUMNS = ["coverage", "matches", "insertions", "deletions"]
#Index of columns
COV  = 0
MAT  = 1
INS  = 2  
DEL  = 3  
#Must not exceed 300,000 data points
CHUNKSHAPE = (4, 70000)

##############################
## --- Helper Functions --- ##
##############################

#{{{ http://code.activestate.com/recipes/511478/ (r1)
import math
import functools

def percentile(N, percent):
    """
    Find the percentile of a list of values.

    @parameter N - is a list of values. Note N MUST BE already sorted.
    @parameter percent - a float value from 0.0 to 1.0.

    @return - the percentile of the values
    """
    if not N:
        return None
    k = (len(N)-1) * percent
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return N[int(k)]
    d0 = N[int(f)] * (c-k)
    d1 = N[int(c)] * (k-f)
    return d0+d1
# median is 50th percentile.
## end of http://code.activestate.com/recipes/511478/ }}}

def expandCigar(cigar):
    """
    Turns the abbreviated cigar into the full array
    0 = M; 1 = I; 2 = D
    -- C translate candidate -- 
    """
    ret = []
    for t,s in cigar:
        if t < 3: #remove non mid (dangerous if blasr changes)
            ret.extend([t]*s)
    return ret

class SpotResult():
    def __init__(self, chrom=None, out_start=None, start=None, in_start=None, \
                                in_end=None, end=None, out_end=None, tags=None):
        self.chrom = chrom
        self.out_start = out_start
        self.start = start
        self.in_start = in_start
        self.in_end = out_end
        self.end = end
        self.out_end = out_end
        self.size = -1
        self.tags = tags if tags is not None else {}
        
    def offset(self, start):
        """
        moves the spot to an offset
        """
        if self.out_start is not None:
            self.out_start += start
        if self.start is not None:
            self.start += start
        if self.in_start is not None:
            self.in_start += start
        if self.in_end is not None:
            self.in_end += start
        if self.end is not None:
            self.end += start
        if self.out_end is not None:
            self.out_end += start
    
    def fetchbounds(self):
        """
        return (start,end) tuple of spot's boundaries
        """
        pnts = [x for x in [self.out_start, self.start, self.in_start, \
                            self.in_end, self.end, self.out_end] \
                            if x is not None]
        return min(pnts), max(pnts)
    
    def estimateSize(self):
        """
        estimate the sv's size by using the mean or bounds
        """
        if 'szMean' in self.tags:
            self.size = self.tags["szMean"]
        else:
            s,e = self.fetchbounds()
            self.size = e-s
    
    def qregstr(sefl):
        """
        returns quick region string chrom:start-end
        """
        return "%s:%d-%d" % (self.chrom, self.start, self.end)
    
    def __str__(self):
        """
        changes a spot named tuple to a svp string
        """
        tag = []
        for key in self.tags:
            if key == 'label':
                self.type = self.tags[key]
            else:
                try:
                    tag.append("%s=%0.3f" % (str(key), self.tags[key]))
                except TypeError:
                    tag.append("%s=%s" % (str(key), str(self.tags[key])))
        
        
        tag = ";".join(tag)
        dat = [self.chrom, self.out_start, self.start, self.in_start, \
               self.in_end, self.end, self.out_end, self.type, self.size, \
               tag]

        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}".format(*dat) \
                .replace("None", ".")
 
def countErrors(reads, offset, size, args, readCount=None):
    """
    Sum the errors over any particular reference base
    """
    container = numpy.zeros( ( len(COLUMNS), size ), dtype=BIGINTY )
    numReads = 0.0 #number finished
    if readCount is not None:
        logging.info("%d reads to parse" % readCount)
    
    for align in reads:
        cigar = expandCigar(align.cigar)
    
        #get starts within region
        regionStart = 0 if align.pos < offset else align.pos - offset
        regionEnd = size if align.aend > (offset + size) else align.aend - offset
        
        #I'm always starting at the beginning of the read,
        # but the beginning may hit before my regionStart
        start = align.pos - offset
        
        #check covering bases
        container[COV,  regionStart : regionEnd] += BIGINTY(1)
        #MAQ?
        
        #previous base was an insert prevent multiple base ins
        pins = False
        pinsz = 0
        pdel = False
        pdels = None
        
        def pinsLoad(start, size):
            if not pins:
                return False, 0
            if size >= args.minIndelErr:
                begin = max(0,start-(size/2))
                end = min(start+(size/2), container.shape[1])
                container[INS, begin:end] += BIGINTY(1)
            return False, 0
        
        def pdelLoad(startPos, curPos):
            if not pdel:
                return False, None
            if curPos-startPos >= args.minIndelErr:
                container[DEL, startPos:curPos] += BIGINTY(1)
            return False, None
        
        for code in cigar:
            if start < regionStart or start >= regionEnd:
                if code != 1: 
                    start += 1
                continue
            elif code == 0:
                container[MAT, start] += BIGINTY(1)
                start += 1
                pins, pinsz = pinsLoad(start, pinsz)
                pdel, pdels = pdelLoad(pdels, start)
            elif code == 1: #ins
                pins = True
                pinsz += 1
                pdel, pdels = pdelLoad(pdels, start)
            elif code == 2: #del
                if pdels is None:
                    pdel = True
                    pdels = start
                start += 1
                pins, pinsz = pinsLoad(start, pinsz)
    
    return container

def preprocessSignal(signal, coverage):
    """
    Normalize and print stats returning data and it's std
    """
    rate = numpy.convolve(signal/coverage, avgWindow, "same")
    rate[numpy.any([numpy.isinf(rate), numpy.isnan(rate)], axis=0)] = 0
    mu = numpy.mean(rate)
    sd = numpy.std(rate)
    logging.info("RateMean %f  -- RateStd  %f" % (mu, sd))
    return rate, mu, sd
    
def callHotSpots(data, offset, args): #threshPct, covThresh, binsize, offset):
    """
    """
    ret = []
    avgWindow = numpy.ones(binsize, dtype=numpy.float16)/float(args.binsize)
    #coverage
    cov = numpy.convolve(data[COV], avgWindow, "same")
    covTruth = numpy.all([cov >= args.minCoverage, cov <= args.maxCoverage], axis=0)
    logging.info("MaxCov:%d MeanCov:%d StdCov:%d MinCov:%d" \
            % (numpy.max(data[COV]), numpy.mean(data[COV]), \
               numpy.std(data[COV]), numpy.min(data[COV])))
    
    #ins
    logging.info("INS processing")
    ins, mu, sd = preprocessSignal(data[INS], data[COV])
    startPoints = makeSpots(numpy.all([ins*cov >= args.threshold, covTruth], axis=0), args.binsize, 's')
    for start, end, s in startPoints:
        mySpot = SpotResult(tags={"label":"INS"})
        mySpot.start = start
        mySpot.end = end
        ret.append(mySpot)
    del(ins)
    
    #dele
    logging.info("DEL processing")
    dele, mu, sd = preprocessSignal(data[DEL], data[COV])
    startPoints = makeSpots(numpy.all([dele*cov >= args.threshold, covTruth], axis=0), args.binsize, 's')
    for start, end, s in startPoints:
        mySpot = SpotResult(tags={"label":"DEL"})
        mySpot.start = start
        mySpot.end = end
        ret.append(mySpot)
    del(dele)
    
    return ret
    
def makeSpots(truth, binsize, label):
    """
    make the points for the truth set made from the data container
    truth = numpy.array() with boolean values
    """
    #prevent weirdness
    truth[-1] = False
    shift = numpy.roll(truth, 1)
    
    starts = truth & ~shift
    ends = ~truth & shift
    
    points = zip(numpy.nonzero(starts)[0], numpy.nonzero(ends)[0])
    npoints = []
    if len(points) == 0:
        return npoints
    curStart, curEnd = points[0]
    
    #compress: <-- Don't need anymore...?
    for start, end in points[1:]:
        #if start - curEnd <= binsize:
            #curEnd = end
        #else:
        npoints.append((curStart, curEnd, label))
        curStart = start
        curEnd = end
    
    npoints.append((curStart, curEnd, label))
    
    return npoints

def supportingReadsFilter(spot, args):
    """
    filters insertions or deletion spots based on errors
    """
    if spot.tags["label"] == "INS":
        errId = 1
        errLab = 'insertion'
    elif spot.tags["label"] == "DEL":
        errId = 2
        errLab = 'deletion'
    else:#don't worry about other types
        return False

    begin, ending = spot.fetchbounds()
    begin -= args.buffer #abs(begin-ending)*.5
    ending += args.buffer #abs(begin-ending)*.5
    #do the hard work
    reads = args.bam.fetch(str(spot.chrom), begin, ending)
    totSizes = []
    coverage = 0
    nReadsErr = 0
    #For tandem
    strandCnt = {True: 0, False: 0}
    
    #count reads and errSizes
    for i in reads:
        mySize = 0
        coverage += 1
        start = i.pos - 1
        cigar = expandCigar(i.cigar)
        curSize = 0
        extraSize = 0
        readHasErr = False
        
        #What if I just intersect any stretches of errors with my boundaries.
        #Then for insertions I'll keep coordinates
        #For deletions I'll user outer bounds?
        for code in cigar: 
            if code != 1:
                start += 1
            #must be in region
            if start < begin:
                continue
            if start >= ending:
                break
            
            if code == errId:
                curSize += 1
            if curSize != 0 and code != errId:
                if curSize >= args.minIndelErr:
                    readHasErr = True
                    mySize += curSize
                elif curSize > 1:#1bp errors will inflate
                    extraSize += curSize
                curSize = 0
            

        if readHasErr and mySize >= args.minIndelSize:
            nReadsErr += 1
            totSizes.append(mySize + extraSize)
            strandCnt[i.is_reverse] += 1
    
    spot.tags["strandCnt"] = "%d,%d" % (strandCnt[False], strandCnt[True])
    if len(totSizes) == 0:
        logging.debug("no %s found!? %s" % (errLab, str(spot)))
        return True # true you should filter
    
    if len(totSizes) < max(math.ceil(coverage * args.minIndelPct), args.minErrReads):
        logging.debug("not large cnt %s found %s " % (errLab, str(spot)))
        return True
    
    totSizes.sort()
    totSizes = numpy.array(totSizes)
    mean = totSizes.mean()
    median = numpy.percentile(totSizes, 50)
    firstQ = numpy.percentile(totSizes, 25)
    thirdQ = numpy.percentile(totSizes, 75)
    
    logging.debug("PassFilt %s" % (str(spot)))   
    logging.debug("cov    %d" % coverage )
    logging.debug("size %d %s" % (len(totSizes), str(totSizes)))
    logging.debug("mean   %d" % mean )
    logging.debug("median %d" % median)
    logging.debug("firstQ %d" % firstQ)
    logging.debug("thirdQ %d" % thirdQ)
    
    spot.tags["szCount"]  = int(nReadsErr)
    spot.tags["szMean"]   = int(mean)
    spot.tags["szMedian"] = int(median)
    spot.tags["sz1stQ"]   = int(firstQ)
    spot.tags["sz3rdQ"]   = int(thirdQ)
    return False

def consensusCalling(spot, args):
    """
    Make a consensus of all the reads in the region and identify all of the SVs in the region
    """
    def readTrim(read, start, end):
        """
        Trims a pysam.AlignedRead to only include the sequence that's aligned (or should be aligned)
        between start and end on reference
        returns the sequence and quality
        """
        score = 0
        if not read.is_unmapped:
            regTrim = 0
            upS = read.cigar[0][1] if read.cigar[0][0] == 4 else 0
            dnS = read.cigar[-1][1] if read.cigar[-1][0] == 4 else 0
            
            trimS = None
            trimE = None
            if start > read.pos:
                for queryPos, targetPos in read.aligned_pairs:
                    if trimS is None and targetPos >= start:
                        trimS = queryPos
            else:
                score += abs(read.pos - start)
            if end < read.aend:
                for queryPos, targetPos in read.aligned_pairs[::-1]:
                    if trimE is None and targetPos <= end:
                        trimE = queryPos
            else:
                score += abs(read.aend-end)
            
            if trimS is not None:
                trimS = max(0, trimS) + upS
            else:
                trimS = 0
                    
            if trimE is not None:
                trimE = min(len(read.seq), trimE)  - dnS
            else:
                trimE = len(read.seq)
            seq = read.seq[trimS:trimE]
            qual = read.qual[trimS:trimE]
            if not read.is_reverse:
                seq = seq.translate(revComp)[::-1]
                qual = qual[::-1]
        
        return seq, qual
    
    #END readTrim
    
    chrom, start, end = spot.chrom, spot.start, spot.end
    buffer = args.buffer
    bam = args.bam
    #work
    supportReads = []
    spanReads = []
    #Fetch reads and trim
    totCnt = 0
    for read in bam.fetch(chrom, start-buffer, end+buffer):
        seq, qual = readTrim(read, start-buffer, end+buffer)
        if read.pos < start-300 and read.aend > end+300:
            spanReads.append((len(seq), seq, qual))
        else:
            supportReads.append((seq, qual))
        totCnt += 1
        
    if len(spanReads) == 0:
        logging.info("noone spans - consensus aborted. %s" % (str(spot)))
        spot.tags["noSpan"] = True
        return [spot]
        
    spanReads.sort(reverse=True)
    refread = spanReads[0]
    logging.debug("%d reads %d support" % (totCnt, len(supportReads)))
    supportReads.extend([(x[1], x[2]) for x in spanReads[1:]])
    #read that spans most of the region goes first
    #use the rest for cleaning
    
    #building consensus sequence
    foutreads = NamedTemporaryFile(suffix=".fastq")
    for id, i in enumerate(supportReads):
        foutreads.write("@%d\n%s\n+\n%s\n" % (id, i[0], i[1]))
    foutreads.flush()
    
    foutref = NamedTemporaryFile(suffix=".fasta")
    foutref.write(">%s:%d-%d\n%s"%("ecoli", start, end, refread[1]))
    foutref.flush()
    
    alignOut = NamedTemporaryFile(suffix=".m5")
    
    blasr(foutreads.name, foutref.name, bestn=1, nproc=1, outname=alignOut.name)
    #shutil.copyfile(foutreads.name, "sup.fastq")
    #shutil.copyfile(foutref.name, "base.fasta")
    #shutil.copyfile(alignOut.name, "align.m5")
    if not args.pbdagcon:
        aligns = M5File(alignOut.name)
        con = ">con\n%s\n" % consensus(aligns).sequence
    else:
        logging.debug("pbdagcon")
        r, con, e = exe("pbdagcon -m 25 -c 1 -t 0 %s" % (alignOut.name))
        logging.debug(str(r) + " - " + str(e))
        con = con[con.index("\n")+1:]
        logging.debug("MySeq: " + con)
        #Check if con is blank
    
    conOut = NamedTemporaryFile(suffix=".fasta")
    conOut.write(con)
    conOut.flush()
    refOut = NamedTemporaryFile(suffix=".fasta")
    refOut.write(">%s:%d-%d\n%s\n" % (chrom, start, end, \
                 args.reference.fetch(chrom, start-buffer, end+buffer)))
    refOut.flush()
    
    #map consensus to refregion
    varSam = NamedTemporaryFile(suffix=".sam")
    cmd = "blasr %s %s -sam -bestn 1 -affineAlign -out %s" % (conOut.name, refOut.name, varSam.name)
    logging.debug(cmd)
    logging.debug(exe(cmd))
    
    foutreads.close()
    foutref.close()
    alignOut.close()

    #convert sam to bam
    input = pysam.Samfile(varSam.name)
    varBam = NamedTemporaryFile(suffix=".bam")
    output = pysam.Samfile(varBam.name, 'wb', template=input)
    nReads = 0
    for read in input:
        output.write(read)
        nReads += 1
    logging.info("%d consensus reads created" % (nReads))
    varSam.close()
    input.close()
    output.close()
    
    #do pileup for sequence
    pysam.sort(varBam.name, varBam.name[:-4])
    pysam.index(varBam.name)
    bam = pysam.Samfile(varBam.name, 'rb')
    
    mySpots = []
    for pos in bam.pileup():
        size = pos.pileups[0].indel
        if abs(size) < args.minIndelSize or size == 0:
            continue
        newspot = copy.deepcopy(spot)
        if size > 0:
            newspot.start = pos.pos + start - buffer
            newspot.end = pos.pos + start - buffer
            align = pos.pileups[0]
            newspot.tags["seq"] = align.alignment.seq[align.qpos : align.qpos + align.indel]
            newspot.size = size
            newspot.tags["label"] = "INS"
            mySpots.append(newspot)
        elif size < 0:
            newspot.start = pos.pos + start - buffer
            newspot.end = pos.pos + abs(size) + start - buffer
            #newspot.tags["seq"] = args.reference.fetch(chrom, pos.pos, pos.pos + abs(size))
            newspot.size = -size
            newspot.tags["label"] = "DEL"
            mySpots.append(newspot)
    bam.close()
    varBam.close()
    logging.debug("%d spots found" % (len(mySpots)))
    return mySpots
    
def parseArgs(argv, established=False):
    parser = argparse.ArgumentParser(prog="Honey.py spots", description=USAGE, \
            formatter_class=argparse.RawDescriptionHelpFormatter)
    
    ioGroup = parser.add_argument_group("I/O Argument")
    ioGroup.add_argument("bam", metavar="BAM", type=pysam.Samfile, \
                        help="BAM containing mapped reads")
    if established:
        ioGroup.add_argument("hon", metavar="HON.H5", type=str, \
                        help="HON.h5 containing signal data")
        
    ioGroup.add_argument("-r", "--region", type=str, default=None,\
                        help="Only call spots in region.bed")
    #ioGroup.add_argument("--chrom", type=str, default=None, \
                        #help="Only call spots on one chromosome (%(default)s)")
    ioGroup.add_argument("-o", "--output", type=str, default=None, \
                        help="Basename for output (BAM.hon)")
    #ioGroup.add_argument("-n", "--nproc", type=int, default=1, \
                        #help="Number of processors to use (only 4 consensus) (%(default)s)")
    
    pGroup = parser.add_argument_group("Spot-Calling Threshold/Filtering Arguments")
    pGroup.add_argument("-s", "--noCallSpots", action="store_true",\
                        help=("Don't call spots where error rates spike (False)"))
    pGroup.add_argument("-b", "--binsize", type=int, default=100, \
                        help="binsize for window averaging (%(default)s)")
    pGroup.add_argument("-e", "--threshold",  type=float, default=2,
                        help="Minimum Spot Threshold (%(default)s)")
    pGroup.add_argument("-c", "--minCoverage", type=int, default=2, \
                        help="Minimum coverage of a region (%(default)s)")
    pGroup.add_argument("-C", "--maxCoverage", type=int, default=BIGINT, \
                        help="Maximum coverage of a region (%(default)s)")
    pGroup.add_argument("-m", "--minIndelErr", type=int, default=5,
                        help="Minimum size of an indel error to be counted (%(default)s)")
    pGroup.add_argument("-i", "--minIndelSize", type=int, default=50, \
                        help="Minimum indel SV size (%(default)s)")
    pGroup.add_argument("-I", "--minIndelPct", type=float, default=0.20, \
                        help="Minimum pct of reads with indel (max(%(default)s*cov, minErrReads)")
    pGroup.add_argument("-E", "--minErrReads", type=int, default=2, \
                        help="Minimum number of reads with indel (max(minIndelPct, %(default)s))")
    pGroup.add_argument("--spanMax", type=int, default=2000, \
                        help="Maximum Size of spot to be called (%(default)s)")
    
    aGroup = parser.add_argument_group("Consensus Arguments")
    aGroup.add_argument("--noConsensus", action="store_true", \
                        help="Turn off consensus calling, just report spots (False)")
    aGroup.add_argument("--buffer", default=1000, \
                        help="Buffer around SV to assemble (%(default)s)")
    aGroup.add_argument("--reference", default=None, type=pysam.Fastafile, \
                        help="Sample reference. Required with consensus calling (None)")
    aGroup.add_argument("--pbdagcon", action="store_true", \
                        help="Use pbdagcon for consensus.")
    aGroup.add_argument("--blasr", default="blasr", \
                        help="Path to blasr if it's not in the env")
    parser.add_argument("--debug", action="store_true", \
                        help="Verbose logging")
    
    args = parser.parse_args(argv)
    if args.maxCoverage > BIGINT:
        logging.error("Max Coverge must be less than %d" % (BIGINT))
        exit(0)
    if args.output is None:
        args.output = args.bam.filename[:-4]+".hon"
    
    if not args.noConsensus:
        if args.reference is None:
            logging.error("Reference is required with consensus calling")
            exit(0)
        #Check pbdagcon
    return args

def run(argv):
    numpy.seterr(all="ignore")
    args = parseArgs(argv)
    setupLogging(args.debug)
    
    try:
        if args.bam.header["HD"]["SO"] != "coordinate":
            logging.warning("BAM is not sorted by coordinates! Performance may be slower")
    except KeyError:
        logging.warning("Assuming BAM is sorted by coordinate. Be sure this is correct")
   
    if not args.noCallSpots:
        hotspots = open(args.output+".spots", 'w')
        hotspots.write("#CHROM\tOUTERSTART\tSTART\tINNERSTART\tINNEREND\tEND\tOUTEREND\tTYPE\tSIZE\tINFO\n")
    
    regions = []
    if args.region:
        fh = open(args.region,'r')
        for line in fh.readlines():
            data = line.strip().split('\t')
            regions.append((data[3], data[0], int(data[1]), int(data[2])))
        fh.close()
    else:
        for chrom, size in zip(args.bam.references, args.bam.lengths):
            regions.append((chrom, chrom, 0, size))
    
    h5Name = args.output + ".h5"
    results = h5py.File(h5Name, 'w')
    results.attrs["version"] = VERSION
    results.attrs["columns"] = COLUMNS
    results.attrs["parameters"] = str(args)
    results.close()

    totReads = 0
    totSpots = 0
    for groupName, chrom, start, end in regions:
        size = end - start
        regName =  "%s:%d-%d" % (chrom, start, end)
        logging.info("Making container for %s (%s %d bp)" % (groupName, regName, size))
        results = h5py.File(h5Name, 'a')
        out = results.create_group(groupName)
        
        out.attrs["reference"] = chrom
        out.attrs["start"] = start
        out.attrs["end"] = end
        
        logging.info("Parsing bam" )
        
        #I'll need to extract this part
        reads = args.bam.fetch(chrom, start, end)
        readCount = args.bam.count(chrom, start, end)
        if readCount == 0:
            logging.warning("No reads found in %s" % groupName)
            continue
        myData = countErrors(reads, start, size, args, readCount)
        
        if size < CHUNKSHAPE[1]:
            chunk = (CHUNKSHAPE[0], size-1)
        else:
            chunk = CHUNKSHAPE
        container = out.create_dataset("data", data = myData, \
                                chunks=chunk, compression="gzip")
        
        if not args.noCallSpots:
            logging.info("Calling Spots")
            spots = callHotSpots(container, start, args)
            fspots = 0
            logging.info("Filtering for read support")
            for spot in spots:
                spot.chrom = chrom
                spot.offset(start)
                
                spot.estimateSize()
                #the sv spans too far
                if spot.size > args.spanMax:
                    continue
                #the sv doesn't have adequate read support
                if supportingReadsFilter(spot, args):
                    continue
                
                if args.noConsensus:
                    spot.estimateSize()#get a better size estimate
                    #Filter based on minimum size expectaions
                    if spot.size < args.minIndelSize or spot.size > args.spanMax:
                        continue
                    fspots += 1
                    if groupName != chrom:
                        spot.tags["RegName"] = groupName
                    hotspots.write(str(spot)+"\n")
                else:##Consensus Calling
                    logging.info("Calling Consensus")
                    s = consensusCalling(spot, args)
                    if len(s) > 0:
                        hotspots.write("\n".join([str(x) for x in s]) + "\n")
                        fspots += len(s)
                    
            logging.info("Found %d spots" % (fspots))
            totSpots += fspots
            hotspots.flush()
            #Down to here for multiprocessing
        results.flush()
        results.close()
        myData = None
        container = None
        gc.collect()
    
    logging.info("Finished %d reads" % (totReads))
    if not args.noCallSpots:
        logging.info("Finished %d spots" % (totSpots))
        hotspots.close()

if __name__ == '__main__':
    run(sys.argv[1:])



###########################
# --- DEPRICATED CODE --- #
###########################
def signalTransform(dat):
    """
    Go through the signal processing changes on a 1d array
    """
    return numpy.convolve(dat, slopWindow, "same") 

def postSigStats(sig):
    mx, mn, std, min = (numpy.max(sig), numpy.mean(sig), numpy.std(sig), \
                        numpy.min(sig))
    logging.info("MaxSig: %f MeanSig: %f StdSig %f MinSig: %f" \
                 % (mx, mn, std, min))
    return mx, mn, std, min

def callHotSpotsOrig(data, offset, args): #threshPct, covThresh, binsize, offset):
    """
    """
    ret = []
    #coverage
    cov = numpy.convolve(data[COV], avgWindow, "same")
    covTruth = numpy.all([cov >= args.minCoverage, cov <= args.maxCoverage], axis=0)
    logging.info("MaxCov:%d MeanCov:%d StdCov:%d MinCov:%d" \
            % (numpy.max(data[COV]), numpy.mean(data[COV]), \
               numpy.std(data[COV]), numpy.min(data[COV])))
    
    #mis
    logging.info("MIS processing")
    mis, mu, sd = preprocessSignal(data[MIS], data[COV])
    mis = signalTransform(mis) 
    mx, mn, sd, min = postSigStats(mis)
    logging.info("Threshold = %.3f" % (sd* args.threshold))
    ret.extend(makeSpotResults(mis, sd, "MIS", cov, covTruth, args))
    del(mis)
    
    #ins
    logging.info("INS processing")
    ins, mu, sd = preprocessSignal(data[INS], data[COV])
    ins = signalTransform(ins)
    mx, mn, sd, min = postSigStats(ins)
    logging.info("Threshold = %.3f" % (sd* args.threshold))
    ret.extend(makeSpotResults(ins, sd, "INS", cov, covTruth, args))
    del(ins)
    
    #dele
    logging.info("DEL processing")
    dele, mu, sd = preprocessSignal(data[DEL], data[COV])
    dele = signalTransform(dele)
    mx, mn, sd, min = postSigStats(dele)
    logging.info("Threshold = %.3f" % (sd* args.threshold))
    ret.extend(makeSpotResults(dele, sd, "DEL", cov, covTruth, args))
    del(dele)
    
    return ret

def markDups(bam):
    """
    Marks all reads that aren't the best read from a zmw as duplicate
    """
    names = {}
    numDups = 0
    for read in bam:
        n = "/".join(read.qname.split('/')[:2])
        try:
            cRead, cScore = names['/'.join(read.qname.split('/')[:2])]
            #myScore = sum([ord(y)-33 for y in read.qqual])/float(len(read.qqual))
            myScore = read.mapq
            if cScore > myScore:
                read.is_duplicate = True
            else:
                cRead.is_duplicate = True
            numDups += 1
        except KeyError:
            #myScore = sum([ord(y)-33 for y in read.qqual])/float(len(read.qqual))
            myScore = read.mapq
            names[n] = (read, myScore)
    logging.info("Marked %d ZMW duplicates" % (numDups))
    del(names)

def expandMd(md):
    """
    Turns abbreviated MD into a full array
    --- C translate candidate --
    """
    ret = []
    for i in re.findall("\d+|\^?[ATCGN]+", md):
        if i[0] == '^':
            d = list(i[1:])
        elif i[0] in ["A","T","C","G","N"]:
            d = list(i)
        else:
            d = xrange(int(i))
        ret.extend(d)
    return ret

def pileupErrors(chrom, start, end, args):
    size = float(end-start)
    container = numpy.zeros( ( len(COLUMNS), int(size) ), dtype=BIGINTY )
    
    size = float(size)
    lastpct = 0
    for base in args.bam.pileup(chrom, start, end):
        if base.pos < start:
            continue
        if base.pos >= end:
            logging.info("Finished region")
            break #short circuit
        
        if lastpct < (base.pos-start)/size:
            logging.info("%.1f%% complete" % (((base.pos-start)/size)*100))
            lastpct += .1
            
        for plup in base.pileups:
            if abs(plup.indel) > args.minIndelErr:
                if plup < 0:
                    s = base.pos
                    e = min(base.pos - plup.indel, container.shape[1])
                    container[DEL, s:e] += BIGINTY(1)
                else:   
                    s = max(0, base.pos - (plup.indel/2))
                    e = min(base.pos + (plup.indel/2), container.shape[1])
                    container[INS, s:e] += BIGINTY(1)
    return container, None


def makeSpotResults(datpoints, sd, label, cov, covTruth, args):
    """
    Find starts and ends ranges, then group the start/end pairs by closest proximity
    Need to explore the edge cases
    need to fix
        V^V ^

    What if the procedure was to:
    1) turn starts into SpotResults and place into sorted list.
    2) for each endPoint, we find the nearest neighbor using insort
        if the start currently doesn't have a nearest neighbor
            if the end is contained within a start, the we throw the end away (hoping another end will be around)
            if the start is contained within an end, we throw away the start and keep searching for the end.
        if the end has a nearest neighbor
            we let the ends compete, the closest end will be put together with the start,
            the losing end will be put through the search 
    This doesnt work
    """
    thresh = sd * args.threshold
    entries = []
    startPoints = makePoints(numpy.all([datpoints <= -thresh, covTruth], axis=0), args.binsize, 's')
    endPoints   = makePoints(numpy.all([datpoints >=  thresh, covTruth], axis=0), args.binsize, 'e')
    sPos = 0
    ePos = 0
        
    points = []
    for i in startPoints:
        bisect.insort(points, i)
    for i in endPoints:
        bisect.insort(points, i)
    
    i = 0
    
    while i < len(points):
        mySpot = SpotResult(tags={"label":label})
        if points[i][2] == 'e':
            mySpot.in_end  = points[i][0]
            mySpot.end     = datpoints[points[i][0]:points[i][1]].argmax() + points[i][0]
            mySpot.out_end = points[i][1]
            mySpot.tags["endCov"] = cov[points[i][0]:points[i][1]].mean()
            mySpot.tags["endSig"] = datpoints[points[i][0]:points[i][1]].mean()
            i += 1
        elif points[i][2] == 's':
            if i+1 < len(points) and points[i+1][2] == 'e':
                mySpot.out_start = points[i][0]
                mySpot.start     = datpoints[points[i][0]:points[i][1]].argmin()+points[i][0]
                mySpot.in_start  = points[i][1]
                mySpot.tags["startCov"] = cov[points[i][0]:points[i][1]].mean()
                mySpot.tags["startSig"] = datpoints[points[i][0]:points[i][1]].mean()
                mySpot.in_end    = points[i+1][0]
                mySpot.end     = datpoints[points[i+1][0]:points[i+1][1]].argmax()+points[i+1][0]
                mySpot.out_end   = points[i+1][1]
                mySpot.tags["endCov"] = cov[points[i+1][0]:points[i+1][1]].mean()
                mySpot.tags["endSig"] = datpoints[points[i+1][0]:points[i+1][1]].mean()
                i += 2
            else:
                mySpot.out_start = points[i][0]
                mySpot.start     = datpoints[points[i][0]:points[i][1]].argmin()+points[i][0]
                mySpot.in_start  = points[i][1]
                mySpot.tags["startCov"] = cov[points[i][0]:points[i][1]].mean()
                mySpot.tags["startSig"] = datpoints[points[i][0]:points[i][1]].mean()
                i += 1
        
        entries.append(mySpot)
        
    logging.info("%d %s entries" % (len(entries), label))
    return entries

def makeKernals(binsize=100):
    """
    My Kernals for Convolution -- push global
    """
    global avgWindow
    global slopWindow
    #slop - downstream average minus upstream average
    slopWindow = numpy.ones(binsize, dtype=numpy.float16) / (binsize/2)
    slopWindow[:binsize/2] = -slopWindow[:binsize/2]
