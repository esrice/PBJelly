#!/usr/bin/env python
import gc
import os
import re
import sys
import copy
import h5py
import math
import time
import numpy
import pysam
import random
import logging
import argparse 
import traceback
import multiprocessing

from contextlib import contextmanager
from tempfile import NamedTemporaryFile
from collections import defaultdict, Counter

from pbsuite.utils.setupLogging import *
from pbsuite.utils.CommandRunner import exe
from pbsuite.banana.Polish import consensus
from pbsuite.utils.FileHandlers import M5File, revComp
from pbsuite.honey.bampie import BLASRPARAMS, EEBLASRPARAMS
from pbsuite.utils.VCFIO import VCFFile, VCFEntry, HONTEMPLATE

VERSION = "15.01.08"
AUTHOR = "Adam English"

USAGE = """\
Detect 'smaller' SVs via measuring discordance between sample and reference in long reads.
"""

#########################
## --- Global Vars --- ##
#########################

#Biggest integer I want to deal with
BIGINT  = 2000 # Greatest Coverage we deal with (limited by BIGINTY)
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
#Size of match stretches that can be ignored
#Explore is for raw reads, confirm is for consenus validation
EXPLORECOLLAPSE = 1
CONFIRMCOLLAPSE = 5
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


#{{{ Adapted fromColby Chiang's svtyper https://github.com/cc2qe/svtyper
#I lower the priors here because it seems reasonable that some reads are
# variant supporting but do not make it through all the filtering I perform
def genotype(spot, avgCov=None, priors=[['0/0',0.05], ['0/1',0.45], ['1/1',0.9]]):
    def log_choose(n, k):
        # swap for efficiency if k is more than half of n
        r = 0.0
        if k * 2 > n:
            k = n - k
            
        for  d in xrange(1,k+1):
            r += math.log(n, 10)
            r -= math.log(d, 10)
            n -= 1
        
        return r
    
    total = int(spot.tags["coverage"]) if avgCov is None else avgCov
    alt = int(spot.tags["szCount"])
    ref = total - alt
    gtList = []
    for gt, p_alt in priors:
        gtList.append(log_choose(total, alt) + alt * math.log(p_alt, 10) + ref * math.log(1 - p_alt, 10))
    gt_idx = gtList.index(max(gtList))
    GL = gtList[gt_idx]
    gtList.remove(GL)
    GT = priors[gt_idx][0]
    GQ = -10 * (GL - sum(10**x for x in gtList))
    
    return GT, GQ
## end of https://github.com/cc2qe/svtyper }}}

def expandCigar(read, minSize, collapse=-1, makeAlt=False):
    """
    Turns the abbreviated cigar into the full array
    0 = M; 1 = I; 2 = D

    collapse set of Ms if they are <= collapse#
    They become 'null' 3s
    and flanks are same error >= collapse
    
    returns hasVar, svs
    where svs equals [(start, size, svtype),...]
    svtypes are MAT INS and DEL
    hasVar says if we found any insertion or deletion. good for short circuiting
    """
    ret = []
    qseq = []
    qpos = 0
    pmat = False
    pmatz = 0
    pdel = False
    pdelz = 0
    pins = False
    pinsz = 0
    
    def pmatLoad(start, size):
        if pmat:
            ret.append((start-size, size, "MAT"))
        return False, 0
        
    def pinsLoad(start, size):
        if pins and size >= minSize:
            if makeAlt:
                ret.append((start, size, "INS", "".join(qseq)))
            else:
                ret.append((start, size, "INS"))
        return False, 0
    
    def pdelLoad(start, size):
        if pdel and size >= minSize:
            if makeAlt:
                ret.append((start-size, size, "DEL", None))
            else:
                ret.append((start-size, size, "DEL"))
        return False, 0
    
    basePosition = read.pos
    for code, size in read.cigar:
        if code == 0:
            pmat = True
            pmatz += 1
            if size > collapse: #break pindel if this isn't collapse noise
                pins, pinsz = pinsLoad(basePosition, pinsz)
                pdel, pdelz = pdelLoad(basePosition, pdelz)
                if makeAlt:
                    qseq = []
            basePosition += size 
            qpos += size
        elif code == 1: #ins
            if size >= minSize:
                pins = True
                pinsz += size
                if makeAlt:
                    qseq.append(read.query[qpos - 1 : qpos - 1 + size])
            pdel, pdelz = pdelLoad(basePosition, pdelz)
            qpos += size
            #pmat, pmatz = pmatLoad(basePosition, pmatz)
        elif code == 2: #del
            if size >= minSize:
                pdel = True
                pdelz += size
                pins, pinsz = pinsLoad(basePosition, pinsz)
                if makeAlt:
                    qseq = []
            #pmat, pmatz = pmatLoad(basePosition, pmatz)
            basePosition += size
        #elif code == 4:
            #qpos += size
        
    return ret

def blasr(query, target, format, nproc = 1, outname = "out.m5", consensus=True):
    """
    Simple mapper
    """
    cmd = ("blasr %s %s %s -nproc %d -bestn 1 -out %s ") \
           % (query, target, format, nproc, outname)
    #need to figure out how to m5-pie it...maybe
    if consensus:
        r, o, e = exe(cmd + " -noSplitSubreads -minMatch 5 " + \
                     "-nCandidates 20 -sdpTupleSize 6 -insertion 1 -deletion 1 -bestn 1")
    else:
        r, o, e = exe(cmd + " -maxAnchorsPerPosition 100 "
               "-affineAlign -affineOpen 100 -affineExtend 0 "
               "-insertion 10 -deletion 10 "
               "-noSplitSubreads -nCandidates 20 ")
        
    logging.debug("blasr - %d - %s - %s" % (r, o, e))

def parseArgs(argv, established=False):
    parser = argparse.ArgumentParser(prog="Honey.py spots", description=USAGE, \
            formatter_class=argparse.RawDescriptionHelpFormatter)
    
    ioGroup = parser.add_argument_group("I/O Arguments")
    ioGroup.add_argument("bam", metavar="BAM", type=str, \
                        help="BAM containing mapped reads")
    ioGroup.add_argument("--hon", metavar="HON.H5", type=str, default=None, \
                        help="HON.h5 containing Error data. Skips ErrorCouting.")
    ioGroup.add_argument("-r", "--region", type=str, default=None,\
                        help="Only call spots in region.bed")
    ioGroup.add_argument("--chrom", type=str, default=None, \
                        help="Only call spots on specified chromosomes (comma-separated) (%(default)s)")
    ioGroup.add_argument("-n", "--nproc", type=int, default=1, \
                        help="Number of processors to use (only for consensus) (%(default)s)")
    ioGroup.add_argument("-o", "--output", type=str, default=None, \
                        help="Basename for output (BAM.hon)")
    ioGroup.add_argument("--readFile", action="store_true", \
                        help="Create a file with what reads support what events (%(default)s)")
                        
    pGroup = parser.add_argument_group("Spot-Calling Threshold/Filtering Arguments")
    pGroup.add_argument("-b", "--binsize", type=int, default=50, \
                        help="Binsize for window averaging (%(default)s)")
    pGroup.add_argument("-e", "--threshold",  type=float, default=3,
                        help="Minimum Spot Threshold (%(default)s)")
    pGroup.add_argument("-c", "--minCoverage", type=int, default=2, \
                        help="Minimum coverage of a region (%(default)s)")
    pGroup.add_argument("-C", "--maxCoverage", type=int, default=BIGINT, \
                        help="Maximum coverage of a region (%(default)s)")
    pGroup.add_argument("-q", "--minMapQ", type=int, default=1, \
                        help="Minimum map quality of reads considered (%(default)s)")
    pGroup.add_argument("-m", "--minIndelErr", type=int, default=5,
                        help="Minimum size of an indel error to be counted (%(default)s)")
    pGroup.add_argument("-i", "--minIndelSize", type=int, default=50, \
                        help="Minimum indel SV size (%(default)s)")
    pGroup.add_argument("-E", "--minErrReads", type=int, default=3, \
                        help="Minimum number of reads with indel (%(default)s)")
    pGroup.add_argument("--spanMax", type=int, default=3000, \
                        help="Maximum Size of spot to be called (%(default)s)")
    #pGroup.add_argument("-I", "--minIndelPct", type=float, default=0.20, \
                        #help="Minimum pct of reads with indel (max(%(default)s*cov,minErrReads)")
    
    aGroup = parser.add_argument_group("Consensus Arguments")
    aGroup.add_argument("--consensus", type=str, default="pbdagcon", choices=["pbdagcon", "pbbanana", "None"], \
                        help="Method for polishing consensus. (%(default)s)")
    aGroup.add_argument("--buffer", default=1000, type=int, \
                        help="Buffer around SV to consense (%(default)s)")
    aGroup.add_argument("--reference", default=None, type=str, \
                        help="Sample reference. Required with consensus calling (None)")
    aGroup.add_argument("--reportContig", action="store_true", \
                        help="Report the contig created that called the spot")
    #aGroup.add_argument("--blasr", default="blasr", \
                        #help="Path to blasr if it's not in the env")
    #aGroup.add_argument("--contig", default="store_false", \
                        #help="Report the full contig sequences and QVs in INFO (False)")
    parser.add_argument("--debug", action="store_true", \
                        help="Verbose logging")
    
    args = parser.parse_args(argv)
    setupLogging(args.debug)
    if args.maxCoverage > BIGINT:
        logging.error("Max Coverge must be less than %d" % (BIGINT))
        exit(0)
    
    #check bam is bamfile
    
    if args.output is None:
        #args.output = args.bam.filename[:-4]+".hon"
        if args.hon is not None:
            args.output = args.hon.rstrip(".h5")
        else:
            args.output = args.bam[:-4]+".hon"
    
    if args.consensus != "None":
        if args.reference is None:
            logging.error("Reference is required with consensus calling")
            exit(0)
        #Check is fastafile
    if args.chrom is not None:
        args.chrom = args.chrom.split(',')
    return args


class SpotResult():
    """
    Represents an SVP style entry with structure
    #CHROM	OUTERSTART	START	INNERSTART	INNEREND	END	OUTEREND	TYPE	SIZE	INFO
    I want to change this to output a .vcf entry
    """
    def __init__(self, chrom=None, start=None,  end=None,  \
                 svtype=None, size=None, tags=None):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.svtype = "UNK" if svtype is None else svtype
        self.size = -1 if size is None else size
        self.tags = tags if tags is not None else {}
        self.varReads = []
        self.varReadsSize = []
        self.refReads = []
   
    @classmethod
    def parseLine(cls, line):
        """
        Turn a line into a spot result
        """
        data = line.strip().split('\t')
        tags = {}
        for i in data[9].split(';'):
            k, v = i.split('=')
            try:
                v = int(v)
            except ValueError:
                try:
                    v = float(v)
                except ValueError:
                    pass
            tags[k] = v
        
        return cls(chrom=data[0], start=int(data[2]), end=int(data[5]), \
                   svtype=data[7], size=int(data[8]), tags=tags)
    
    def offset(self, start):
        """
        moves the spot to an offset
        """
        if self.start is not None:
            self.start += start
        if self.end is not None:
            self.end += start
    
    def fetchbounds(self):
        """
        return (start,end) tuple of spot's boundaries
        """
        return self.start, self.end
    
    def estimateSize(self):
        """
        estimate the sv's size by using the mean or bounds
        """
        if 'szMean' in self.tags:
            self.size = self.tags["szMean"]
        else:
            s,e = self.fetchbounds()
            self.size = e-s
    
    def qregstr(self):
        """
        returns quick region string chrom:start-end
        """
        return "%s:%d-%d" % (self.chrom, self.start, self.end)
    
    def toVCFRecord(self):
        #make
        #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
        CHROM = self.chrom
        POS = self.pos
        ID = '.'
        REF = '.' if self.svtype == "INS" else self.tags["seq"]
        ALT = '.' if self.svtype == "DEL" else self.tags["seq"]
        FILTER = "PASS" if "noSpan" not in self.tags else "NoSpan" #0/0?
        INFO = OrderedDict()
        #Based on what we have, populate info
        
        infos = OrderedDict({"SVTYPE":self.svtype, \
                             "SVLEN":len(self.tags["seq"])})
        
        #myRecord = vcf.model._Record(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO,
        #FORMAT, sample_indexes, samples=None)
   
    def __cmp__(self, other):
        if self.chrom == other.chrom and self.start == other.start and \
          self.end == other.end and self.svtype == other.svtype and \
          self.size == other.size:
            return 0
        return self.start - other.start
        
    def __str__(self):
        """
        changes a spot named tuple to a svp string
        """
        tag = []
        for key in self.tags:
            try:
                tag.append("%s=%0.3f" % (str(key), self.tags[key]))
            except TypeError:
                tag.append("%s=%s" % (str(key), str(self.tags[key])))
        
        
        tag = ";".join(tag)
        dat = [self.chrom, self.start, self.end, self.svtype, self.size, tag]

        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(*dat).replace("None", ".")

class SpotH5():
    """
    Holds a HoneySport h5 file.
    Keeps a multiprocess.Lock to help prevent multiple access
    """
    def __init__(self, filename, version=None, columns=None, parameters=None, mode='r'):
        if mode not in ['r', 'w', 'a']:
            logging.error("Mode must be r(ead), w(rite), or (a)ppend!")
            logging.error("Forcing 'r'")
        self.filename = filename
        if mode == 'w' and os.path.exists(self.filename):
            logging.error("Ouput H5 %s already exists!" % self.filename)
            exit(1)
        if mode == 'r' and not os.path.exists(self.filename):
            logging.error("Output H5 %s doesn't exist!" % self.filename)
            exit(1)
        if mode != 'r':
            self.__reopen(mode)
            self.__results.attrs["version"] = version
            self.__results.attrs["columns"] = columns
            self.__results.attrs["parameters"] = parameters
            self.__close()   
            
        self.__lock = multiprocessing.Lock()

    @contextmanager
    def acquireH5(self, mode='r'):
        with self.__lock:
            self.__reopen(mode)
            yield self.__results
            self.__close()
    
    def lockH5(self):
        """
        You can only lock the h5 in read-only
        """
        self.__lock.acquire()
        self.__results = h5py.File(self.filename, 'r')
        return self.__results

    def releaseH5(self):
        """
        You can only lock the h5 in read-onely
        """
        self.__close()
        self.__lock.release()
        
    def __reopen(self, mode='r'):
        self.__results = h5py.File(self.filename, mode)
        
    def __close(self):
        self.__results.close()
        
###################################
## --- Consumer/Task Objects --- ##
###################################

class Consumer(multiprocessing.Process):
    """
    Basic Consumer. Follow the two queues with your *args and **kwargs that should be sent
    to the task when __call__ 'd

    NOTE! args can't hold anything that isn't pickle-able for the subprocess
    """
    def __init__(self, task_queue, result_queue, bamName, referenceName, honH5): #*args, **kwargs):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.bam = pysam.Samfile(bamName)
        if referenceName is not None:
            self.reference = pysam.Fastafile(referenceName)
        else:
            self.reference = None
        self.honH5 = honH5
        #self.args = args
        #self.kwargs = kwargs
        
    def run(self):
        try: 
            proc_name = self.name
            while True:
                next_task = self.task_queue.get()
                if next_task is None:
                    # Poison pill means shutdown
                    logging.info('Thread %s: Exiting' % proc_name)
                    self.task_queue.task_done()
                    break
                try:
                    next_task(self.bam, self.reference, self.honH5)
                except Exception as e:
                    logging.error("Exception raised in task %s - %s" % (next_task.name, str(e)))
                    
                    exc_type, exc_value, exc_traceback = sys.exc_info()
                    logging.error("Dumping Traceback:")
                    traceback.print_exception(exc_type, exc_value, exc_traceback, file=sys.stdout)

                    next_task.failed = True
                    next_task.errMessage = str(e)
                
                self.result_queue.put(next_task)
                self.task_queue.task_done()
            
            return
        except Exception as e:
            logging.error("Consumer %s Died\nERROR: %s" % (self.name, e))
            return

#I can probably have a generic task class...
class ErrorCounter(): 
    """
    Counts all the errors in the region
    """
    def __init__(self, groupName, chrom, start, end, args):
        self.groupName = groupName
        self.chrom = chrom
        self.start = start
        self.end = end
        self.args = args #I think I'm sending these twice.. whoops
        self.name = groupName + ":ErrorCounter"
        self.failed = False
        self.errMessage = ""
    
    def countErrors(self, reads, offset, regsize, args):
        """
        Sum the errors over any particular reference base
        """
        container = numpy.zeros( ( len(COLUMNS), regsize ), dtype=BIGINTY )
        mqFilt = 0
        progress = .05
        for align in reads:
            if align.mapq < args.minMapQ:
                mqFilt += 1
                continue
            regionStart = 0 if align.pos < offset else align.pos - offset
            regionEnd = regsize if align.aend > (offset + regsize) else align.aend - offset
            container[COV, regionStart:regionEnd] += BIGINTY(1)
            #progress monitoring
            pct = float(regionStart)/container.shape[1]
            if pct >= progress:
                logging.info("%d%% complete" % (pct*100))
                progress += 0.05
                
            for svstart, svsize, svtype in expandCigar(align, args.minIndelErr, EXPLORECOLLAPSE):
                if svstart >= offset and svstart <= offset + regsize:
                    if svtype == "MAT":#depricated
                        begin = max(0, svstart - svsize - offset)
                        end = min(begin - offset, container.shape[1])
                        container[MAT, begin:end] += BIGINTY(1)
                    elif svtype == "DEL":
                        begin = max(0, svstart - offset)
                        end = min(begin + svsize, container.shape[1])
                        container[DEL, begin:end] += BIGINTY(1)
                    elif svtype == "INS":
                        begin = max(0, svstart - (svsize/2) - offset)
                        end = min(begin + (svsize/2), container.shape[1])
                        container[INS, begin:end] += BIGINTY(1)
                     
        logging.debug("%d reads filtered for low mapq in %s" % (mqFilt, self.groupName))
        return container

    def __call__(self, bam, reference, honH5):
        """
        Takes a pysam.Samfile
        """
        logging.info("Starting %s" % (self.name))
        size = self.end - self.start
        regName =  "%s:%d-%d" % (self.chrom, self.start, self.end)
        
        logging.debug("Making container for %s (%s %d bp)" % (self.groupName, regName, size))
                        
        logging.debug("Parsing bam" )
        st = max(0, self.start)
        readCount = bam.count(self.chrom, st, self.end)
        reads = bam.fetch(self.chrom, st, self.end)
        if readCount == 0:
            logging.warning("No reads found in %s" % self.groupName)
            self.failed = True
            self.errMessage = "No reads found in %s" % self.groupName
            return
        else:
            logging.info("%d reads to parse in %s" % (readCount, self.groupName))
        
        myData = self.countErrors(reads, self.start, size, self.args)
        
        #request loc on honH5 and flush results
        with honH5.acquireH5('a') as h5dat:
            out = h5dat.create_group(self.groupName)
            
            out.attrs["reference"] = self.chrom
            out.attrs["start"] = self.start
            out.attrs["end"] = self.end
            
            if size < CHUNKSHAPE[1]:
                chunk = (CHUNKSHAPE[0], size-1)
            else:
                chunk = CHUNKSHAPE
            
            container = out.create_dataset("data", data = myData, \
                                    chunks=chunk, compression="gzip")
            h5dat.flush()

class SpotCaller():
    """
    Takes a full matrix from ErrorCounter and calls/filters spots
    """
    def __init__(self,  groupName, chrom, start, end, args):
        self.groupName = groupName
        self.chrom = chrom
        self.start = start
        self.end = end
        self.args = args
        self.name = groupName + ":SpotCaller"
        self.failed = False
        self.errMessage = ""
        
        #Signal stats are None by default
        self.maxCov = None
        self.avgCov = None
        self.stdCov = None
        self.minCov = None
        
    def preprocessSignal(self, signal, coverage):
        """
        Normalize and print stats returning data and it's std
        """
        #rate = numpy.convolve(signal/coverage, self.avgWindow, "same")
        rate = numpy.convolve(signal, self.avgWindow, "same")
        rate[numpy.any([numpy.isinf(rate), numpy.isnan(rate)], axis=0)] = 0
        mu = numpy.mean(numpy.any([rate != 0], axis=0))
        sd = numpy.std(numpy.any([rate != 0], axis=0))
        logging.info("|%s|RateMean %f  -- RateStd  %f" % (self.name, mu, sd))
        return rate, mu, sd
        
    def callHotSpots(self, data, offset, bam, args): #threshPct, covThresh, binsize, offset):
        """
        """
        ret = []
        self.avgWindow = numpy.ones(args.binsize, dtype=numpy.float16)/float(args.binsize)
        
        #coverage
        cov = numpy.convolve(data[COV], self.avgWindow, "same")
        covTruth = numpy.all([cov >= args.minCoverage, cov <= args.maxCoverage], axis=0)
        
        self.maxCov = numpy.max(data[COV])
        self.avgCov = numpy.mean(data[COV])
        self.stdCov = numpy.std(data[COV])
        self.minCov = numpy.std(data[COV])
        logging.info("|%s|COV processing" % (self.name))
        logging.info("|%s|MaxCov:%d MeanCov:%d StdCov:%d MinCov:%d" \
                % (self.name, self.maxCov, self.avgCov, self.stdCov, self.minCov))
        del(cov)
        
        #ins
        logging.info("|%s|INS processing" % self.name)
        ins, mu, sd = self.preprocessSignal(data[INS], data[COV])
        startPoints = self.makeSpots(ins, mu, sd, covTruth, args.threshold, args.buffer)
        for start, end, zscore in startPoints:
            mySpot = SpotResult(chrom=self.chrom, start=start, end=end, svtype="INS")
            mySpot.tags["groupName"] = self.groupName
            mySpot.offset(self.start)
            if self.supportingReadsFilter(mySpot, bam, args):
                mySpot.tags["zscore"] = zscore
                ret.append(mySpot)
        del(ins)
        
        #dele
        logging.info("|%s|DEL processing" % self.name)
        dele, mu, sd = self.preprocessSignal(data[DEL], data[COV])
        startPoints = self.makeSpots(dele, mu, sd, covTruth, args.threshold, args.buffer)
        for start, end, zscore in startPoints:
            mySpot = SpotResult(chrom=self.chrom, start=start, end=end, svtype="DEL")
            mySpot.tags["groupName"] = self.groupName
            mySpot.offset(self.start)
            if self.supportingReadsFilter(mySpot, bam, args):
                mySpot.tags["zscore"] = zscore
                ret.append(mySpot)
        del(dele)
        
        return ret
        
    def makeSpots(self, data, mu, sd, covTruth, threshold, buffer):
        """
        make the points for the truth set made from the data container
        truth = numpy.array() with boolean values
        return start,end and zscore
        """
        #prevent weirdness
        #truth = numpy.all([data >= mu+threshold*sd, covTruth], axis=0)
        truth = numpy.all([data >= threshold, covTruth], axis=0)
        truth[-1] = False
        shift = numpy.roll(truth, 1)
        
        starts = truth & ~shift
        ends = ~truth & shift
        
        points = zip(numpy.nonzero(starts)[0], numpy.nonzero(ends)[0])
        npoints = []
        if len(points) == 0:
            return npoints
        
        curStart, curEnd = points[0]
        for start, end in points[1:]:
            if start - curEnd <= buffer:
                curEnd = end
            else:
                zscore = (mu - data[curStart:curEnd].mean())/sd
                npoints.append((curStart, curEnd, zscore))
                curStart = start
                curEnd = end
        
        zscore = (mu - data[curStart:curEnd].mean())/sd
        npoints.append((curStart, curEnd, zscore))
        return npoints
    
    def supportingReadsFilter(self, spot, bam, args):
        """
        filters insertions or deletion spots based on errors
        """
        if spot.svtype == "INS":
            errId = 1
            errLab = 'insertion'
        elif spot.svtype == "DEL":
            errId = 2
            errLab = 'deletion'
        else:#don't worry about other types
            return True
        
        begin, ending = spot.fetchbounds()
        buff = abs(begin-ending) * 2#area we allow read errors to exist
        
        reads = bam.fetch(str(spot.chrom), max(0, begin - 1000), ending + 1000)
        totSizes = []
        coverage = 0
        nReadsErr = 0
        #For tandem
        strandCnt = {True: 0, False: 0}
        
        #count reads and errSizes
        mqFilt = 0
        for read in reads:
            #must span -- still worried about this
            if not (read.pos < begin and read.aend > ending):
                continue
            if read.mapq < args.minMapQ: #mq filt
                mqFilt += 1
                continue

            coverage += 1
            readHasErr = False
            totErrSize = 0
            for svstart, svsize, svtype in expandCigar(read, args.minIndelErr, EXPLORECOLLAPSE):
                if svtype != spot.svtype:
                    continue
                #contained
                if  svstart >= begin - buff and svstart < ending + buff:
                    totErrSize += svsize
                #overlaps
                elif svtype == 'DEL' and begin-buff <= svstart + svsize < ending+buff:
                    totErrSize += svsize
                
                
            if totErrSize >= args.minIndelSize:
                readHasErr = True
                nReadsErr += 1
                totSizes.append(totErrSize)
                strandCnt[read.is_reverse] += 1
                spot.varReads.append(read.qname)
                spot.varReadsSize.append(totErrSize)
            else:
                spot.refReads.append(read.qname)
                
        spot.tags["mqfilt"] = mqFilt
        spot.tags["strandCnt"] = "%d,%d" % (strandCnt[False], strandCnt[True])
        if len(totSizes) < args.minErrReads:
            return False # false - you should filter
        
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
        
        spot.tags["coverage"] = coverage
        spot.tags["szCount"]  = int(nReadsErr)
        spot.tags["szMean"]   = int(mean)
        spot.tags["szMedian"] = int(median)
        spot.tags["sz1stQ"]   = int(firstQ)
        spot.tags["sz3rdQ"]   = int(thirdQ)
        return True
    
    def __call__(self, bam, reference, honH5):
        """
        """
        logging.info("Starting %s" % (self.name))
        with honH5.acquireH5('r') as h5dat:
            #if self.end - self.start != len(h5dat[self.groupName]["data"]
            #I don't know how to subregion a full h5 anymore...
            #if self.start == 0 and self.end == len(h5dat[self.groupName]["data"][0,:]):
            #myData = numpy.array(h5dat[self.groupName]["data"])[:,self.start:self.end])
            myData = numpy.array(h5dat[self.groupName]["data"])#[:,self.start:self.end])
        
        self.calledSpots = self.callHotSpots(myData, self.start, bam, self.args)
        
        with honH5.acquireH5('a') as h5dat:
            if "max_coverage" not in h5dat[self.groupName].attrs:
                #I don't want to overwrite incase I'm recalling over a region
                h5dat[self.groupName].attrs["max_coverage"] = self.maxCov
                h5dat[self.groupName].attrs["avg_coverage"] = self.avgCov
                h5dat[self.groupName].attrs["std_coverage"] = self.stdCov
                h5dat[self.groupName].attrs["min_coverage"] = self.minCov

class ConsensusCaller():
    """
    For any particular spot, create a consensus
    """
    def __init__(self, spot, args):
        self.spot = spot
        self.args = args
        self.name = spot.qregstr() + ":ConsensusCaller"
        self.failed = False
        self.errMessage = ""
        
    def readTrim(self, read, start, end):
        """
        Trims a pysam.AlignedRead to only include the sequence that's aligned (or should be aligned)
        between start and end on reference
        returns the sequence and quality
        """
        score = 0
        if not read.is_unmapped:
            regTrim = 0
            #upS = read.cigar[0][1] if read.cigar[0][0] == 4 else 0
            #dnS = read.cigar[-1][1] if read.cigar[-1][0] == 4 else 0
            
            trimS = None
            trimE = None
            
            if start > read.pos:
                for queryPos, targetPos in read.aligned_pairs:
                    
                    if trimS is None and targetPos is not None and targetPos >= start:
                        trimS = queryPos
            else:
                score += abs(read.pos - start)
            if end < read.aend:
                for queryPos, targetPos in read.aligned_pairs[::-1]:
                    if trimE is None and targetPos is not None and targetPos <= end:
                        trimE = queryPos
            else:
                score += abs(read.aend-end)
            
            if trimS is not None:
                #trimS = max(0, trimS) + upS
                trimS = max(0, trimS)
            else:
                trimS = 0
                    
            if trimE is not None:
                #trimE = min(len(read.seq), trimE)  - dnS
                trimE = min(len(read.seq), trimE)  
            else:
                trimE = len(read.seq)
            seq = read.query[trimS:trimE]
            qual = read.qqual[trimS:trimE]
            if not read.is_reverse:
                seq = seq.translate(revComp)[::-1]
                qual = qual[::-1]
        
        return seq, qual
    
    def consensusCalling(self, spot, bam, reference, args):
        """
        Make a consensus of all the reads in the region and identify all of the SVs in the region
        """
        #
        MAXNUMREADS = 100 #I don't think we'll need more than this many reads
        MAXATTEMPTS = 5   #MAXNUMREADS/2 #I don't feel like trying 100 times
        SPANBUFFER  = 100 #number of bases I want a read to span
        
        chrom, start, end = spot.chrom, spot.start, spot.end
        buffer = args.buffer
        
        supportReads = []
        spanReads = []
        #Fetch reads and trim
        totCnt = 0
        for read in bam.fetch(chrom, max(0, start-buffer-SPANBUFFER), end+buffer+SPANBUFFER):
            if read.qname not in spot.varReads:
                continue
            seq, qual = self.readTrim(read, start-buffer, end+buffer)
            if read.pos < start-SPANBUFFER and read.aend > end+SPANBUFFER:
                sz = spot.varReadsSize[spot.varReads.index(read.qname)]
                spanReads.append((abs(sz - spot.tags["szMedian"]), seq, qual))
            else:
                supportReads.append((seq, qual))
            totCnt += 1
            
        if len(spanReads) == 0:
            logging.debug("noone spans - consensus aborted. %s" % (str(spot)))
            spot.tags["noSpan"] = True
            return [spot]
            
        #spanReads.sort(reverse=True)
        spanReads.sort()
        if len(spanReads) > MAXNUMREADS:
            origSupportReads = [(x[1], x[2]) for x in spanReads[:MAXNUMREADS]]
        elif len(spanReads) + len(supportReads) > MAXNUMREADS:
            origSupportReads = [(x[1], x[2]) for x in spanReads] + supportReads[:MAXNUMREADS-len(spanReads)]
        else:
            origSupportReads = [(x[1], x[2]) for x in spanReads] + supportReads
        mySpots = []
        refReadId = 0
        haveVar = False
        #Attempt each spanRead until we get one that passes
        while refReadId < len(spanReads) and not haveVar and refReadId < MAXATTEMPTS:
            refread = spanReads[refReadId]
            supportReads = origSupportReads[:refReadId] + origSupportReads[refReadId+1:] 
            refReadId += 1
            
            #read that spans most of the region goes first
            #use the rest for cleaning
            
            #building consensus sequence
            foutreads = NamedTemporaryFile(suffix=".fastq")
            for id, i in enumerate(supportReads):
                foutreads.write("@%d\n%s\n+\n%s\n" % (id, i[0], i[1]))
            foutreads.flush()
            foutref = NamedTemporaryFile(suffix=".fasta")
            foutref.write(">%s:%d-%d\n%s" % (spot.chrom, start, end, refread[1]))
            foutref.flush()
            
            alignOut = NamedTemporaryFile(suffix=".m5")
            logging.debug("making the contig....")
            #run it through phrap

            #then run it through consensus
            blasr(foutreads.name, foutref.name, format="-m 5", nproc=1, outname=alignOut.name)
            if args.consensus == "pbbanana":
                aligns = M5File(alignOut.name)
                con = ">con\n%s\n" % consensus(aligns).sequence
                conName = "pbbanana"
            elif args.consensus == "pbdagcon":
                logging.debug("pbdagcon is running")
                #using minerreads - 1 because one f them is already being used as seed!
                #I want to be sure I get something out... so just require somebody on there
                #r, con, e = exe("pbdagcon -c %d -t 0 %s" % (1, alignOut.name), timeout=1)
                #r, con, e = exe("pbdagcon -m 100 -c %d -t 0 %s" % (max(args.minErrReads - 1, 0), alignOut.name), timeout=1)
                r, con, e = exe("pbdagcon -m 100 -c %d -t 0 %s" % (3, alignOut.name), timeout=1)
                logging.debug("back from pbdagcon")
                logging.debug((r,e))
                #raw_input("press ent")
                if con is not None:
                    con = con[con.index("\n")+1:]
                else:
                    con = ""
                conName = "pbdagcon"
            alignOut.close()
            foutref.close()
            foutreads.close()
            #we don't have a consensus - retry
            if len(con) == 0:
                logging.debug("Trying another seed read for consensus")
                continue
            logging.debug("%s %d bp seq" % (conName, len(con.split('\n')[1])))
            
            #try improving consensus
            conOut = NamedTemporaryFile(suffix=".fasta")
            conOut.write(con)
            #conOut.close()
            conOut.flush()
            
            refOut = NamedTemporaryFile(suffix=".fasta")
            #j = reference.fetch(chrom, max(0, start-buffer), end+buffer)
            #fout = open("fuckme.ref.fasta",'w')
            #fout.write(j)
            #fout.close()
            refOut.write(">%s:%d-%d\n%s\n" % (chrom, start, end, \
                        reference.fetch(chrom, max(0, start-(buffer*2)), end+(buffer*2))))
            refOut.flush()
            
            #map consensus to refregion
            varSam = NamedTemporaryFile(suffix=".sam")
            blasr(conOut.name, refOut.name, format="-sam", outname=varSam.name,\
                consensus=False) #-- would this help?
                #or what if I fed it through leftalign?
            #os.system("cp %s ." % (refOut.name))
            #os.system("cp %s ." % (varSam.name))
            sam = pysam.Samfile(varSam.name)
            
            matches = 0.0
            bases = 0.0
            nReads = 0
            minVarDiff = 10000
            for read in sam:
                localSpots = []
                nReads += 1
                spot.tags["consensusCreated"] = True
                for svstart, svsize, svtype, altseq in expandCigar(read, args.minIndelSize, CONFIRMCOLLAPSE, True):
                    newspot = copy.deepcopy(spot)
                    
                    if spot.svtype == svtype and svtype == "INS":
                        #haveVar = True
                        newspot.start = svstart + start - (buffer*2)
                        newspot.end = svstart + start - (buffer*2)
                        newspot.tags["seq"] = altseq
                        newspot.size = svsize
                        gt, gq = genotype(newspot)
                        newspot.tags["GT"] = gt
                        newspot.tags["GQ"] = gq
                        if abs(spot.tags["szMedian"] - newspot.size) < minVarDiff:
                            minVarDiff = abs(spot.tags["szMedian"] - newspot.size)
                        if args.reportContig:
                            newspot.tags["contigseq"] = read.seq
                            newspot.tags["contigqual"] = read.qual
                        localSpots.append(newspot)
                    
                    elif spot.svtype == svtype and svtype == "DEL":
                        #haveVar = True
                        newspot.start = svstart + start - (buffer*2)
                        newspot.end = svstart + svsize + start - (buffer*2)
                        newspot.size = svsize
                        gt, gq = genotype(newspot)
                        newspot.tags["GT"] = gt
                        newspot.tags["GQ"] = gq
                        newspot.tags["seq"] = reference.fetch(chrom, newspot.start, newspot.end)
                        if abs(spot.tags["szMedian"] - newspot.size) < minVarDiff:
                            minVarDiff = abs(spot.tags["szMedian"] - newspot.size)
                        if args.reportContig:
                            newspot.tags["contigseq"] = read.seq
                            newspot.tags["contigqual"] = read.qual
                        localSpots.append(newspot)
                if len(localSpots) > 0:
                    mySpots.append((minVarDiff, localSpots))
            
            #identity = matches/bases
            #If no var, nothing is returned.
            #for newspot in mySpots:
                #newspot.tags["alnIdentityEstimate"] = identity
                #Keep reporting the actual contigs out until we 
                #find a reason to need it (and also we can get quals...)
                #vbam.reset()
                #for id, read in enumerate(vbam):
                    #newspot.tags["contigSeq%d" % (id)] = read.seq 
                    #newspot.tags["contigQual%d" % (id)] = read.qual 
            
            #vbam.close()
            #varBam.close()
            refOut.close()
            
            #logging.debug("%d consensus reads created %d spots" % (nReads, len(localSpots)))

        if len(mySpots) == 0:
            return []
        
        mySpots.sort()
        return mySpots[0][1]
    
    def __call__(self, bam, reference, honH5):
        """
        """
        logging.info("Starting %s" % (self.name))
        self.newSpots = self.consensusCalling(self.spot, bam, reference, self.args)

############################
## --- Execution Code --- ##
############################

def run(argv):
    numpy.seterr(all="ignore")
    args = parseArgs(argv)
    
    bam = pysam.Samfile(args.bam)
    try:
        if bam.header["HD"]["SO"] != "coordinate":
            logging.warning("BAM is not sorted by coordinates! Performance may be slower")
    except KeyError:
        logging.warning("Assuming BAM is sorted by coordinate. Be sure this is correct")
   
    hotSpots = open(args.output+".spots", 'w')
    hotSpots.write("#CHROM\tOUTERSTART\tSTART\tINNERSTART\tINNEREND\tEND\tOUTEREND\tTYPE\tSIZE\tINFO\n")
    if args.readFile:
        readFile = open(args.output + ".reads", 'w')
        readFile.write("#CHROM\tOUTERSTART\tSTART\tINNERSTART\tINNEREND\tEND\tOUTEREND\tTYPE\tSIZE\tINFO\n")
        
    regions = {}
    if args.region:
        fh = open(args.region,'r')
        for line in fh.readlines():
            data = line.strip().split('\t')
            if args.chrom is None or data[0] in args.chrom:
                regions[data[3]] = (data[3], data[0], int(data[1]), int(data[2]))
                #regions.append((data[3], data[0], int(data[1]), int(data[2])))
        fh.close()
    else:
        for chrom, size in zip(bam.references, bam.lengths):
            if args.chrom is None or chrom in args.chrom:
                regions[chrom] = (chrom, chrom, 0, size)
    
        
    tasks = multiprocessing.JoinableQueue();
    results = multiprocessing.Queue();
    
    if args.hon is not None:
        honH5 = SpotH5(args.hon)
        gotoQueue = results
    else:
        honH5 = SpotH5(args.output + ".h5", \
                       VERSION, \
                       COLUMNS, \
                       str(args), mode='w')
        gotoQueue = tasks

    #ErrorCounting
    errConsumers = [ Consumer(tasks, results, args.bam, args.reference, honH5) for i in xrange(args.nproc) ]
    for w in errConsumers:
        w.start()
           
    num_jobs = 0
    for groupName, chrom, start, end in regions.values():
        gotoQueue.put(ErrorCounter( groupName, chrom, start, end, args))
        num_jobs += 1
               
    while num_jobs:
        result = results.get()
        num_jobs -= 1
        
        if result.failed:
            logging.error("Task %s Failed (%s)" % (result.name, result.errMessage))
        
        elif result.name.endswith("ErrorCounter"):#counted -> call
            logging.info("%s -> spot" % (result.name))
            groupName, chrom, start, end = regions[result.groupName]
            tasks.put(SpotCaller( groupName, chrom, start, end, args ))
            num_jobs += 1
            
        elif result.name.endswith("SpotCaller"):#call -> consensus
            nspot = 0
            for spot in result.calledSpots:
                spot.estimateSize()
                if spot.size < args.minIndelSize or spot.size > args.spanMax:
                    continue 
                
                nspot += 1
                if args.consensus == "None":
                    hotSpots.write(str(spot) + '\n')
                    if args.readFile:
                        readFile.write(str(spot) + '\n')
                        if len(x.varReads) > 0:
                            readFile.write("\n".join(['var ' + j for j in x.varReads]) + '\n')
                        if len(x.refReads) > 0:
                            readFile.write("\n".join(['ref ' + j for j in x.refReads]) + '\n')
                else:
                    tasks.put(ConsensusCaller(spot, args))
                    num_jobs += 1
            if nspot > 0:
                logging.info("%s -> %d consensus" % (result.name, nspot))
            
        elif result.name.endswith("ConsensusCaller"):#consensus -> finish
            logging.info("Task %s -> finish" % (result.name))
            for x in result.newSpots:
                hotSpots.write(str(x)+'\n')
                if args.readFile:
                        readFile.write(str(x) + '\n')
                        if len(x.varReads) > 0:
                            readFile.write("\n".join(['var ' + j for j in x.varReads]) + '\n')
                        if len(x.refReads) > 0:
                            readFile.write("\n".join(['ref ' + j for j in x.refReads]) + '\n')
    
    #Poison the Consumers.. I'm done with them
    for i in xrange(args.nproc):
        tasks.put(None)
    
    logging.info("Finished")

def test(argv):
    numpy.seterr(all="ignore")
    args = parseArgs(argv)
    setupLogging(True)#keep debug on.. you're testing!
    logging.critical(("Running HSpots.py directly implements testing mode. "
                      "If you're trying to run the full, actual program, use "
                      "Honey.py spots"))
       
    bam = pysam.Samfile(args.bam)
    reference = pysam.Fastafile(args.reference)
    try:
        if bam.header["HD"]["SO"] != "coordinate":
            logging.warning("BAM is not sorted by coordinates! Performance may be slower")
    except KeyError:
        logging.warning("Assuming BAM is sorted by coordinate. Be sure this is correct")
    logging.info("Running in test mode")
    
    #do what you will.. from here
    # This is what I need to start with
    #spot = SpotResult(chrom="7", start=138402727, end=138402830, svtype="INS", size=113)
    chrom="3"      
    start,end = (195498264, 195498609)
    start -=200
    end +=200
    spot = SpotResult(chrom=chrom, start=start, end=end, svtype="DEL", size=100)
    
    #fh = open("possible.bed")
    #for line in fh.readlines():
        #data = line.strip().split('\t')
        #spot = SpotResult(chrom=data[0], start=int(data[8]), end = int(data[9]), \
                          #size=int(data[5]), svtype=data[4])
        
    j = SpotCaller('group', spot.chrom, spot.start, spot.end, args)
    if j.supportingReadsFilter(spot, bam, args):
        consen = ConsensusCaller(spot, args)
        consen(bam, reference, 'none')
        for i in consen.newSpots:
            i.tags["seqmade"] = True
            print i
        if len(consen.newSpots) == 0:
            spot.tags["noseq"] = True
            print str(spot)
    else:
        spot.tags["filtfail"] = True
        print str(spot)
    #done with test code
    logging.info("Finished testing")
    
if __name__ == '__main__':
    #run(sys.argv[1:]) 
    test(sys.argv[1:])

