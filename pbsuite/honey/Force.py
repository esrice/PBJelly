#!/usr/bin/env python
import sys, argparse, re
from collections import defaultdict, namedtuple, Counter
import pysam
import pbsuite.honey.TGraf as tails
import pbsuite.honey.HSpots as spots
from pbsuite.utils.setupLogging import *

vtypes = ['INS', 'DEL', 'MIS', 'UNK', 'CTX', 'ITX'] 
USAGE = """
Checks if there are any reads in the given .bam that support predicted SVs

Takes a .bed with the first 6 columns being:
    chrom  start  end  name  svtype size

- svtype must be one of DEL, INS, MIS
- size is what the SV's size is estimated to be.  

If you have single breakpoint events (such as translocations) specify --bedPE
Your input.bed's 9 columns become:
    chrom1  start1  end1 orient1  chrom2  start2  end2 orient2  name  svtype size 

- orient is the directionality of the sequence leading upto the breakpoint (+/-)

RegionBuffer is the +- space in which you consider reads for support 
    around predicted sv

SizeBuffer is the +-  percent of predicted size the read needs to support the SV

Results are an extra column appended to the end of in the format REF[TAILS|SPOTS]
    REF:
        True if we found evidence of the reference over the region
        False if we had the opportunity to support the reference, but didn't.
        ? if we didn't have the opportunity to support the reference
    TAILS/SPOTS:
        A comma-separated list of chrBPSchr:start-end(svtype)size*cnt coordinates for
        reads that have interrupted or discordant mapping support of the SV.

TAILS/SPOTS:
    chr         The chromosome
    BPS         Breakpoint string showing orientations
    start/end   The breakpoints start and end coordinates
    svtype      One of %s
    cnt         Number of reads that support this
""" % ("/".join(vtypes))



class Variant():
    def __init__(self, chrom, start, end, svtype, size, read=None):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.svtype = svtype
        self.size = size
        self.read = read
    
    def __str__(self):
        return "%s:%d-%d(%s)%d" % (self.chrom, self.start, self.end, self.svtype, self.size)
    
class BedEntry():
    def __init__(self, chrom, start, end, name, svtype, size, *args):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.svtype = svtype
        self.size = int(size)
        self.rest = args
        if self.start > self.end:
            self.start, self.end = self.end, self.start
    
    def __str__(self):
        return "\t".join([str(x) for x in [self.chrom, self.start, \
                                           self.end, self.name, \
                                           self.svtype, self.size, \
                                           "\t".join(self.rest)]])

    def __repr__(self):
        return "<BedEntry '%s'>" % (str(self).replace('\t',' '))

class BedPEEntry():
    def __init__(self, chrom1, start1, end1, orient1, chrom2, start2, end2, orient2, name, svtype, size, *args):
        self.chrom1 = chrom1
        self.start1 = int(start1)
        self.end1 = int(end1)
        self.orient1 = orient1
        self.chrom2 = chrom2
        self.start2 = int(start2)
        self.end2 = int(end2)
        self.orient2 = orient2
        self.name = name
        self.svtype = svtype
        self.size = int(size)
        self.rest = args

    def __str__(self):
        return "\t".join([str(x) for x in [self.chrom1, self.start1, \
                                           self.end1, self.orient1, \
                                           self.chrom2, self.start2, \
                                           self.end2, self.orient2, \
                                           self.name, self.svtype, \
                                           self.size, "\t".join(self.rest)]])
    
    def __repr__(self):
        return "<BedPEEntry '%s'>" % (str(self).replace('\t',' '))

def parseArgs(args):
    parser = argparse.ArgumentParser(prog="Honey.py force", description=USAGE, \
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bam", metavar="BAM", type=str, \
                        help="Assembled Contigs Bam")
    parser.add_argument("bed", metavar="BED", type=str, \
                        help="Bed of locations to force SV Calls")
    parser.add_argument("-s", "--sizebuffer", type=float, default=0.35, \
                        help=("Buffer of estimated sv size to "
                              "create match (%(default)s)"))
    parser.add_argument("-d", "--maxDelta", type=int, default=500, \
                        help="Max distance between predicted and discovered variant (%(default)s)")
    parser.add_argument("-f", "--fetchbuffer", type=int, default=1000, \
                        help="Buffer for fetching reads from .bam (%(default)s)")
    #parser.add_argument("-o", "--overlapbuffer", type=float, default=0.50, \
                        #help="Percent overlap required from calls to tails (%(default)s)")
    parser.add_argument("-q", "--minMapq", type=int, default=100, \
                        help="Minimum mapping quality of a read and it's tail to consider (%(default)s)")
    parser.add_argument("-m", "--minErr", type=int, default=5, \
                        help="Minimum ins/del error size to consider (%(default)s)")
    #parser.add_argument("-a", "--asm", action="store_true", \
                        #help="Input reads are high-quality contigs")
    parser.add_argument("-p", "--bedPE", action="store_true", \
                        help="Input bed file is bedPE - only tails searching will be performed")
    parser.add_argument("--debug", action="store_true", \
                        help="Verbose logging")
    args = parser.parse_args(args)
    setupLogging(args.debug)
    
    return args

#MonkeyWrenching of sorts
class FakeBread(tails.Bread):
    
    def __init__(self, svtype, refKey, uRef, dRef, uBreak, dBreak, uDir, \
                 dDir, is_reverse, uTail, dTail, read):
        self.svtype = svtype
        self.refKey = refKey
        self.uRef = uRef
        self.dRef = dRef
        self.uBreak = uBreak
        self.dBreak = dBreak
        self.uDir = uDir
        self.dDir = dDir
        self.is_reverse = is_reverse
        self.uTail = uTail
        self.dTail = dTail
        self.read = read
        self.isInverted = False # This is a problem
        self.remainSeq = 1000
    
    def annotate(self):
        return tails.Bread.annotate(self)

Fr = namedtuple("fakeread", "is_reverse")
def bedpeToFakey(bed):
    """
    Make a fakeBread for PE bed entry
    if bed.svtype is unk, we have to make all orientations
    """
    # These are all the element I'll be expected to make
    # refKey, uBreak, dBreak, uDir, dDir, is_reverse, uTail, dTail
    uRef = bed.chrom1; dRef = bed.chrom2
    j = [uRef, dRef]; j.sort(); refKey = "_".join(j)
    ret = []
    if bed.svtype == 'UNK':
        #Need to make all orientations
        #3i 3e    +/-
        #+/- to pick the breakpoint locations
        bp1 = bed.start1 if True else bed.end1
        bp2 = bed.start2 if True else bed.end2
        #choose orientation of the first read and automatically flip the second
        orient1 = '5' if bed.orient1 == '-' else '3'
        orient2 = '5' if bed.orient2 == '+' else '3'
        j = [(bp1, orient1, bed.chrom1), (bp2, orient2, bed.chrom2)]; j.sort(); uBP1, dBP2 = j
        uBreak, uDir, uRef = uBP1
        dBreak, dDir, dRef = dBP2
        j = [uRef, dRef]; j.sort(); refKey = "_".join(j)
        is_reverse = bed.orient1 == bed.orient2 
        fr = Fr(False)
        ret.append(FakeBread(bed.svtype, refKey, uRef, dRef, uBreak, dBreak, uDir, dDir, is_reverse, 'i', 'e', fr))
        
        bp1 = bed.start1 if True else bed.end1
        bp2 = bed.start2 if False else bed.end2
        #choose orientation of the first read and automatically flip the second
        orient1 = '5' if bed.orient1 == '-' else '3'
        orient2 = '5' if bed.orient2 == '+' else '3'
        j = [(bp1, orient1, bed.chrom1), (bp2, orient2, bed.chrom2)]; j.sort(); uBP1, dBP2 = j
        uBreak, uDir, uRef = uBP1
        dBreak, dDir, dRef = dBP2
        j = [uRef, dRef]; j.sort(); refKey = "_".join(j)
        is_reverse = bed.orient1 == bed.orient2 
        fr = Fr(False)
        ret.append(FakeBread(bed.svtype, refKey, uRef, dRef, uBreak, dBreak, uDir, dDir, is_reverse, 'i', 'e', fr))

        bp1 = bed.start1 if False else bed.end1
        bp2 = bed.start2 if True else bed.end2
        #choose orientation of the first read and automatically flip the second
        orient1 = '5' if bed.orient1 == '-' else '3'
        orient2 = '5' if bed.orient2 == '+' else '3'
        j = [(bp1, orient1, bed.chrom1), (bp2, orient2, bed.chrom2)]; j.sort(); uBP1, dBP2 = j
        uBreak, uDir, uRef = uBP1
        dBreak, dDir, dRef = dBP2
        j = [uRef, dRef]; j.sort(); refKey = "_".join(j)
        is_reverse = bed.orient1 == bed.orient2 
        fr = Fr(False)
        ret.append(FakeBread(bed.svtype, refKey, uRef, dRef, uBreak, dBreak, uDir, dDir, is_reverse, 'i', 'e', fr))

        bp1 = bed.start1 if False else bed.end1
        bp2 = bed.start2 if False else bed.end2
        #choose orientation of the first read and automatically flip the second
        orient1 = '5' if bed.orient1 == '-' else '3'
        orient2 = '5' if bed.orient2 == '+' else '3'
        j = [(bp1, orient1, bed.chrom1), (bp2, orient2, bed.chrom2)]; j.sort(); uBP1, dBP2 = j
        uBreak, uDir, uRef = uBP1
        dBreak, dDir, dRef = dBP2
        j = [uRef, dRef]; j.sort(); refKey = "_".join(j)
        is_reverse = bed.orient1 == bed.orient2 
        fr = Fr(False)
        ret.append(FakeBread(bed.svtype, refKey, uRef, dRef, uBreak, dBreak, uDir, dDir, is_reverse, 'i', 'e', fr))
    else:
        #PE -> <-
        #+/- to pick the breakpoint locations
        bp1 = bed.start1 if bed.orient1 == '-' else bed.end1
        bp2 = bed.start2 if bed.orient2 == '-' else bed.end2
        #choose orientation of the first read and automatically flip the second
        orient1 = '5' if bed.orient1 == '-' else '3'
        orient2 = '5' if bed.orient2 == '+' else '3'
                
        j = [(bp1, orient1, bed.chrom1), (bp2, orient2,bed.chrom2)]; j.sort(); uBP1, dBP2 = j
        uBreak, uDir, uRef = uBP1
        dBreak, dDir, dRef = dBP2
        
        is_reverse = bed.orient1 == bed.orient2 
        fr = Fr(is_reverse)
        
        #My main read I'll fake cluster with everything
        ret.append(FakeBread(bed.svtype, refKey, uRef, dRef, uBreak, dBreak, uDir, dDir, is_reverse, 'i', 'e', fr))
    
    #logging.debug(fakey)
    return ret

def bedToFakey(bed):
    """
    Make a fakeBread for bed entry
    returns a list.
    If bed.svtype == 'UNK',
        then a fakey for every orientation
        combination will be created
    """
    ret = []
    uRef = bed.chrom; dRef = bed.chrom
    j = [uRef, dRef]; j.sort(); refKey = "_".join(j)
    
    #PE -> <-
    #+/- to pick the breakpoint locations
    fr = Fr(False)
    if bed.svtype == 'DEL' or bed.svtype == 'UNK':
        ret.append(FakeBread(bed.svtype, refKey, uRef, dRef, bed.start, bed.end, '3', '3', False, 'i', 'e', fr))
    elif bed.svtype == 'INS' or bed.svtype == 'UNK':
        ret.append(FakeBread(bed.svtype, refKey, uRef, dRef, bed.start, bed.end, '3', '3', False, 'e', 'i', fr))
    elif bed.svtype == 'MIS' or bed.svtype == 'UNK':
        ret.append(FakeBread(bed.svtype, refKey, uRef, dRef, bed.start, bed.end, '3', '5', False, 'i', 'e', fr))
    
    return ret

def removeRedundantReads(reads):
    """
    check readname, chrom/pos, flag to see if you're looking at the same
    reads
    """
    ret = []
    while len(reads) > 0:
        cur = reads.pop()
        rm = []
        for i in range(len(reads)):
            cmp = reads[i]
            if cur.compare(cmp) == 0:
                rm.append(cmp)
        ret.append(cur)
        for i in rm:
            reads.remove(i)
    return ret
            
#refactored
def tailsSearch(bam, fakey, args):
    """
    Look in bam for reads that support the predicted structural variant
    """
    reads = []
    #First breakpoint - 
    fetchS = max(0, fakey.uBreak - args.fetchbuffer)
    fetchE = min(fakey.uBreak + args.fetchbuffer, bam.lengths[bam.references.index(fakey.uRef)])
    if fetchS > fetchE:#don't know why this happened
        fetchS, fetchE = fetchE, fetchS
    reads.extend([x for x in bam.fetch(fakey.uRef, fetchS, fetchE)])
    
    #Second breakpoint -
    fetchS = max(0, fakey.dBreak - args.fetchbuffer)
    fetchE = min(fakey.dBreak + args.fetchbuffer, bam.lengths[bam.references.index(fakey.dRef)])
    if fetchS > fetchE:#don't know why this happened
        fetchS, fetchE = fetchE, fetchS
    reads.extend([x for x in bam.fetch(fakey.dRef, fetchS, fetchE)])
    
    #It's possible that the same reads are fetched from both breakpoints
    #We need to remove those redundant reads
    reads = removeRedundantReads(reads)
    
    #This is the evidence we have
    points, tlocs = tails.parseBreakReads(reads, getrname = bam.getrname)
    
    #now we need to keep track of who is near our fakey
    nears = []
    for key in points:
        for read in points[key]:
            #search for the original
            if read.near(fakey):
                nears.append(read.toBriefString())
    
    return len(reads) > 1, nears
    
def oldTailsSearch(bam, bed, args):
    """
    Populate the answer dictionary by looking for tails 
    through the bam

    Returns a list of pbsuite.honey.TGraf.Bnode that support
    """
    fetchS = max(0, bed.start - args.fetchbuffer)
    fetchE = min(bed.end + args.fetchbuffer, bam.lengths[bam.references.index(bed.chrom)])
    
    points = tails.makeBreakReads(bam.fetch(bed.chrom, fetchS, fetchE), getrname = bam.getrname)
    
    reads = []
    anyCoverage = False
    
    for key in points:
        anyCoverage = True
        #eventually will need tloc work
        if key.split('_')[0] != bed.chrom:
            continue
        #eventually will need a reference allele check for tails
        for read in points[key]:
            anno = read.annotate()
            if anno in ['TLOC', 'INV']:
                anno = 'MIS'
            
            #TLOCs...
            if bed.chrom != bam.getrname(read.read.tid):
                continue
            
            if anno != bed.svtype:
                #Not perfect..
                continue
            
            #within reciprocal ovl
            maxStart = max(bed.start, read.uBreak)
            minEnd   = min(bed.end, read.dBreak)
            if minEnd <= maxStart: #No overlap
                continue
            maxSpan  = max(bed.end-bed.start, read.dBreak - read.uBreak)
            recipOvl = abs(maxStart-minEnd) / float(maxSpan)
            logging.debug("predictVar [%d:%d] - tailRead [%d:%d]" \
                          % (bed.start, bed.end, read.uBreak, read.dBreak))
            if recipOvl < args.overlapbuffer:#not enough overlap
                continue
            
            anno  = read.annotate()    
            reads.append(read.toBriefString())
    
    #ret = ",".join(['t[%s]' % (str(x)) for x in reads])
    return anyCoverage, reads

def spotsSearch_asm(bam, bed, args):
    """
    find spots in high-accuracy contigs

    I'm going to have a problem with Insertion offsets before the variant
    if there are too many of them, I'm going to be effed
    """
    
    #MIS types can't be resolved from spots
    # EXCEPT, however, if they're actually INS/DEL 
    # except a single bp, they could
    if bed.svtype == 'MIS':
        return False, '?',[]
    fetchS = max(0, bed.start - args.fetchbuffer)
    fetchE = min(bed.end + args.fetchbuffer, bam.lengths[bam.references.index(bed.chrom)])
    
    leeway = bed.size * args.sizebuffer
    ref = '?'
    vars = []
    anyCoverage = False
    for read in bam.fetch(bed.chrom, fetchS, fetchE):
        anyCoverage = True
        if read.pos > bed.start or read.aend < bed.end:
            #Not spanning our region
            logging.debug("%s doesn't span region" % (read.qname))
            continue
        ref = False
        if bed.size + leeway > 50000:
            logging.debug("Variant is too long (%dbp), we are assuming reference" % (read.qname, len(read.seq)))
            continue
        
        cigar = spots.expandCigar(read.cigar)
        
        regionStart = max(read.pos, bed.start - args.maxDelta)
        regionEnd = min(read.aend, bed.end + args.maxDelta)
        readPosition = read.pos
        
        c = "".join([str(x) for x in cigar])
        logging.debug(c)
        logging.debug(c[regionStart-read.pos : regionEnd-read.pos])
        if bed.svtype == 'INS':
            match = re.search("(^|[^1])1{%d,%d}([^1]|$)" % (bed.size-leeway, bed.size+leeway), c[regionStart-read.pos : regionEnd-read.pos])
        elif bed.svtype == 'DEL':
            match = re.search("(^|[^2])2{%d,%d}([^2]|$)" % (bed.size-leeway, bed.size+leeway), c[regionStart-read.pos : regionEnd-read.pos])
        
        if match is None:
            ref = True
        else:
            #subtract insertion errors to correct the offset
            #Insertion offset subraction
            subtract = cigar[:(regionStart-read.pos) + match.start()].count(1)
            if bed.svtype == "INS":
                pos = match.start() + read.pos + (regionStart - read.pos) - subtract
                s,e = match.span(); size = e-s
                var = Variant(bed.chrom, pos, pos + 1, bed.svtype, size)
            if bed.svtype == "DEL":
                subtract = cigar[:(regionStart-read.pos)+ match.start()].count(1)
                spos = match.start() + read.pos + (regionStart -read.pos) - subtract
                epos = match.end() + read.pos + (regionStart - read.pos) - subtract
                s,e = match.span(); size = e-s
                var = Variant(bed.chrom, spos, epos, bed.svtype, size)

            vars.append(var)
    
    if len(vars) > 0:
        vars = str(vars[0]) + ("*%d" % len(vars))
    else:
        vars = ""
    
    return anyCoverage, ref, vars
    
def spotsSearch(bam, bed, args):
    """
    take a pysam.Samfile and fetch reads in chrom/start/end region
    see if any reads support the call

    But this doesn't take into account that I have specific groupIds to use...
    """
    leeway = bed.size * args.sizebuffer
    
    fetchS = max(0, bed.start - args.fetchbuffer)
    fetchE = min(bed.end + args.fetchbuffer, bam.lengths[bam.references.index(bed.chrom)])
    if fetchS > fetchE:
        fetchS, fetchE = fetchE, fetchS
    vars = []
    ref = '?'
    anyCoverage = False
    for read in bam.fetch(bed.chrom, fetchS, fetchE):
        if read.pos > bed.start or read.aend < bed.end:
            #Not spanning our region
            logging.debug("%s doesn't span region" % (read.qname))
            continue
        logging.debug("looking at %s" % (read.qname))
        
        anyCoverage = True
        ref = False #we now have the opportunity to find the reference
         
        #I'm going to need md if I get good a MIS
        regionStart = bed.start - args.maxDelta
        regionEnd = bed.end + args.maxDelta
        leeway = bed.size * args.sizebuffer
        
        logging.debug(("regionStart, regionEnd, estSize, leeway, "
                       "estSize+leeway, estSize-leeway"))
        logging.debug("%d %d %d %d %d %d" % (regionStart, regionEnd, \
                        bed.size, leeway, bed.size+leeway, bed.size-leeway))
        foundVar = False
        for svstart, svsize, svtype in spots.expandCigar(read, args.minErr, collapse=3, makeAlt=False):
            if svstart >= regionStart and svstart <= regionEnd and \
                bed.svtype == svtype and \
                bed.size - leeway <= svsize <= bed.size + leeway:
                #check overlap -- I like this logic
                #if (var.start <= bed.start and bed.end <= var.end) \
                     #or (bed.start <= var.start and var.end <= bed.end) \
                     #or (bed.start <= var.end and var.end <= bed.end) \
                     #or (bed.start <= var.start and var.start <= bed.end):
                    #vars.append(var)
                    #foundVar = True
                #else:
                    #maxS = max(var.start, bed.start)
                    #minE = min(var.end, bed.end)
                    #if abs(maxS-minE) <= args.maxDelta:
                        #vars.append(var)
                        #foundVar = True
                foundVar = True
                if svtype == "DEL":
                    vars.append(Variant(bed.chrom, svstart, svstart + svsize, svtype, svsize))
                elif svtype == "INS":
                    vars.append(Variant(bed.chrom, svstart, svstart, svtype, svsize))
                    
        if not foundVar:#this might be broken
            ref = True
                    
    #why only the first and not an average? HOMAlt... would suck
    if len(vars) > 0:
        vars = str(vars[0]) + ("*%d" % len(vars))
    else:
        vars = ""
    
    return anyCoverage, ref, vars

#Here are some helper methods for parsing force annotation results
forceRe = re.compile("(?P<ref>True|False|\?)\[(?P<tails>.*)\|(?P<spots>.*)\]")
def parseForce(data):
    """
    turns the force output into a dict
    """
    if data == 'no_cov' or data == '.':
        return None

    #will fail on malformed entries
    search = forceRe.search(data)
    if search is not None:
        d = search.groupdict()
    else:
        print "problem parsing", data
        return 'prob'

    if d["ref"] == '?':
        d["ref"] = None
    elif d["ref"] == 'True':
        d["ref"] = True
    elif d["ref"] == 'False':
        d["ref"] = False
    d["tails"] = [x for x in d["tails"].split(',') if x != '']
    d["spots"] = [x for x in d["spots"].split(',') if x != '']
    return d

def genoTyper(data):
    """
    """
    if data["ref"] is not None:
        if data["ref"]:
            if len(data["tails"]) > 0 or len(data["spots"]) > 0:
                genoType = "0/1"
            else:
                genoType = "0/0"
            
        elif not data["ref"]:
            if len(data["tails"]) > 0 or len(data["spots"]) > 0:
                genoType = "1/1"
            else:
                genoType = "./."
    else:
        if len(data["tails"]) > 0 or len(data["spots"]) > 0:
            genoType = "./1"
        else:
            genoType = "./."
    
    return genoType

def run(args):
    args = parseArgs(args)
    bam = pysam.Samfile(args.bam)

    #putative caller
    fh = open(args.bed)
    
    tails.BUFFER = args.maxDelta
    
    #CTX and ITX are for breakpoints that have orientations
    #UNK we'll try to find anything that matches (good debugging because we should only be finding
    #support for one of the things we make most always)
    numEntries = 0
    for line in fh.readlines():
        if line.startswith("#"):
            continue
        myentry = line.strip().split('\t')
        
        if not args.bedPE:
            myentry = BedEntry(*myentry)
        else:
            myentry = BedPEEntry(*myentry)
        
        if myentry.svtype not in vtypes:
            #if myentry.svtype == 'UNK':
                #logging.warning("Bed Entry %s is UNK and can't be forced... skipping" % (repr(myentry)))
                #sys.stdout.write(line.strip() + "\t.\n")
                #continue
            #else:
            logging.error("Bed Entry %s svtype column isn't one of %s" % (repr(myentry), str(vtypes)))
            exit(1)
        
        if not args.bedPE:
            if myentry.chrom not in bam.references:
                logging.error("Invalid Chromosome %s" % myentry.chrom)
                continue
            
            mySupport = bedToFakey(myentry)
            tailVars = []
            anyCoverage1 = False
            for i in mySupport:
                ancov1, t = tailsSearch(bam, i, args)
                anyCoverage1 = anyCoverage1 or ancov1
                tailVars.extend(t)
            
            if False:#args.asm:
                anyCoverage2, foundRef, spotVars = spotsSearch_asm(bam, myentry, args)
            else:
                anyCoverage2, foundRef, spotVars = spotsSearch(bam, myentry, args)
        
        else:
            anyCoverage1 = False
            tailVars = []
            mySupport = bedpeToFakey(myentry)
            for i in mySupport:
                ancov1, t = tailsSearch(bam, i, args)
                anyCoverage1 = anyCoverage1 or ancov1
                tailVars.extend(t)
            foundRef = False
            anyCoverage2 = False
            spotVars = ""
        #Eventually, I can get rid of the True/False/? once the genotyper is finished
        #I'm outputing the variant reads and if we found the ref (True, False, ?) where ? means no evidence for or
        # against
        if not anyCoverage1 and not anyCoverage2:
            annot = "no_cov"
            logging.info("no coverage")
        else:
            annot = "%s[%s|%s]" % (foundRef, ",".join(tailVars), spotVars)
            #this is wrong. It's the string's length
            logging.info("Found %d tailed, %d spotted reads" % (len(tailVars), len(spotVars)))
        
        sys.stdout.write(line.strip() + "\t" + annot + "\n")
        numEntries += 1
        if numEntries % 250 == 0:
            sys.stdout.flush()

if __name__ == '__main__':
    run(sys.argv[1:])
