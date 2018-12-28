#!/usr/bin/env python
import re, os, sys, logging, json, copy
from bisect import bisect_left, bisect_right
from optparse import OptionParser
from collections import defaultdict

from pbsuite.utils.setupLogging import *
from pbsuite.utils.FileHandlers import *

from networkx import Graph

USAGE = """Support.py <alignmentFile> <gapInfo> <outFile>

Parses an .m4 or .m5 alignment file and determines which gaps in gapInfo are
supported by reads. As well as Results reported in outFile"""

"""
TODO

Need ability to place scaffolds within captured gaps
    >scaf1
    ====.....====
    >scaf2
    ---

    Actual Genome

    ====.---.====
    That plus the pacbio evidence could polish off that entire gap
"""
#Max Distance a read's alignment can stop from a gap and still be considered 
#for support of the gap
MAXFLANK = 50
#Min number of bases a read needs to reach into a gap to be considered for 
#support of the gap
MINCOVERS = 25

SUPPORTFLAGS = enum(none =    0, \
                    left =    1, \
                    right =   2, \
                    span =    3, \
                    contain = 8  )

ALIGNFLAGS = enum(uniqueMapping = 1, \
                  multiMapping  = 2, \
                  lqpa          = 4, \
                  mostAccurate  = 8, \
                  bestScore     = 16 )

class AlignmentConnector():
    """
    Helper object that will group hits, connect alignments, identify 
    missed adapters, test if a read maps concordantly or discordantly, 
    and test if alignments are within regions
    
    Takes a set of alignments, and parses multi-mapping information
    in order to produce a single alignment.
    """
    
    def __init__(self):
        """
        Only has static methods
        """
        pass
    
    def parseAlignments(self, alignments, minMapq, idAdapters = True):
        """
        Given a set of alignments (from one or more read) put them through the paces
        and return a list of lists of alignments.
        """
        groups = self.groupReadHits(alignments, minMapq)
        
        ret = {}
        for key in groups:
            r = self.connect(groups[key])
            if idAdapters and self.idAdapters(r):
                # if it has adapters, I'll need to 
                # separate the pieces -- idAdapters renames them
                for i in r:
                    logging.debug("About to set qname on %s" % (i))
                    logging.debug(r)
                    try:
                        i.qname = i.qname + "##%d#%d##" % i.trim
                    except TypeError:
                        #didn't get trimmed (happens with 3 subreads per
                        #a single read - only one missed adapter)
                        pass
                    ret[i.qname] = [i]
            else:
                ret[key] = r
        
        return ret.values()
            
    def connect(self, hits, sameStrand=False, sameTar=False):
        """
        Standard procedure to take a group of alignments and 
        connect concordant hits while also using ALIGNFLAGS to characterize
        each read

        Input: a list of hits from a read
        Return: Concordant Reads classified

        Use sameStrand and sameTar(get) to restrict the definition of 
        concordant
        """
        ret = copy.deepcopy(hits)
        con = True
        for h in ret:
            if h.flag != 0:
                con = False
        
        if con:#Don't group compare twice
            self.groupComparison(ret)
        
        ret = self.untangle(ret)
        return ret
        
    def groupComparison(self, hits, flag=True):
        """
        Find the best LQPA, ACC, and HISCORE
        if flag == False:
            we don't change the flag, we just return a list of lists
            [ [bestScore], [mostAccurate], [lqpa] ]
        """
        #The number behind the best; the list of bests (Cause there can be ties)
        LQPA = hits[0].queryPctAligned; lqpas = [hits[0]]
        ACC = hits[0].pctsimilarity; accs = [hits[0]]
        SCORE = hits[0].score; scores = [hits[0]]
        #Find the bests
        for align in hits[1:]:
            if align.queryPctAligned == LQPA:
                lqpas.append(align)
            elif align.queryPctAligned > LQPA:
                LQPA = align.queryPctAligned
                lqpas = [align]
                
            if align.pctsimilarity == ACC:
                accs.append(align)
            elif align.pctsimilarity > ACC:
                ACC = align.pctsimilarity
                accs = [align]
            
            if align.score == SCORE:
                scores.append(align)
            elif align.score < SCORE:
                SCORE = align.score
                scores = [align]
        #Set the apropriate flags
        if not flag:
            return [scores, accs, lqpas]
        for align in lqpas:
            align.flag += ALIGNFLAGS.lqpa 
        for align in accs:
            align.flag += ALIGNFLAGS.mostAccurate
        for align in scores:
            align.flag += ALIGNFLAGS.bestScore
    
    def mappingType(self, hits):
        """
        Sets Reads flag based on Alignment 
        "UniqueMapping", "MultiMapping"
        """
        logging.debug("Determining Read's mapping type from hits")
        fInc = 0
        if len(hits) == 1:
            fInc = ALIGNFLAGS.uniqueMapping
        else:
            fInc = ALIGNFLAGS.multiMapping
        for i in hits:
            i.flag += fInc

    def isQueryConcordant(self, A, B, sameStrand=False):
        """
        Takes two reads (A,B) and checks their mapping concordancy
        meaning the query alignment positions are reasonable 
        given a sequencing read
        
        for example:
            concordant
                q ===---   ---===
                t ===============
            !concordant
                q ---===   ---===
                t ===============

        sameStrand == True:
            read's hits need to be on the sameStrand

        Return  0 for False
        Return  1 for A is upstream (5' - 3') of B
        Return -1 for B is upstream 
        """
        if sameStrand and A.tstrand != B.tstrand:
            return 0
        
        if A.qstart < A.qend < B.qstart < B.qend:#
            return 1
        if B.qstart < B.qend < A.qstart < A.qend:
            return -1
        
        return 0
        
    def isTargetConcordant(self, A, B, sameStrand=False, sameTar=False):
        """
        Takes two alignments and checks if they're concordant. 
        if sameStrand:
            hits must be on the same strand
        if sameTar:
            hits must be on the same target

        Return  0 for False
        Return  1 for A is upstream (5' - 3') of B
        Return -1 for B is upstream 

        Question:
            Can you have hits that are target concordant but not 
            query concordant? 

            queryA ==-
            queryB  -==
            target ----
            A deletion?
        """
        if sameStrand and a.tstrand != B.tstrand:
            return 0
            
        if sameTar and A.tname != B.tname:
            return 0
            
        if A.tend < B.tstart and A.qend < B.qstart:
            return 1
        if B.tend < A.tstart and B.qend < A.qstart:
            return -1
        
        return 0
    
    def extendsTarget(self, alignment, maxFlank=50, minCovers=25):
        """
        Checks to see if a read extends it's target
        Returns the direction to which the target is 
        extended using SUPPORTFLAGS
        """
        ret = SUPPORTFLAGS.none
        logging.debug("maxFlank %d - minCovers %d" % (maxFlank, minCovers))
        logging.debug("Checking 5End of Scaff %s %s" % (alignment.tname, alignment.qname))
        ret += self.supportsRegion(alignment, alignment.tname, \
                                         -sys.maxint, 0, maxFlank, minCovers)
        logging.debug("Checking 3End of Scaff %s %s" % (alignment.tname, alignment.qname))
        ret += self.supportsRegion(alignment, alignment.tname, \
                                      alignment.tseqlength, sys.maxint, maxFlank, minCovers)
        
        #Orientation correct - how we're extending target
        if ret == SUPPORTFLAGS.right:
            ret = SUPPORTFLAGS.left
        elif ret == SUPPORTFLAGS.left:
            ret = SUPPORTFLAGS.right
        
        return ret

    def supportsRegion(self, alignment, rName, rStart, rEnd, maxFlank=50, minCovers=25):
        """
        Checks if read's alignment supports a region. 
        if r(egion)Name/Start/End is not specified the target sequence's boundaries are used
        
        The region of support can be extended using:
            maxFlank  = Max distance from end query alignment can stop to still be a candiate
            minCovers = Min amount of sequence that extends the target
        
        Uses the SUPPORTFLAGS eunum for return values
        Returns "left" for covering the region's left most point
        Returns "right" for covering the region's right most point
        Returns "span" for covering left and right
        Returns "contain" for being within left and right point
        Returns "none" if no relation (different targets - contained - otherwise)
        """

        ret = SUPPORTFLAGS.none
        
        if rName != alignment.tname:
            return ret
        
        if alignment.tstrand == "0":
            #Moving into Region from left
            #Meaning we extend the region into the left (to the right)
            distanceFromEnd = rStart - alignment.tend
            remainingReadSeq3 = (alignment.qseqlength - alignment.qend) - minCovers
            logging.debug("+ Strand on " + alignment.qname)
            logging.debug("LeftDist %d remainSeq %d" % (distanceFromEnd, remainingReadSeq3))
            if distanceFromEnd >= 0 and \
               distanceFromEnd < remainingReadSeq3 and \
               distanceFromEnd <= maxFlank :
                #Positive Strand Maps on Left Contig and enters gap     
                ret += SUPPORTFLAGS.left
                logging.debug("left")
            
            #moving out of Region to right
            #meaning we extend into the right (to the left)
            distanceFromBeginning = alignment.tstart - rEnd
            remainingReadSeq5 = alignment.qstart - minCovers
            logging.debug("RightDist %d remainSeq %d" % (distanceFromBeginning,remainingReadSeq5))
            if distanceFromBeginning >= 0 and \
               distanceFromBeginning < remainingReadSeq5 and \
               distanceFromBeginning <= maxFlank :
                #Positive Strand Maps on Right Contig and Exits Gap
                ret += SUPPORTFLAGS.right
                logging.debug("right support")
                 
            if alignment.tstart <= (rStart - maxFlank) and alignment.tend >= (rEnd + maxFlank):
                ret = SUPPORTFLAGS.span
                logging.debug("span support")
            
            elif alignment.tstart >= rStart and alignment.tend <= rEnd:
                ret = SUPPORTFLAGS.contain
                logging.debug("contain support")

        elif alignment.tstrand == "1":
            #Moving into region from left  on - strand
            #Meaning we extend to the right on + strand 
            distanceFromBeginning = alignment.tstart - rEnd
            remainingReadSeq3 = (alignment.qseqlength - alignment.qend) - minCovers
            
            logging.debug("- Strand on "+alignment.qname)
            logging.debug("RightDist %d remainSeq %d" % (distanceFromBeginning,remainingReadSeq3))
            if distanceFromBeginning >= 0 and \
               distanceFromBeginning < remainingReadSeq3 and \
               distanceFromBeginning <= maxFlank :
                ret += SUPPORTFLAGS.right
                logging.debug("right support")
            
            #Moving out of region to the right on - strand
            #Meaning we extend to the left on + strand
            distanceFromEnd = rStart - alignment.tend
            remainingReadSeq5 = alignment.qstart - minCovers
            logging.debug("LeftDist %d remainSeq %d" % (distanceFromEnd,remainingReadSeq5))
            if distanceFromEnd >= 0 and \
               distanceFromEnd < remainingReadSeq5 and \
               distanceFromEnd <= maxFlank :
                ret += SUPPORTFLAGS.left
                logging.debug("left support")
                
            if alignment.tstart <= rStart - minCovers and alignment.tend >= rEnd + minCovers:
                ret = SUPPORTFLAGS.span
                logging.debug("span support")
                
            elif alignment.tstart >= rStart and alignment.tend <= rEnd:
                ret = SUPPORTFLAGS.contain
                logging.debug("contain support")
        logging.debug("")
        return ret
        
    def groupReadHits(self, alignments, minMapq):
        logging.debug("Grouping Read Hits")
        reads = defaultdict(list)#readname: [hit hit hit]
        
        for line in alignments:
            if line.mapqv >= minMapq:
                reads[line.qname].append(line)
            else:
                logging.debug("Hit for %s has mapq %d - below threshold %d" % (line.qname, line.mapqv, minMapq))
        
        return reads    
    
    def getBestScore(self, reads):
        """
        Gets the read with the best score.
        parses all of the reads and returns the best hit of the group
        
        best hit is the alignment with the best score
        and if there is a tie, the tie break is based
        on alignment score
        
        and the to just take the first of the remaining alignments
        """
        if len(reads) == 1:
            return reads[0]
        #Get who we're using
        bestScore = []
        mostAccurate = []
        for read in reads:
            if read.flag & ALIGNFLAGS.bestScore:
                bestScore.append(read)
            if read.flag & ALIGNFLAGS.mostAccurate:
                mostAccurate.append(read)
        
        if len(bestScore) == 0:
            #Actually have to hunt it down, we're looking at some weiners
            bestScore, mostAccurate, lqpa = self.groupComparison(reads, False)
            #return None
        
        #Tie breaking
        anchor = None
        if len(bestScore) != 1:
            if len(mostAccurate) == 1:
                for i in bestScore:
                    if i == mostAccurate[0]:
                        anchor = mostAccurate[0]
        
        #Just arbiturarily take one
        if anchor is None:
            anchor = bestScore[0]
        
        return anchor

    def untangle(self, reads):
        """
        Given a group of subread's hits, see if we can eliminate the 
        spurious repeat matches
        
        Find LQPA 
        recursively call(give me all the concordant neighbors to the left and to the right)
        continue until there isn't exactly one neighbor
        """
        if len(reads) < 2:
            return reads
        
        anchor = self.getBestScore(reads)
        #I don't like this happening.. something is strange
        if anchor == None:
            logging.warning("Read %s doesn't have a best hit" % (reads[0].qname))
            return []
        
        newReads = self.layout(anchor, list(reads))
        
        newReads.sort(cmp=lambda x,y: x.qstart < y.qstart)
        return newReads
    
    """
    after I have my anchor, I want to find queryConcordant reads.
    I'm going to find all the reads that are queryConcordant outside of the anchor.
    If there is more than one read per side, I'll tie break those by 
     getting the one that is considered the best hit
        
    1 - Anchor solution is above.
    2 - layout 1
    3 - layout -1
    """
    def layout(self, anchor, allReads):
        """
        uses the anchor to find all unique concordant hits from the anchor 
        side = -1, check UpStream(5') of anchor
        side = 1, check DownStream (3') of anchor
        """
        foundNew = True
        #holds the min/max used qstart qend
        regionHolder = copy.deepcopy(anchor)
        ret = [anchor]
        while foundNew:
            foundNew = False
            reads = filter(lambda x: x not in ret, allReads)
            if len(reads) == 0:
                continue
            sides = { 1: [],
                     -1: []}
            
            for read in reads:
                side = self.isQueryConcordant(regionHolder, read)
                if side != 0:
                    foundNew = True
                    sides[side].append(read)
            
            if len(sides[1]) == 1:
                regionHolder.qstart = sides[1][0].qstart
                #region update (upstream)
                ret.append(sides[1][0])
            elif len(sides[1]) > 1:#Let all these guys duke it out
                bestHit = self.getBestScore(sides[1])
                if bestHit is None:
                    logging.warning("Here")
                    logging.warning("\n"+"\n".join([str(x) for x in sides[1]]))
                newGuys = self.layout(bestHit, sides[1])
                regionHolder.qstart = newGuys[0].qstart
                ret.extend(newGuys)
                #region update (upstream)
                
            if len(sides[-1]) == 1:
                bestHit = sides[-1]
                #region update (downstream)
                regionHolder.qend = sides[-1][0].qend
                ret.append(sides[-1][0])
            elif len(sides[-1]) > 1:
                bestHit = self.getBestScore(sides[-1])
                newGuys = self.layout(bestHit, sides[-1])
                regionHolder.qend = newGuys[-1].qend
                ret.extend(newGuys)
        
        ret.sort(cmp=lambda x,y: x.qstart < y.qstart)
        return ret

    def isDiscordant(self, alignment, tailAllowed=150):
        """
        If a read maps with long tails that don't map off the target's ends, call it discordant
        
        tailAllowed = max tail length allowed to not map to the reference
        
        Note - This doesn't check a reaason for the discordantcy (sp) -
             - Only checks for long tails.
             
        Deprecated again
        """
        orientation = self.extendsTarget(alignment)
        if orientation == SUPPORTFLAGS.span:
            #Not much to prove here.
            return False
            
        isDiscord = True
        if alignment.tstrand == "0":
            threeLen = alignment.qseqlength - alignment.qend
            fiveLen = alignment.qstart
            
            if orientation == SUPPORTFLAGS.left and fiveLen <= tailAllowed:
                isDiscord = False
            elif orientation == SUPPORTFLAGS.right and threeLen <= tailAllowed:
                isDiscord = False
            
        elif alignment.tstrand == "1":
            
            threeLen = alignment.qstart
            fiveLen = alignment.qseqlength - alignment.qend
            
            if orientation == SUPPORTFLAGS.left and fiveLen <= tailAllowed:
                isDiscord = False
            elif orientation == SUPPORTFLAGS.right and threeLen <= tailAllowed:
                isDiscord = False
            
        return isDiscord
    
    def idAdapters(self, reads):
        """
        Using the alignments for a particular read,
        Identify adapters and return a new list of reads that have been split up
        I'll rename the reads that have adapters
        
        ##split those subreads into two alignments
        Consider a subread that overlaps with itself as a missed adapter, split it.
            example - for a single read
            hit1:
                qstart 0: qend 100: tstart = 400: tend=500: strand = 0
            hit2:
                qstart = 150: qend = 250: tstart = 400: tend = 500: strand = 1
        this is highly indicitive of a missed adapter
        
        reads is a list of alignments that must be sorted in query concordant order 
            (meaning the hits are in order of how they should have been sequenced...
             return of AlignmentConnector.untangle makes this sort)
        
        """
        split = False
        for i in range(1, len(reads)):
            a = reads[i-1]
            b = reads[i]
            
            #Major chacteristics of missed adapter
            #Same target - different strand
            #Aligns in the same area
            if a.tname == b.tname and a.tstrand != b.tstrand \
               and (abs(a.tstart - b.tstart) <= 75 or abs(a.tend - b.tend) <= 75)\
               and (b.qstart - a.qend) <= 50:
                
                #Do some more tests to make sure that they overlap in the
                #Way we expect missed addapters to overlap
                ovl = self.supportsRegion(a, b.tname, b.tstart, b.tend)
                if a.tstrand == '0' and (ovl & SUPPORTFLAGS.right or \
                                         ovl & SUPPORTFLAGS.contain or\
                                         ovl & SUPPORTFLAGS.span):
                    split = True
                    logging.debug("Missed adapter in read! %s" % a.qname)
                    #We have a missed adapter
                    a.trim = (0, a.qend)
                    a.qseqlength = a.qend
                    #a.qname += "##%d#%d##" % a.trim

                    b.trim = (b.qstart, b.qseqlength)
                    shift = b.qend - b.qstart
                    b.qstart = 0
                    b.qend -= shift
                    b.qseqlength -= shift
                    #b.qname += "##%d#%d##" % b.trim

                elif a.tstrand == '1' and (ovl & SUPPORTFLAGS.left or \
                                           ovl & SUPPORTFLAGS.contain or\
                                           ovl & SUPPORTFLAGS.span):
                    split = True
                    logging.debug("Missed adapter in read! %s" % a.qname)
                    b.trim = (0, b.qend)
                    b.qseqlength = b.qend
                            
                    a.trim = (a.qstart, a.qseqlength)
                    shift = a.qend - a.qstart
                    a.qstart = 0
                    a.qend -= shift
                    a.qseqlength -= shift
            
        return split
        
class GapSupporter():
    """
    Holds a gapInfo file, as you feed self.support reads,
    it'll figure out the support and then keep track of it.
    
    you can add an existing alignCon to the GapSupporter if you've
    made one elsewhere and want to save space. But another will
    be automatically made if you don't

    Also, this does the bookkeeping to track what gaps are supported
    by what reads and how.
    
    Can run summary statistics on how many gaps are supported and such
    """
    def __init__(self,  gapInfo, alignCon=None):
        self.gapInfo = gapInfo
        
        try:
            self.gapIndex = self.gapInfo.getSortedGaps()
        except Exception:#For any reason
            self.gapIndex = None
        
        self.alignCon = AlignmentConnector() if alignCon is None else alignCon
        self.gapGraph= GapGraph()
        
        #{ readName: [(alignment, [nodeName, nodeName]), ] }
        #Same as that weird name manipulation as before.
        #see if it is a span
        #self.readSupport = {}
        
    def classifyRead(self, alignmentGroup, capturedOnly=False):
        """
        For each alignment in a read, add information about read supporting
        gaps
        alignmentGroup must be in order

        spanOnly == True:
            only add span support
        capturedOnly == False
            ignore scaffold gaps
        """
        self.capturedGapSupport(alignmentGroup)
        if not capturedOnly:
            self.scaffoldGapSupport(alignmentGroup)
            #self.consolidate spport
     
    def capturedGapSupport(self, alignmentGroup):
        """
        returns a list of nodes in the graph the read supports
        """
        # gapName : SUPPORTFLAG
        ret = defaultdict(int)
        #Finding gaps to support and classifying
        for alignment in alignmentGroup:
            candidates = []
            scaffold = alignment.tname.split('|')[-1]
            if self.gapIndex is not None:
                try:
                    #Find range of gaps we can potentially support 
                    #I play it safe and get up to 2 extra gaps we could support
                    startIndex = max(0, bisect_left(self.gapIndex[scaffold], alignment.tstart) - 1)
                    endIndex = bisect_right(self.gapIndex[scaffold], alignment.tend) + 1
                    candidates = self.gapIndex[scaffold][startIndex:endIndex]
                except KeyError:
                    pass#No gaps
            else:#never happens?
                gaps = filter(lambda x: x.startswith(scaffold), self.gapInfo.keys())
                for key in gaps:
                    candidates.append(self.gapInfo[key])
            
            for gap in candidates:
                logging.debug("gapSup")
                supType = self.alignCon.supportsRegion(alignment,  \
                                            rName=gap.scaffold, \
                                            rStart=gap.start,\
                                            rEnd = gap.end)
                if supType == SUPPORTFLAGS.none:
                    continue
                ret[gap.name] += supType
        #Adding support to graph
        readName = alignmentGroup[0].qname
        for gapName in ret:
            #New way
            gap = self.gapInfo[gapName]
            supType = ret[gapName]
            lftNode = gap.leftContig + "e3"
            rhtNode = gap.rightContig + "e5"
            
            #Extends gap to left or rightContig to right
            if supType == SUPPORTFLAGS.left:
                self.gapGraph.add_extend(lftNode, readName)
            elif supType == SUPPORTFLAGS.right:
                self.gapGraph.add_extend(rhtNode, readName)
            elif supType == SUPPORTFLAGS.span:
                self.gapGraph.add_evidence(lftNode, rhtNode, readName)
        
        return ret

    def scaffoldGapSupport(self, alignmentGroup):
        """ 
        Checks to see if and how an alignment supports between scaffolding gaps
        
        alignmentGroup must be sorted by query portion used
        sourceA is alignmentGroup[0] support
        sourceB is alignmentGroup[1] support
        add a link between A and B
        sourceC is alignmentGroup[2] support
        """
        readName = alignmentGroup[0].qname
        logging.debug("looking at: " + readName + " for scaffold extension/unification")
        
        alignmentGroup.sort(cmp=lambda x,y: x.qstart - y.qstart)
        anchor = self.alignCon.getBestScore(alignmentGroup)
        
        flags = []
        
        logging.debug("Building flags table")
        logging.debug("%d - %s" % (len(alignmentGroup), " ".join(map(str, alignmentGroup))))
        for alignment in alignmentGroup:
            #[(Flag, lftNode, rhtNode),..]
            base = alignment.tname.split('|')[-1]
            flag = self.alignCon.extendsTarget(alignment)
            flags.append((flag,                 #SupFlag
                          base+"e5",            #LftNode
                          base+"e3",            #RhtNode
                          alignment.tstrand))   #Strand   
            logging.debug(str(base) + " " + str(flag) + " " + str(alignment))
        #Solo..
        if len(alignmentGroup) == 1:
            logging.debug("Solo alignment comparison")
            flag1, lftNode1, rhtNode1, strand1 = flags[0]
            if flag1 == SUPPORTFLAGS.span:
                #This can't happen...
                self.gapGraph.add_extend(lftNode1, readName)
                self.gapGraph.add_extend(rhtNode1, readName)
            elif flag1 == SUPPORTFLAGS.right:
                self.gapGraph.add_extend(rhtNode1, readName)
            elif flag1 == SUPPORTFLAGS.left:
                self.gapGraph.add_extend(lftNode1, readName)
            return
        
        index = alignmentGroup.index(anchor)
        
        #go upstream of anchor
        pairs = []
        for i in range(index, 0, -1):
            #pairs.append((alignmentGroup[i], flags[i], alignmentGroup[i-1], flags[i-1]))
            pairs.append((alignmentGroup[i-1], flags[i-1], alignmentGroup[i], flags[i]))
            logging.debug("made pair upstream")
            logging.debug(" ".join(map(str,[alignmentGroup[i], flags[i], alignmentGroup[i-1], flags[i-1]])))
        self.__scaffRangeDo__(pairs)
        #Go downstream of anchor
        pairs = []
        for i in range(index, len(alignmentGroup) - 1, 1):
            pairs.append((alignmentGroup[i], flags[i], alignmentGroup[i+1], flags[i+1]))
            logging.debug("made pair downstream")
            logging.debug(" ".join(map(str,[alignmentGroup[i], flags[i], alignmentGroup[i+1], flags[i+1]])))
        self.__scaffRangeDo__(pairs)
        
    
    def __scaffRangeDo__(self, pairs):
        """
        Go through the pairs and see if there is a link between the them
        """
        if len(pairs) == 0:
            return
        readName = pairs[0][0].qname
        i = 1
        for align1, flags1, align2, flags2 in pairs:
            logging.debug(flags1)
            logging.debug(flags2)
            flag1, lftNode1, rhtNode1, strand1 = flags1
            flag2, lftNode2, rhtNode2, strand2 = flags2
            #Super Conservative, I'm not going to let
            #ends of multi-mappers extend scaffolds
            if strand1 == '0':
                if flag1 == SUPPORTFLAGS.right or flag1 == SUPPORTFLAGS.span:#extend correctly
                    if strand2 == '0':
                        if flag2 == SUPPORTFLAGS.left or flag2 == SUPPORTFLAGS.span:
                            self.gapGraph.add_evidence(rhtNode1, lftNode2, readName)
                    elif strand2 == '1':
                        if flag2 == SUPPORTFLAGS.right or flag2 == SUPPORTFLAGS.span:
                            self.gapGraph.add_evidence(rhtNode1, rhtNode2, readName)
                    else:
                        return; #We've broken from the anchor
                else:
                    return; #we've broken the anchor chain
            
            elif strand1 == '1':
                #I never strand corrected? left means right, right means left
                if flag1 == SUPPORTFLAGS.left or flag1 == SUPPORTFLAGS.span:#extend correctly
                    if strand2 == '0':
                        if flag2 == SUPPORTFLAGS.left or flag2 == SUPPORTFLAGS.span:
                            self.gapGraph.add_evidence(lftNode1, lftNode2, readName)
                    elif strand2 == '1':
                        if flag2 == SUPPORTFLAGS.right or flag2 == SUPPORTFLAGS.span:
                            self.gapGraph.add_evidence(lftNode1, rhtNode2, readName)
                    else:
                        return; #We've broken from the anchor           
                else:
                    return; #we've broken the anchor chain

class Support():
    """
    Worker for this script
    Takes Reads, Connects the Alignments, Classifies the Gap Support
    """
    
    def __init__(self):
        self.parseArgs()
        setupLogging(self.options.debug)
    
    def parseArgs(self):
        parser = OptionParser(USAGE)
        
        parser.add_option("-m", "--minMapq", default=200, type=int, \
                          help=("Minimum MapQ of a read to be considered "
                                "for support (200)"))
        parser.add_option("--spanOnly", action="store_true", \
                          help=("Only allow support by reads that span an"
                                " entire gap. i.e. no contig extension."))
        parser.add_option("--capturedOnly", action="store_true", \
                  help=("Only find support for captured gaps. "\
                        " i.e. no between-scaffold gap-filling"))
        parser.add_option("--debug", action="store_true", \
                          help="Increases verbosity of logging" )

        self.options, args = parser.parse_args()
        
        if len(args) != 3:
            parser.error("Error! Incorrect number of arguments")
        
        if not os.path.isfile(args[0]):
            parser.error("Error! Alignment File Does Not Exist")
        self.alignmentFileName = args[0]
        
        if not os.path.isfile(args[1]):
            parser.error("Error! Gap Info File Does Not Exist") 
        self.gapFileName = args[1]
        
        if os.path.isfile(args[2]):
            sys.stderr.write("[WARNING] Output File Being Overwritten!")
        self.outputFileName = args[2]
        
        self.gapInfo = GapInfoFile(self.gapFileName)
        if os.path.splitext(self.alignmentFileName)[1] == '.m4':
            self.alignments = M4File(self.alignmentFileName)
        elif os.path.splitext(self.alignmentFileName)[1] == '.m5':
            self.alignments = M5File(self.alignmentFileName)
        else:
            parser.error("Error! Alignment File Extension (%s) not recognized." \
                         % os.path.splitext(self.alignmentFileName)[1])
    
    def run(self):
        """
        Given an alignment file, put it through the paces that will figure out what
        gaps each read supports
        """
        logging.info("Building Helper Objects")
        connector = AlignmentConnector()
        supporter = GapSupporter(self.gapInfo, alignCon = connector)
        logging.info("Connecting Alignments")
        alignments = connector.parseAlignments(self.alignments, self.options.minMapq)
        
        logging.info("Classifying Alignments' Support")
        for readGroup in alignments:
            supporter.classifyRead(readGroup, self.options.capturedOnly)
        
        logging.info("Saving Support Graph")
        supporter.gapGraph.saveGraph(self.outputFileName, self.options.spanOnly)
        logging.info("Finished")

if __name__ == '__main__':
    main = Support()
    main.run()
    pass
