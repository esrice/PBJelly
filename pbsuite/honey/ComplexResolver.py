#!/usr/bin/env python
import itertools, copy, argparse
from collections import defaultdict, Counter
from intervaltree import bio
from pbsuite.utils.setupLogging import *
import networkx as nx

USAGE = ("""Parses a hon.tails file to find overlapping tails that may "
indicate a complex (multi-breakpoint) event and tries to resolve them""")

class BreakPoints():
    def __init__(self,id, start, end, dir1, dir2, data):
        self.id = id
        self.start = start
        self.end = end
        self.dir1 = dir1
        self.dir2 = dir2
        self.data = data
        
    def __eq__(self, other):
        return self.id == other.id
        if self.start == other.start and self.end == other.end and \
           self.dir1 == other.dir1 and self.dir2 == other.dir2:
            return True
        elif self.start == other.end and self.end == other.end and \
           self.dir1 == other.dir2 and self.dir2 == other.dir1:
            return True
        return False
    
    def __hash__(self):
        return hash(self.id)
        j = [self.start, self.end, self.dir1, self.dir2]
        j.sort()
        return hash(str(j))
    
    def __str__(self):
        return "id%s s%d e%d 1%s 2%s " % (self.id, self.start, self.end, self.dir1, self.dir2)
    
    def __repr__(self):
        return "BreakPoint(%s)" % (str(self))

class Block():
    """
    Holds a block of genomic sequence
    """
    def __init__(self, id, start, end, strand):
        self.id = id
        self.start = start
        self.end = end
        self.strand = strand
    def __str__(self):
        return "%s %r %r %s" % (self.id, self.start, self.end, self.strand)
        
    def __eq__(self, other):
        return self.id == other.id
    
    def __ne__(self, other):
        return self.id != other.id

    def __cmp__(self, other):
        if self == other:
            return 0
        if self.id < other.id:
            return -1
        if self.id > other.id:
            return 1
        
    def __hash__(self):
        return hash(self.id)

class ComplexResolver():
    
    def __init__(self, fileName, minBlock=500, maxSpan=100000, maxOvl=10, maxRefBlocks=10, outFile=None):
        self.fileName = fileName
        self.minBlock = minBlock
        self.maxSpan = maxSpan
        self.maxOvl = maxOvl
        self.maxRefBlocks = maxRefBlocks
        if outFile is None:
            self.outFile = sys.stdout
        else:
            self.outFile = open(outFile, 'w')
        
    def run(self):
        self.readTails()
        for chrom in self.points:
            #cluster tails -- need the points that overlap
            for cluster in self.makeClusters(chrom):
                #Need to compress the points here
                refBlocks, gInterval, blockLookup = self.makeReferenceBlocks(cluster)
                
                if refBlocks is None:
                    continue
                
                self.outFile.write("RefBlocks: %s\n" % chrom)
                for i in refBlocks:
                    self.outFile.write("\t" + str(i) + "\n")
                
                self.outFile.write(" -> ".join([x.id for x in refBlocks]) + '\n')
                hasValid  = False
                for hypo in self.makeHypothesis(cluster):
                    valid, sampleBlocks = self.testHypothesis(hypo, gInterval, blockLookup)
                    if valid:
                        hasValid = True
                        self.outFile.write("Events:\n")
                        self.outFile.write("\t" + "\n\t".join(self.annotate(refBlocks, sampleBlocks)) + "\n")
                        self.outFile.write("SampleBlocks:\n")
                        self.outFile.write("\t" + "\n\t".join([str(x) for x in sampleBlocks]) + '\n')
                        d = []
                        for x in sampleBlocks:
                            if x.strand == '-':
                                x.id = x.id.lower()
                        self.outFile.write(" -> ".join([x.id for x in sampleBlocks]) + '\n')
                if not hasValid:
                    self.outFile.write("Events:\nSampleBlocks:\n")
                    self.outFile.write("No Resolution Found\n")
                #Need to report the events
                self.outFile.write(str("#"*10) + '\n')
    
    def makeClusters(self, chrom):
        """
        yield all non-solo clusters
        returns all of the BreakPoints that comprise the events within the
        cluster
        """
        #PARAM
        MAXSPAN = self.maxSpan
        MAXOVLS = self.maxOvl
        
        intervalTree = self.genome[chrom]
        overlaps = []
        
        graph = nx.Graph()
        #Find clusters on the chromosome
        #By building a graph
        for interval in intervalTree:
            if abs(interval.end - interval.begin) > MAXSPAN:
                continue
            ovls = intervalTree.search(interval)
            if len(ovls) == 1 or len(ovls) > MAXOVLS:
                continue
            for a,b in itertools.combinations(ovls, 2):
                if abs(a.end - a.begin) > MAXSPAN \
                   or abs(b.end - b.begin) > MAXSPAN:
                    continue
                graph.add_edge(a, b)
            
        #clusters are sub graphs
        for subG in nx.connected_component_subgraphs(graph):
            if len(subG.edges()) == 0:
                continue
            
            lookup = [x.data for x in subG.nodes()]
            ret = [x for x in self.points[chrom] if x.id in lookup]
            
            #I need to merge points that are too close
            
            if len(ret) <= 10:
                yield ret

    def readTails(self):
        """
        Opens up a .tails file and returns a list of all the redundant points
        and your genome interval tree
        """
        fh = open(self.fileName, 'r')
        
        self.points = defaultdict(list)
        self.genome = bio.GenomeIntervalTree()
        
        for line in fh.readlines():
            if line.startswith("#"):
                continue
            
            #Need a tails parser
            data = line.strip().split('\t')
            id = data[0]
            chrom = data[2]
            start = int(data[3])
            end = int(data[6])
            svtype = data[9]
            
            directions = [tuple([x[:3], x[-3:]]) for x in data[12].split(';')]
            d1s = [x[0] for x in directions]
            d2s = [x[1] for x in directions]
            
            #I make two, one for using the first bp first and one for second bp first
            self.points[chrom].append(BreakPoints(id, start, end, d1s, d2s, data))
            self.points[chrom].append(BreakPoints(id, end, start , d2s, d1s, data))
            
            self.genome[chrom].addi(start, end, data=id)
        
    def makeReferenceBlocks(self, cluster):
        """
        #Make my blocks -- These are the pieces of the reference
        """
        #PARAM max num region blocks
        MAXREFBLOCKS = self.maxRefBlocks
        
        blockLookup = {}
        refBlocks = []
        #Need to reconstruct the genome... I know this is wasteful 
        #but I've fucked up somewhere
        region = bio.IntervalTree()
        for i in cluster:
            region.addi(i.start, i.end)
        
        #Give the source and sink regions
        region.addi(region.begin() - (self.minBlock + 1), \
                    region.end() +   (self.minBlock + 1) )
        region.split_overlaps()
        
        #remove blocks that are too small
        for r in list(region):
            size = abs(r.end - r.begin)
            if size > 0 and size <= self.minBlock:
                region.remove(r)
        
        i = ord('A')
        for interval in sorted(region):
            if interval.begin == region.begin():
                start = 'source'
            else:
                start = interval.begin
            
            if interval.end == region.end():
                end = 'sink'
            else:
                end = interval.end
            #need a blocklookp for finding from keys .. why?
            #and I need the refBlock list
            blockLookup[interval] = Block(chr(i), start, end, '+')
            refBlocks.append(Block(chr(i), start, end, '+'))
            i += 1
            if i - ord('A') > MAXREFBLOCKS:
                return None, None, None
        
        #logging.debug("ReferenceBlockStructure " +  " -> ".join([str(x) for x in refBlocks]))
        return refBlocks, region, blockLookup
            
    def makeHypothesis(self, myBps):
        """
        This is a generator that will yield permutations of breakpoints 
        for testing
        """
        #need to create the non-redundant set of every permutation
        nonRedun = set()
        for combo in itertools.combinations(myBps, len(myBps)/2):
            if len(combo) == len(set(combo)):
                nonRedun.add(combo)
            
        #for each set in the nonRedund set, I'll need to get the permutations
        for i in nonRedun:
            for h in itertools.permutations(i):
                yield h
        
    def __fetchBlocks(self, start, end, gInterval, blockLookup, strand='+'):
        """
        Returns the blocks betwen start and end on chromosome
        ... might need some strandedness..
        """
        a, b = (start, end) if start < end else (end, start)
        ret = []
        for interval in sorted(gInterval.search(a,b)):
            b = copy.deepcopy(blockLookup[interval])
            #b.start = start
            #b.end = end
            b.strand = strand
            ret.append(b)
        return ret

        
    def testHypothesis(self, hypo, gInterval, blockLookup):
        """
        Given a hypothesis of how to move through the SV, test if it's possible
        hypothesis need chromosome information
        Returns: 
        sampleBlocks -- The potential structure of the Sample
        """
        #I need to include the sets of not using everything because I can move through
        #FPs
        transPnt = 0
        transDir = '->'
        moveThroughCnt = 0
        output = ["Hypothesis Result:"]
        sampleBlocks = []
        
        for breakpoint in hypo:
            if transDir == '->' and transPnt < breakpoint.start  and \
            len([x for x in breakpoint.dir1 if x.startswith(transDir)]) != 0:
                output.append((transPnt, transDir, breakpoint.start, "JUMP", breakpoint.id, breakpoint.end))
                
                sampleBlocks.extend(self.__fetchBlocks(transPnt, breakpoint.start, gInterval, blockLookup))
                
                odirs = [(x,y) for x,y in zip(breakpoint.dir1, breakpoint.dir2) if x.startswith(transDir)]
                transDir = odirs[0][1].replace('e','').replace('i','').replace('p','')
                transPnt = breakpoint.end
                
                moveThroughCnt += 1
            
            elif transDir == '->' and transPnt < breakpoint.end and \
                len([x for x in breakpoint.dir2 if x.startswith(transDir)]) != 0:
                
                output.append((transPnt, transDir, breakpoint.end, "JUMP", breakpoint.id, breakpoint.start))
                
                sampleBlocks.extend(self.__fetchBlocks(transPnt, breakpoint.end, gInterval, blockLookup))
                
                odirs = [(x,y) for x,y in zip(breakpoint.dir1, breakpoint.dir2) if y.startswith(transDir)]
                transDir = odirs[0][0].replace('e','').replace('i','').replace('p','')
                transPnt = breakpoint.start
                
                moveThroughCnt += 1
                
            elif transDir == '<-' and  transPnt > breakpoint.start and \
                len([x for x in breakpoint.dir1 if x.endswith(transDir)]) != 0:
                output.append((transPnt, transDir, breakpoint.start, "JUMP", breakpoint.id, breakpoint.end))
                
                sampleBlocks.extend(self.__fetchBlocks(transPnt, breakpoint.start, gInterval, blockLookup, strand='-'))
                
                odirs = [(x,y) for x,y in zip(breakpoint.dir1, breakpoint.dir2) if x.endswith(transDir)]
                transDir = odirs[0][1].replace('e','').replace('i','').replace('p','')
                transPnt = breakpoint.end
                
                moveThroughCnt += 1
                
            elif transDir == '<-' and transPnt > breakpoint.end and \
                len([x for x in breakpoint.dir2 if x.endswith(transDir)]) != 0:
                output.append((transPnt, transDir, breakpoint.end, "JUMP", breakpoint.id, breakpoint.start))
                    
                sampleBlocks.extend(self.__fetchBlocks(transPnt, breakpoint.end, gInterval, blockLookup, strand='-'))
                
                odirs = [(x,y) for x,y in zip(breakpoint.dir1, breakpoint.dir2) if y.endswith(transDir)]
                transDir = odirs[0][1].replace('e','').replace('i','').replace('p','')
                transPnt = breakpoint.end
                
                moveThroughCnt += 1
            else:
                #We can't get into the next breakpoint, Fail
                #I don't know why I'm scared of this
                return False, None
        
        output.append((transPnt, transDir, 'inf'))
        
        sampleBlocks.extend(self.__fetchBlocks(transPnt, gInterval.end(), gInterval, blockLookup))
        
        #Do we use all the breakpoints and ENTER(A) and EXIT(A+nblocks-1) the SV
        #this is too strict, though
        #I don't think I need this because I have the False, None short circuit
        #AND I'm only testing a single hypothesis
        #if moveThroughCnt == len(hypo)/2 and sampleBlocks[0]=='A' \
        if len(sampleBlocks) > 1 and sampleBlocks[0].id == 'A' \
           and sampleBlocks[-1].id == chr(ord('A') + len(sampleBlocks)-1):
            logging.debug("\n\t".join([str(x) for x in output]))
            #logging.debug("SampleBlockStructure " +  " -> ".join([str(x) for x in sampleBlocks]))
            logging.debug('bps used %d' % moveThroughCnt)
            sampleBlocks[0].start = 'source'
            sampleBlocks[-1].end = 'sink'
            return True, sampleBlocks
        return False, None
        
    def annotate(self, refBlocks, sampleBlocks):
        """
        Given the reference and the sample, 
        See if you can explain what's happening in a simpler manner
        Grouped DEL/INS/INV/TLOC events
        """
        #Annotation:
        #Then I should be able to put the full coordinates. A reverse lookup of the blocks and intervals
        nMoves = Counter()
        anno = []
        
        for blockPos, blockId in enumerate(refBlocks):
            if sampleBlocks.count(blockId) == 0:#LOSS
                anno.append('DEL %s' % blockId)
                nMoves[blockPos] -= 1
            
            elif sampleBlocks.count(blockId) == refBlocks.count(blockId):
                totMoves = sum([nMoves[x] for x in nMoves if x < blockPos])   
                #position correction upto this point. if not, i's tloc
                if sampleBlocks.index(blockId) != refBlocks.index(blockId) + totMoves:
                    
                    ref = refBlocks[refBlocks.index(blockId)]
                    sam = sampleBlocks[sampleBlocks.index(blockId)]
                    prev = sampleBlocks[sampleBlocks.index(blockId)-1]
                    
                    pos = refBlocks[refBlocks.index(prev)].end
                    #the s isn't putting in WHERE it's inserted correctly...
                    #then I 
                    anno.append('CNVNeutral %s DEL -> INS %s' \
                                % (sam, pos)) 
                    
                    nMoves[refBlocks.index(blockId)] -= 1
                    nMoves[sampleBlocks.index(blockId)] += 1
            
            elif sampleBlocks.count(blockId) > 1:#GAIN
                lastIndex = 0
                #for every occurance of the refblock
                txt = 'CNVGain %s' % (blockId)
                for j in range(sampleBlocks.count(blockId)):
                    #get it's index (offset by the lastIndex)
                    idx = sampleBlocks[lastIndex:].index(blockId) + lastIndex
                    #if it's in the same position as in refBlocks.. it's no change
                    #this might not actually be true, through. need to check neighbors
                    if idx == refBlocks.index(blockId):
                        txt += ' NoChange %s' % (refBlocks[idx])
                    else:
                        #It's been inserted in this position 
                        txt +=  ' INS %d' % (sampleBlocks[idx-1].end)
                        nMoves[idx] += 1
                    lastIndex = idx + 1
                anno.append(txt)
        
        return anno

def parseArgs(argv):
    parser = argparse.ArgumentParser(prog="Honey.py cpxres", description=USAGE, \
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("tails", metavar="TAILS", type=str, \
                        help="Input hon.tals file")
    parser.add_argument("-o", "--output", type=str, default=None, \
                        help="Output file (<tails>.cpx)")
    parser.add_argument("-c", "--minBlock", type=int, default=500, \
                        help=("To prevent 'tiny' reference bocks, remove "
                              "those with a size less than (%(default)s)"))
    parser.add_argument("-s", "--maxSpan", type=int, default=100000, \
                        help=("Max Span of a breakpoint to be considered"
                              " (%(default)s)"))
    parser.add_argument("-l", "--maxOvl", type=int, default=10, \
                        help=("Max number of overlaps in a cluster"
                              " (%(default)s)"))
    parser.add_argument("-r", "--maxRefBlocks", type=int, default=10, 
                        help=("Max number of reference blocks to consider"
                              " (%(default)s)"))
    parser.add_argument("--debug", action="store_true", \
                        help="Verbose logging")
    args = parser.parse_args(argv)
    setupLogging(args.debug)
    if args.output is None:
        args.output = args.tails + '.cpx'
    return args

def run(argv):
    args = parseArgs(argv)
    cpxres = ComplexResolver(args.tails, args.minBlock, args.maxSpan, args.maxOvl, \
                              args.maxRefBlocks, args.output)
    cpxres.run()

    
if __name__ == '__main__':
    run(sys.argv[1:]) 
