#!/usr/bin/env python

import sys, json
from glob import glob
from optparse import OptionParser

from pbsuite.utils.setupLogging import *
from pbsuite.utils.FileHandlers import *
from pbsuite.jelly.Jelly import JellyProtocol
from pbsuite.jelly.Support import SUPPORTFLAGS

import networkx as nx

USAGE= "Collection.py <protocol.xml>"

def makeFilMetName(a,b):
    """
    
    """
    a = a.replace('/','.')
    b = b.replace('/','.')
    j = filter(lambda x: x != "", [a,b])
    j.sort()
    return "_".join(j)

class FillingMetrics():
    """
    Hey smart guy... put all of the work you do to split pieces of seed names 
    and order seeds and whatever inside of FillingMetics... easier to ship around
    """
    def __init__(self, data, gapName, minReads = 1):
        self.data = data
        self.gapName = gapName
        self.minReads = minReads
        self.__parseData()

    def __parseData(self):
        data = self.data
        gapName = self.gapName
        g = gapName.split('_')
        if len(g) == 1:
            a = g[0]
            b = "-"
        else:
            a,b = g
        self.span = False
        self.singleExtend1 = False
        self.singleExtend2 = False
        self.doubleExtend = False
        self.seed1Name = a
        self.seed1Strand = None
        self.seed1Trim = 0
        
        self.seed2Name = b
        self.seed2Strand = None
        self.seed2Trim = 0
        
        self.fillLength = 0
        self.seed1ExtendSeq = FastqEntry(None, "", "")
        self.seed2ExtendSeq = FastqEntry(None,"", "")
        self.sameStrand = self.seed1Name[-1] != self.seed2Name[-1]
        
        if data.has_key("predictedGapSize"):
            self.predictedGapSize = data["predictedGapSize"]
        else:
            self.predictedGapSize = None
        
        self.failed = data["fillSeq"] is None and data["extendSeq1"] is None and data["extendSeq2"] is None
            
        #Setting support types
        #Setting Strands 
        if SUPPORTFLAGS.span in data["support"][0] and data["fillSeq"] is not None \
           and data["spanCount"] >= self.minReads:
            self.span = True
            self.seed1Strand = '+' if data["spanSeedStrand1"] == '0' else '-'
            self.seed2Strand = '+' if data["spanSeedStrand2"] == '0' else '-'
        elif data["extendF1Count"] >= self.minReads and data["extendF2Count"] >= self.minReads:
            self.doubleExtend = True
            self.seed1Strand = '+' if data["extendF1SeedStrand"] == '0' else '-'
            self.seed2Strand = '+' if data["extendF2SeedStrand"] == '0' else '-'
            if data["extendSeq1"] != None:
                self.seed1ExtendSeq = \
                FastqEntry(b, data["extendSeq1"], "?"*len(data["extendSeq1"]))
            else:
                self.seed1ExtenSeq = FastqEntry(b, "", "")
            if data["extendSeq2"] != None:
                self.seed2ExtendSeq = \
                    FastqEntry(b, data["extendSeq2"], "?"*len(data["extendSeq2"]))
            else:
                self.seed2ExtendSeq = FastqEntry(b, "", "")
        elif data["extendF1Count"] >= self.minReads and data["extendSeq1"] != None:
            self.singleExtend1 = True
            self.seed1ExtendSeq = \
                FastqEntry(b, data["extendSeq1"], "?"*len(data["extendSeq1"]))
            self.seed1Strand = '+' if data["extendF1SeedStrand"] == '0' else '-'
        elif data["extendF2Count"] >= self.minReads and data["extendSeq2"] != None:
            self.singleExtend2 = True
            self.seed2ExtendSeq = \
                FastqEntry(b, data["extendSeq2"], "?"*len(data["extendSeq2"]))
            self.seed2Strand = '+' if data["extendF2SeedStrand"] == '0' else '-'
        else:
            self.failed = 1
        self.seed1Trim = data["seed1Trim"]
        self.seed2Trim = data["seed2Trim"]
        
        self.fillLength  = data["fillBases"]
        self.contribBases = data["contribBases"]
        self.contribSeqs = data["contribSeqs"]
        
        self.spanCount = data["spanCount"]
        self.expandF1Count = data["extendF1Count"]
        self.expandF2Count = data["extendF2Count"]
        
        self.spanSeedScore = abs(int(data["spanSeedScore"]))
        #Is a self circle
        self.isSelfCircle = False
        if self.seed2Name is not None and  self.seed1Name[:10] == self.seed2Name[:10] and not self.isCapturedGap():
            self.isSelfCircle = True
        
        
        
    def isCapturedGap(self):
        """
        checks if assembly folder's name holds captured gap information.
        returns False if it doesn't. 
        returns the gap's name (in a gapinfofile) if it does
        """
        name = self.gapName
        nodes = name.split('_')
        if len(nodes) != 2:
            return False
    
        a, b = nodes
        try:
            aref, aend = a.split('.')
            bref, bend = b.split('.')
        except ValueError:
            return False
        if aref != bref:
            return False
        if not aend.endswith("e3") or not bend.endswith("e5"):
            return False
        c = int(aend[:-2])
        d = int(bend[:-2])
    
        if not c + 1 == d:
            return False
        # return "%s_%d_%d" % (aref, c, d)
        return True
            
    def getSequence(self):
        """
        Returns the full sequence this metric holds
        If it's a span, you'll get a full sequence.
        If it's a gap reduced, You'll get a single sequence
        with a gap ('N') between the pieces 

        Note:
            all gaps < 25bp (overfilled or otherwise) will be inflated
            to 25bp

        returns None if there is any issue with the metrics
        """
        logging.debug("Getting fill sequence for %s" % (self.gapName))
        if self.span:
            logging.debug('fill span')
            return FastqEntry(self.gapName, self.data["fillSeq"].lower(), "?"*len(self.data["fillSeq"]))
        
        if self.predictedGapSize is None:
            #We can't reduce a gap of unknown size
            #one would just get the extend seq
            logging.debug("No predicted gap size")
            return None
        
        gapLen = self.predictedGapSize - self.fillLength
        
        #No improvement at all
        if not (self.span or self.singleExtend1 or self.singleExtend2 or self.doubleExtend):
            logging.debug("Unimproved Gap - %s" % (self.gapName))
            return FastqEntry(self.gapName, ('N'*gapLen), '!'*gapLen)
            
        if gapLen < GAPINFLATE:
            gapLen = GAPINFLATE
        
        #Single end extension from seed 1
        if self.singleExtend1:
            if self.seed1Strand == '-':
                self.seed1ExtendSeq.reverseCompliment()

            logging.debug('single extend seed 1')
            return FastqEntry(self.gapName, 
                        self.seed1ExtendSeq.seq.lower() + \
                        ('N'*gapLen), \
                        self.seed1ExtendSeq.qual + \
                        ('!'*gapLen))

        #Single end extension from seed 2
        if self.singleExtend2:
            if self.seed2Strand == '-':
                self.seed2ExtendSeq.reverseCompliment()

            logging.debug('single extend seed 2')
            return FastqEntry(self.gapName, 
                        self.seed2ExtendSeq.seq.lower() + \
                        ('N'*gapLen), \
                        self.seed2ExtendSeq.qual + \
                        ('!'*gapLen))
 
        #Stick them together!
        #Fill Sequence is on same strand as it should be
        if self.sameStrand and self.seed1Strand == self.seed2Strand:
            logging.debug('same strand success')
            if self.seed1Name.endswith('e3'):
                logging.debug('seed 1 first')
                return FastqEntry(self.gapName,
                            self.seed1ExtendSeq.seq.lower() + \
                            ('N'*gapLen) + \
                            self.seed2ExtendSeq.seq.lower(), \
                            self.seed1ExtendSeq.qual + \
                            ('!'*gapLen) + \
                            self.seed2ExtendSeq.qual)
            elif self.seed2Name.endswith('e3'):
                logging.debug('seed 2 first')
                return FastqEntry(self.gapName,
                            self.seed2ExtendSeq.seq.lower() + \
                            ('N'*gapLen) + \
                            self.seed1ExtendSeq.seq.lower(), \
                            self.seed2ExtendSeq.qual+ \
                            ('!'*gapLen) + \
                            self.seed1ExtendSeq.qual)
            else:
                logging.error(("Huge Problem! This Should Never Happen!  "
                               "sameStrand strandsEqual"))
                exit(10)
        
        #one fill sequence needs to be flipped.
        elif self.sameStrand and self.seed1Strand != self.seed2Strand:
            logging.debug('same strand success - flip one')
            if self.seed1Name.endswith('e5'):
                logging.debug('seed 1 first, flip 2')
                self.seed2ExtendSeq.reverseCompliment()
                return FastqEntry(self.gapName,
                            self.seed1ExtendSeq.seq.lower() + \
                            ('N'*gapLen) + \
                            self.seed2ExtendSeq.seq.lower(), \
                            self.seed1ExtendSeq.qual+ \
                            ('!'*gapLen) + \
                            self.seed2ExtendSeq.qual)
            elif self.seed2Name.endswith('e5'):
                logging.debug('seed 2 first, flip 1')
                self.seed1ExtendSeq.reverseCompliment()
                return FastqEntry(self.gapName,
                            self.seed2ExtendSeq.seq.lower() + \
                            ('N'*gapLen) + \
                            self.seed1ExtendSeq.seq.lower(), \
                            self.seed2ExtendSeq.qual+ \
                            ('!'*gapLen) + \
                            self.seed1ExtendSeq.qual)
            else:
                logging.error(("Huge Problem! This Should Never Happen!  "
                               "sameStrand !strandsEqual"))
                exit(10)
        
        #One Sequence may be None and not sameStrand
        elif not self.sameStrand and (self.seed1Strand is None or self.seed2Strand is None):
            logging.debug("not same strand - one Null")
            if self.seed1Strand is None:
                logging.debug("seed1 None")
                #5' needs to be extended upstream -- 
                if self.seed2Name.endswith('e5'):
                    logging.debug("5' extend upstream")
                    return FastqEntry(self.gapName, \
                                ('N'*gapLen) + \
                                self.seed2ExtendSeq.seq.lower(), \
                                ('!'*gapLen) + \
                                self.seed2ExtendSeq.qual)
                else:
                    #3' needs to be extended downstream
                    logging.debug("3' extend upstream")
                    return FastqEntry(self.gapName, \
                                self.seed2ExtendSeq.seq.lower() + \
                                ('N'*gapLen), \
                                self.seed2ExtendSeq.qual + \
                                ('!'*gapLen)) 
            
            elif self.seed2Strand is None:
                logging.debug("seed2 None")
                if self.seed1Name.endswith('e5'):
                    logging.debug("5' extend upstream")
                    return FastqEntry(self.gapName, \
                                ('N'*gapLen) + \
                                self.seed1ExtendSeq.seq.lower(), \
                                ('!'*gapLen) + \
                                self.seed1ExtendSeq.qual)
                else:
                    logging.debug("3' extend upstream")
                    return FastqEntry(self.gapName, \
                                self.seed1ExtendSeq.seq.lower() + \
                                ('N'*gapLen), \
                                self.seed1ExtendSeq.qual + \
                                ('!'*gapLen)) 
        
        #one fill sequence may need to be filpped
        elif not self.sameStrand:
            logging.debug("not same strand")
            if self.seed1Strand == self.seed2Strand:
                logging.debug('strand equal, flip 2')
                self.seed2ExtendSeq.reverseCompliment()
                if self.seed2Strand == '-':#minus becomes plus
                    self.seed2Strand = '+'
                else:#Plus becomes minus
                    self.seed1Strand = '-'
            if self.seed1Strand == '+':
                logging.debug('strand nequal, 1 +')
                return FastqEntry(self.gapName,
                            self.seed1ExtendSeq.seq.lower() + \
                            ('N'*gapLen) + \
                            self.seed2ExtendSeq.seq.lower(), \
                            self.seed1ExtendSeq.qual+ \
                            ('!'*gapLen) + \
                            self.seed2ExtendSeq.qual)
            elif self.seed2Strand == '+':
                logging.debug('strand nequal, 2 +')
                return FastqEntry(self.gapName,
                            self.seed2ExtendSeq.seq.lower() + \
                            ('N'*gapLen) + \
                            self.seed1ExtendSeq.seq.lower(), \
                            self.seed2ExtendSeq.qual+ \
                            ('!'*gapLen) + \
                            self.seed1ExtendSeq.qual)
            elif self.seed1Strand == '+' and self.seed2Strand == '-':
                logging.error(("Huge Problem! I haven't done this, yet "
                               "Not sameStrand (+/-)"))
                exit(10)
            elif self.seed1Strand == '-' and self.seed2Strand == '+':
                logging.error(("Huge Problem! I haven't done this, yet "
                               "Not sameStrand (-/+)"))
                exit(10)
                return FastqEntry
            else:
                logging.error(("Huge Problem! This Should Never Happen!  "
                               "Not sameStrand"))
                exit(10)
            #I'm just going to take the directStrand sequence 
            # wherever I stich things together later will need to 
            #ensure this is correct
    
    def getTrim(self, contigEnd):
        """
        Get the trim for a specific node
        """
        if contigEnd == self.seed1Name:
            return self.seed1Trim
        if contigEnd == self.seed2Name:
            return self.seed2Trim
        return 0
        
    def getExtendSequence(self, contigEnd):
        """
        Get the sequence that extends the specified contig end. Returns None if 
        this metric doesn't hold extending sequence for the contig end
        """
        if contigEnd == self.seed1Name:
            self.seed1ExtendSeq.seq = self.seed1ExtendSeq.seq.lower()
            return self.seed1ExtendSeq
        if contigEnd == self.seed2Name:
            self.seed2ExtendSeq.seq = self.seed2ExtendSeq.seq.lower()
            return self.seed2ExtendSeq
        return None
        
    def getSeedStrand(self, name):
        if name == self.seed1Name:
            return self.seed1Strand
        elif name == self.seed2Name:
            return self.seed2Strand
        
    def __str__(self):
        return json.dumps(self.data,indent=4)
    
class Collection():

    def __init__(self):
        self.parseOpts()
        #setupLogging(self.debug)
        setupLogging(True)
    
    def parseOpts(self):
        global GAPINFLATE
        parser = OptionParser(USAGE)
        parser.add_option("-m", "--minReads", default=1, type=int, \
                        help=("Minimum number of reads required to fill a "
                              "gap (1))"))
        parser.add_option("-g", "--gapInflate", default=25, type=int, \
                        help=("Minimum size a gap is allowed to be reduced "
                              "or overfilled down to. Smaller gaps will be"
                              " inflated to this size (25)"))
        parser.add_option("--debug", action="store_true", 
                        help="Increases verbosity of logging")
        opts, args = parser.parse_args()
        self.debug = opts.debug
        
        if opts.minReads < 1:
            logging.warning("MinReads set to 1")
            opts.minReads = 1
        
        self.minReads = opts.minReads
        GAPINFLATE = opts.gapInflate
        
        if len(args) != 1:
            parser.error("Error! Incorrect number of arguments")
        
        self.protocol = JellyProtocol(args[0])
    
    
    def loadReference(self):
        """
        Get the reference information 
        """
        myReference = FastaFile(self.protocol.reference)
        self.reference = {}
        if self.protocol.scaffoldQualName is not None:
            myQual = QualFile(self.protocol.scaffoldQualName)
            for key in myReference:
                self.reference[key.split('|')[-1]] = FastqEntry(key.split('|')[-1], myReference[key], myQual[key])
        else:
            for key in myReference:
                self.reference[key.split('|')[-1]] = FastqEntry(key.split('|')[-1], myReference[key], 'l'*len(myReference[key]))
        
        #Graph gml and what not.
        bmlFile = os.path.join(self.protocol.outDir, "assembly", "masterSupport.bml")
        if not os.path.exists(bmlFile):
            logging.error("Consolidated support graph file %s not found!" \
                            % (bmlFile))
            exit(20)
        
        self.gapInfo = GapInfoFile(self.protocol.gapTable)
        self.gapGraph = GapGraph(self.gapInfo)
        self.gapGraph.loadBML(bmlFile)
        self.inputGml = self.gapGraph.graph

    def metricsCollector(self):
        folder = glob(os.path.join(self.protocol.outDir,"assembly","ref*"))
        gapStats = open(os.path.join(self.protocol.outDir, "gap_fill_status.txt"), 'w')
        self.allMetrics = {}
        numGapsAddressed = len(folder)
        noFillingMetrics = 0
        minReadFail = 0
        nfilled = 0
        filled = 0
        overfilled = 0
        reduced = 0
        extended = 0
        trimmed = 0
        
        for f in folder:
            gapName = f.split('/')[-1]
            try:
                fh = open(os.path.join(f,"fillingMetrics.json"),'r')
            except IOError:
                noFillingMetrics += 1
                gapStats.write("%s\tnofillmetrics\n" % gapName)
                continue
            
            try:
                myMetrics = FillingMetrics(json.load(fh), gapName, self.minReads)
                if myMetrics.failed:
                    if myMetrics.failed == 1:
                        minReadFail += 1
                        gapStats.write("%s\tminreadfail\n" % gapName)
                    else:
                        noFillingMetrics += 1
                        gapStats.write("%s\tnofillmetrics\n" % gapName)
                    continue
                else:
                    self.allMetrics[gapName] = myMetrics
            except ValueError:
                logging.error("WARNING! "+f+" didn't produce a valid JSON output in " + \
                            "fillingMetrics.json. Go check this file and if it is " + \
                            "not a plain text JSON file, try re-running the " + \
                            "assembly Process on this folder")
                exit(1)
            fh.close()
            
            if myMetrics.span:
                if myMetrics.data["spanSeedName"] == "tooShortNs":
                    gapStats.write("%s\tn_filled\n" % gapName)
                    nfilled += 1
                else:
                    gapStats.write("%s\tfilled\n" % gapName)
                    filled += 1
            elif myMetrics.predictedGapSize is not None and myMetrics.predictedGapSize < myMetrics.fillLength:
                gapStats.write("%s\toverfilled\n" % gapName)
                overfilled += 1
            elif myMetrics.singleExtend1 or myMetrics.singleExtend2:
                gapStats.write("%s\tsingleextend\n" % gapName)
                extended += 1
            elif myMetrics.doubleExtend:
                gapStats.write("%s\tdoubleextend\n" % gapName)
                reduced += 1        
            if myMetrics.seed1Trim + myMetrics.seed2Trim > 0:
                trimmed += 1
        gapStats.close()
        logging.info("Number of Gaps Addressed %d" % (numGapsAddressed))
        logging.info("No Filling Metrics %d" % (noFillingMetrics))
        logging.info("Min Read Failed %d" % (minReadFail))
        logging.info("Filled %d" % (filled))
        logging.info("NFilled %d" % (nfilled))
        logging.info("Single-End Reduced %d" % (extended))
        logging.info("Double-End Reduced %d" % (reduced))
        logging.info("Overfilled %d" % (overfilled))
        logging.info("Gaps with trimmed edges %d" % (trimmed))
    
    def cleanGraph(self):
        #Need fully connected graphs to get the diameter
        subG = nx.connected_component_subgraphs(self.inputGml)
        #Just in case subG isn't of list type (networkx 1.7)
        try:
            logging.info("PreFilter: %d subGraphs" % (len(subG)))
        except TypeError:
            pass
        
        self.subGraphs = []
        
        for myGraph in subG:
            #Break any edge that didn't make a correct assembly
            # or doesn't span (excluding scaffold/contig evidence)
            for a,b in myGraph.edges():
                name = makeFilMetName(a,b)
                if "Contig" in myGraph.edge[a][b]['evidence']:
                    if name in self.allMetrics.keys():
                        #circles...
                        del(self.allMetrics[name])
                elif "Scaffold" in myGraph.edge[a][b]['evidence']:
                    if name in self.allMetrics.keys():
                        myGraph.node[a]['trim'] = self.allMetrics[name].getTrim(a)
                        myGraph.node[b]['trim'] = self.allMetrics[name].getTrim(b)
                elif name not in self.allMetrics.keys():
                    logging.debug("Breaking %s due to assembly failure" %  name)
                    myGraph.remove_edge(a, b)
                elif not self.allMetrics[name].span:
                    logging.debug("Breaking %s due to non-span" %  name)
                    myGraph.remove_edge(a, b)
                elif self.allMetrics[name].isSelfCircle:
                    logging.debug("Breaking %s due to self-circularity" % name)
                    myGraph.remove_edge(a, b)
    
            #Resolving "forked" nodes
            # Usually caused by some repeat. Right now, it's all about the fill quality
            for node in myGraph.nodes_iter():
                if len(myGraph.edge[node]) > 2:
                    best = None
                    bestScoreSpan= 0
                    bestScoreSeqs = 0
                    bestSpanScore = 0
                    for edge in myGraph.edge[node]:
                        name = makeFilMetName(node,edge)
                        if name in self.allMetrics.keys():
                            data = self.allMetrics[name]
                            seq = data.getSequence()
                            if seq is None:
                                #I fixed this
                                logging.debug("About to Fail on %s (node: %s)" % (name, node))
                            
                            #if len(seq.qual) == 0 or seq is None:
                            if len(seq.qual) == 0:
                                logging.info("NoFilling %s " % name)
                                myScoreSpan = 0
                                myScoreSeqs = 0
                            else:
                                #myScore = data.contribBases / float(data.contribSeqs)
                                myScoreSpan = data.spanCount
                                myScoreSeqs = data.contribBases
                                #sum([ord(y)-33 for y in seq.qual])/float(len(seq.qual))
                                
                            if myScoreSpan > bestScoreSpan:
                                bestScoreSpan = myScoreSpan
                                bestScoreSeqs = myScoreSeqs
                                best = name
                            elif myScoreSeqs == bestScoreSeqs:
                                if myScoreSpan > bestScoreSpan:
                                    bestScoreSpan = myScoreSpan
                                    bestScoreSeqs = myScoreSeqs
                                    best = name
                            
                    logging.debug("Resolved fork to be %s" % best)
                    if best is None:
                        #Again, think I fixed
                        logging.debug("I don't know how this doesn't get set")
                        logging.debug("Node %s" % (node))
                        logging.debug(json.dumps(myGraph.edge[node], indent=4))
                        
                    for edge in list(myGraph.edge[node]):
                        name = makeFilMetName(node, edge)
                        if "Contig" in myGraph.edge[node][edge]['evidence'] \
                            or "Scaffold" in myGraph.edge[node][edge]['evidence']:
                            pass
                        elif name != best:
                            logging.debug("Removed edge %s" % name)
                            myGraph.remove_edge(node,edge)
                
                #add trim everywhere--everything left is either with or without metrics
                for edge in myGraph.edge[node]:
                    if "Contig" in myGraph.edge[node][edge]['evidence']:
                        continue
                    name = makeFilMetName(node,edge)
                    if name not in self.allMetrics.keys():
                        myGraph.node[node]['trim'] = 0
                        myGraph.node[edge]['trim'] = 0
                        continue
                    myGraph.node[node]['trim'] = self.allMetrics[name].getTrim(node)
                    myGraph.node[edge]['trim'] = self.allMetrics[name].getTrim(edge)
                if node in self.allMetrics.keys() and 'trim' not in myGraph.node[node].keys():
                    myGraph.node[node]['trim'] = self.allMetrics[node].getTrim(node)
                    
            #Getting the contig paths
            for i,s in enumerate(nx.connected_component_subgraphs(myGraph)):
                #print "prefilter diameter of testSub piece %d == %d" % (i, nx.diameter(s))
                #I may get an error here if my above cleaning work isn't good enough 
                #for every case
                try:
                    #I don't understand the error I'm getting here...
                    #Doesn't happen anymore
                    ends = nx.periphery(s)#Singletons...
                except AttributeError:
                    logging.debug("Weird error!")
                    logging.debug("Nodes " + json.dumps(myGraph.node))
                    logging.debug("edges " + json.dumps(myGraph.edge))
                    logging.debug("types " + str(type(s)) + " " + str(i) + " " + str(type(i)))
                    logging.debug("Trying again??")
                    ends = nx.periphery(s)
                    
                if len(ends) > 2 and len(ends) == s.number_of_nodes():
                    logging.warning("Circular graph detected. Breaking weakest edge")
                    worstScoreSpan = 0
                    worstScoreSeqs = 0
                    worstEdge = None
                    for a,b in s.edges():
                        if "Contig" in myGraph.edge[a][b]['evidence'] \
                           or "Scaffold" in myGraph.edge[a][b]['evidence']:
                            continue
                        
                        logging.debug( "HERE" )
                        logging.debug(str(a)+" "+str(b))
                        name = makeFilMetName(a,b)
                        data = self.allMetrics[name]
                        myScoreSpan = data.spanCount
                        logging.debug(myScoreSpan)
                        myScoreSeqs = data.contribBases
                        logging.debug(myScoreSeqs)
                            
                        if myScoreSpan > worstScoreSpan:
                            worstScoreSpan = myScoreSpan
                            worstScoreSeqs = myScoreSeqs
                            worstEdge = (a,b)
                        if myScoreSeqs == worstScoreSeqs:
                            if myScoreSpan > worstScoreSpan:
                                worstScoreSpan = myScoreSpan
                                worstScoreSeqs = myScoreSeqs
                                worstEdge = (a,b)
                        
                        
                    logging.info("breaking at %s" % (str(worstEdge)))
                    s.remove_edge(*worstEdge)
                
                #if the above didn't if didn't fix periphery, we'll get
                #a value error and a problem parsing the graph...
                try:
                    a,b = nx.periphery(s)
                except ValueError:
                    logging.error("Graph doesn't have ends. Check it's repeats in collectionErr.gml")
                    logging.error(nx.periphery(s))
                    nx.write_gml(s, "collectionError.gml")
                    exit(1)
                #What I should default to here is just breaking all non captured links

                #nx.shortest_path(s,a,b)
                self.subGraphs.append(s)
    
        logging.info("PostFilter: %d subGraphs" % (len(self.subGraphs)))
    
    def grabContig(self, nodeA, nodeB, graph):
        """
        grabs the contig from the reference that exists
        between nodes A and B
        """
        nodeA, nodeB = orderSeeds(nodeA, nodeB)
        logging.debug("who? %s %s" % (nodeA, nodeB))
        
        try:
            trimA = graph.node[nodeA]['trim']
        except KeyError:
            trimA = 0
        try:
            trimB = graph.node[nodeB]['trim']
        except KeyError:
            trimB = 0

        logging.debug("Grabbing contig between nodes %s & %s - [trim %d, %d]"\
                       % (nodeA, nodeB, trimA, trimB))
        
        scafName = nodeA[:10]
        seq = self.reference[scafName]
        
        #let's get the start
        if nodeA.count('.') == 1:
            #find gap with /0 name
            gid = int(nodeA[nodeA.rindex('.')+1:-2])
            if nodeA.endswith('e3'):
                gapName = "%s_%d_%d" % (scafName, gid, gid+1)
                #trimA = -trimA
            else:
                gapName = "%s_%d_%d" % (scafName, gid-1, gid)
            gap = self.gapInfo[gapName]
            start = gap.end 
        else:#no / means it's got to be the beginning
            start = 0 

        if nodeB.count('.') == 1:
            gid = int(nodeB[nodeB.rindex('.')+1:-2])
            if nodeB.endswith('e3'):
                gapName = "%s_%d_%d" % (scafName, gid, gid+1)
                #trimB = -trimB
            else:
                gapName = "%s_%d_%d" % (scafName, gid-1, gid)
            gap = self.gapInfo[gapName]
            end = gap.start 
        else:# no/ means it's got to be the end
            end = len(seq.seq)
            
        #need a preventer here
        logging.debug("contig %s to %s" % (str(start), str(end)))
        logging.debug("trimming %d and %d" % (trimA, trimB))
        return seq.subSeq(start + trimA, end - trimB)
            
    def outputContigs(self):
        """
        output all the contigs, use the span and stuff get
        """
        fout = open(os.path.join(self.protocol.outDir, "jelly.out.fasta"), 'w')
        qout = open(os.path.join(self.protocol.outDir, "jelly.out.qual" ), 'w')
        liftOverTable = {}#ContigName: [(piece, strand, size), ]
        for part,graph in enumerate(self.subGraphs):
            logging.debug("Contig %d" % part)
            liftTracker = []
            end, start = nx.periphery(graph)
            path = nx.shortest_path(graph, start, end)
            #Change 1 For Filps -- Try moving down normal
            #if not start.endswith("e5"):
                #start, end = start, end
            
            curFasta = []
            curQual = []
            name = makeFilMetName(start, "")
            #First guy's extender
            if name in self.allMetrics.keys():
                logging.debug("First Extender %s" % name)
                data = self.allMetrics[name]
                seq = data.getExtendSequence(name)
                strand = '+'
                if data.seed1Strand == '-':
                    seq.reverseCompliment()
                    strand = '-'
                liftTracker.append((name, strand, len(seq.seq)))
                curFasta.append(seq.seq)
                curQual.append(seq.qual)
            
            #Did we filp the previous sequence
            pFlip = 1  #No = 1
            #Are we putting this together backwards from the start
            FirstFlip = name.endswith('e3')
            
            for i, nodeA in enumerate(path[:-1]):
                nodeB = path[i+1]
                logging.debug("Moving from %s to %s (p=%d)" % (nodeA, nodeB, pFlip))
                name = makeFilMetName(nodeA, nodeB)
                #Existing sequence -- put in A
                if "Contig" in graph.edge[nodeA][nodeB]['evidence']:
                    logging.debug("contig")
                    #need to output the contig seq
                    #Trim Note 1
                    seq = self.grabContig(nodeA, nodeB, graph)
                    strand = '+'
                    if pFlip == -1:
                        strand = '-'
                        seq.reverseCompliment()
                    curFasta.append(seq.seq)
                    curQual.append(seq.qual)
                    liftTracker.append((name, strand, len(seq.seq)))
                #We have to, at the very least, keep a gap in the sequence
                elif "Scaffold" in graph.edge[nodeA][nodeB]['evidence'] and \
                                name not in self.allMetrics.keys():
                    logging.debug("unimproved gap")
                    #keep mat orientation the same
                    a = nodeA[nodeA.rindex('.')+1:-2]
                    b = nodeB[nodeB.rindex('.')+1:-2]
                    j = [int(a),int(b)]; j.sort(); a,b = j
                    gapName = "%s_%d_%d" % (nodeA[:10], a, b)
                    gapSize = self.gapInfo[gapName].length
                    curFasta.append("N"*gapSize)
                    curQual.append("!"*gapSize)
                    liftTracker.append((name, '?', gapSize))
                elif "Scaffold" in graph.edge[nodeA][nodeB]['evidence'] and \
                                name in self.allMetrics.keys():
                    logging.debug("improved gap")
                    data = self.allMetrics[name]
                    seq = data.getSequence()
                    
                    if not data.sameStrand:
                        logging.error(("Gap %s has opposite strand "
                                         "fillers even though they're "
                                         "within scaffold gaps") % name)
                        exit(10)#never happens
                    strand = '+'
                    if (pFlip == -1 and data.getSeedStrand(nodeA) == '+') or \
                       (pFlip == 1 and data.getSeedStrand(nodeA) == '-'):
                        strand = '+'
                        seq.reverseCompliment()
                    liftTracker.append((name, strand, len(seq.seq)))
                    curFasta.append(seq.seq)
                    curQual.append(seq.qual)
                else:
                    #Else we have new sequence joining Scaffolds
                    logging.debug("new sequence")
                    data = self.allMetrics[name]
                    
                    seq = data.getSequence()
                    
                    a = 1 if data.getSeedStrand(nodeA) == '+' else -1
                    b = 1 if data.getSeedStrand(nodeB) == '+' else -1
                    
                    #I've put the first contig in backwards
                    #if FirstFlip is None:
                        #FirstFlip = nodeA.endswith('e3') 
                        #logging.debug("FirstFlip internal %s - %s" % (str(FirstFlip), nodeA))
                    
                    if pFlip == a:
                        m = 1
                    else:
                        m = -1
                    
                    strand = '+'
                    if m == -1:
                        strand = '-'
                        seq.reverseCompliment()
                    liftTracker.append((name, strand, len(seq.seq)))
                    
                    curFasta.append(seq.seq)
                    curQual.append(seq.qual)
                    
                    pFlip = b * m
                    
            name = makeFilMetName(end, "")
            #Final guy's extender
            if name in self.allMetrics.keys():
                data = self.allMetrics[name]
                seq = data.getExtendSequence(name)
                strand = '+'
                #Here -- if we're in the filp and this is on opposite strand
                if pFlip == -1 and data.seed1Strand == '+':
                    strand = '-'
                    seq.reverseCompliment()
                liftTracker.append((name, strand, len(seq.seq)))
                curFasta.append(seq.seq)
                curQual.append(seq.qual)
            
            #We may have been assembling - strand this whole time and we need
            # revcomp it
            if FirstFlip:
                tF = []
                tQ = []
                tL = []
                logging.debug("FirstFlipping %d" % (part))
                """
                for i in curFasta[::-1]:
                    #tF.append(i.translate(revComp)[::-1])
                    tF.append(revComp)
                for i in curQual[::-1]:
                    #tQ.append(i[::-1])
                    tQ.append(i)
                for i in liftTracker[::-1]:
                    name, strand, size = i
                    #if strand == '+':
                        #strand = '-'
                    #elif strand == '-':
                        #strand = '+'
                    tL.append((name, strand, size))
                """
                #""Change 2 -- Should I be iterating this backwards?
                for i in curFasta:
                    tF.append(i.translate(revComp)[::-1])
                for i in curQual:
                    tQ.append(i[::-1])
                for i in liftTracker:
                    name, strand, size = i
                    if strand == '+':
                        strand = '-'
                    elif strand == '-':
                        strand = '+'
                    tL.append((name, strand, size))
                #"
                curFasta = tF
                curQual = tQ
                liftTracker = tL
                
            fout.write(">Contig%d\n%s\n" %  (part, "".join(curFasta)))
            qout.write(">Contig%d\n%s\n" % (part, "".join(curQual)))
            liftOverTable["Contig%d" % (part)] = liftTracker
        
        fout.close()
        qout.close()
        lout = open(os.path.join(self.protocol.outDir, 'liftOverTable.json'),'w')
        json.dump(liftOverTable, lout)
        lout.close()
    
    def run(self):
        logging.info("Grabbing Filling Metrics")
        self.metricsCollector()
        logging.info("Loading Reference")
        self.loadReference()
        logging.info("Removing Poorly Supported Edges")
        self.cleanGraph()
        logging.info("Outputting new contigs")
        self.outputContigs()
        logging.info("Finished!")
        #g = nx.Graph()
        #for i in self.subGraphs:
            #for n in i.nodes():
                #g.add_node(n)
            #for e in i.edges():
                #g.add_edge(*e)
        #nx.write_gml(g,"output.gml")

def orderSeeds(a, b):
    #contig end inside of gap
    if a.count('.') == 1 and b.count('.') == 0:
        if a.endswith('e3'):
            return b, a
        else:
            return a, b
    #contig end inside of gap?
    if b.count('.') == 1 and a.count('.') == 0:
        if b.endswith('e5'):
            return b, a
        else:
            return a, b
    #interscaffold gap
    if a.count('.') == 0 and b.count('.') == 0:
        return a, b
    #captured gap
    if a.endswith('e5'):
        return a, b
    else:
        return b, a
    #doin'g get here
    aint = int(a.split('.')[1][:-2])
    bint = int(b.split('.')[1][:-2])
    if aint < bint:
        return a, b
    return b,a

if __name__ == '__main__':
    c = Collection()
    c.run()
    """
    Load the GapGraph.gml and the reference sequence
    
    for every successful ref_ref, add a link.
    
    resolve redundant links
        
    try to get the maximum diameter and whatever 
    
    traverse from end to end - marking where each piece is/was
        
    out
    """
