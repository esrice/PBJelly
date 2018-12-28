#!/usr/bin/env python

import sys, glob, os, random
from optparse import OptionParser
from collections import defaultdict, namedtuple


from pbsuite.utils.FileHandlers import *
from pbsuite.jelly.Jelly import JellyProtocol

import networkx

USAGE = """
Consolidates all the reads that support a gap into a folder for assembly.

Extraction.py <Protocol.xml>"""

# Max number of gaps to hold in memory until all output is flushed
# Tweak this based on the amount of memory you have available.
MAXGAPHOLD = 20000;
# Benchmarking needs to be done to guide this number

# Amount of flank to put in from the reference for "Guided" Assembly
FLANKAMT = 2500;
MAXFILES = 10000

TrimInfo = namedtuple("ReadwithTrim", "name start end")

class Extraction():
    
    def __init__(self, args):
        #ArgParse
        self.__parseArgs__(args)
        self.__initLog()
        protocol = JellyProtocol(sys.argv[1])
        self.gapInfo = GapInfoFile(protocol.gapTable)
        self.gapGraph = GapGraph(self.gapInfo)
    
    def __parseArgs__(self, argv):
        parser = OptionParser(USAGE)
        parser.add_option("--debug",action="store_true",default=False)
        self.options, args = parser.parse_args(argv)
            
        if len(args) < 2:
            parser.error("Expected one argument, the Protocol.xml")
        self.protocol = JellyProtocol(args[1])
    
    def __initLog(self):
        """Logging"""
        logLevel = logging.DEBUG if self.options.debug else logging.INFO
        logFormat = "%(asctime)s [%(levelname)s] %(message)s"
        logging.basicConfig( stream=sys.stderr, level=logLevel, format=logFormat )
        logging.info("Running %s" % " ".join(sys.argv) )
    
    def __cleanReadName__(self, readName):
        """
        Build a TrimInfo namedtuple for this read
        """
        getSub = re.compile("(.*)##(\d+)#(\d+)##$")
        #Check for trims
        if readName.endswith('##'):#just endswith? not safe..
            logging.debug("%s -- Missed Adapter" % readName)
            name, start, end = getSub.match(readName).groups()
            start = int(start)
            end = int(end)
        else:
            name = readName
            start = 0
            end = None
        
        return TrimInfo(name, start, end)
    
    def __loadGMLFiles__(self):
        """
        Iterates through all of the files inside of support directory
        merges these graphs with the whole gapGraph
        creates a dictionary of {"readName":[gapName or nodeName,...] }
        """
        self.readSupport = defaultdict(list)
        for i in glob.glob(os.path.join(self.protocol.outDir,"support","*.gml")):
            truncNxVer = float('.'.join(networkx.__version__.split('.')[:2]))
            if truncNxVer >= 1.7:
                try:
                    inputGml = networkx.read_gml(i, relabel=True)
                except ValueError:
                    logging.warning("GML file %s is empty" % i)
                    continue
            elif truncNxVer == 1.1:
                inputGml = networkx.read_gml(i)
            else:
                logging.warning("It is unknown if networkx version %s will work." % networkx.__version__)
                logging.warning("If you get an error here, please report it!!!")
                inputGml = networkx.read_gml(i)
                           
            for node in inputGml.nodes_iter():
                for readName in inputGml.node[node]['extenders'].split(':'):
                    if readName == '':
                        continue
                    trimInfo = self.__cleanReadName__(readName)
                    self.gapGraph.add_extend(node, trimInfo.name)
                    self.readSupport[trimInfo.name].append((node, trimInfo))
                    
            for source, target, evidence in inputGml.edges_iter(data=True):
                for readName in evidence['evidence'].split(':'):
                    if readName == '':
                        continue
                    trimInfo = self.__cleanReadName__(readName)
                    self.gapGraph.add_evidence(source, target, trimInfo.name)
                    self.readSupport[trimInfo.name].append((source, target, trimInfo))
        self.readSupport = dict(self.readSupport)
        
    def __loadReference__(self):
        """
        Get the reference information into memory
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
        
    def __extractReads__(self):
        """
        readsForExtraction = set(readSupport.keys())
        for fastaFile to process:
            readsToCollect = readsForExtraction & set(fastaFile.keys())
            for usedRead in readsToCollect:
                #check the gapFile holder from previous extraction
                write to the file handle of the gaps in each used read
                It might make it progressively more efficient to create set(readSupport.keys())^usedReads
                readsForExtraction = 
            readsForExtraction = readsForExtraction^usedReads
        done.
        """
        self.gapOutputs = {}
        
        self.supportFolder = os.path.join(self.protocol.outDir,"assembly")
        try:
            os.mkdir(self.supportFolder)
        except OSError:
            pass
        
        outputQueue = defaultdict(list)# gapName: [readDataAsString]
        numReads = 0
        for inputFile in self.protocol.inputs:
            logging.info("Parsing %s" % inputFile)
            
            #Check if it is fasta/qual or fastq
            if inputFile.lower().endswith(".fastq"):
                inputReads = FastqFile(inputFile)
                #spaces in the read names...
                for i in inputReads.keys(): #protecting for spaces
                    inputReads[i.split(' ')[0]] = inputReads[i]
                    inputReads[i.split(' ')[0]].name = i.split(' ')[0]
                    #This messes up the names a tiny bit permanently since we're pointing
                    #at the object, but I'm okay with that because we don't need spaces anyway
            elif inputFile.lower().endswith(".fasta"):
                inputReads = self.selectiveMergeFastaQual(inputFile, \
                                inputFile[:inputFile.rindex('.')]+".qual")
            else:
                logging.error("Input Reads file %s doesn't end with .fasta or .fastq!")
                exit(0)
             
            logging.info("Loaded %d Reads" % (len(inputReads.keys())))
            parsed = 0
            for usedRead in inputReads.keys():
                try:
                    gaps = self.readSupport[usedRead]
                except KeyError:
                    continue
                parsed += 1
                for gap in gaps:
                    #We have an extender, add the read into all edges files
                    if len(gap) == 2:
                        contigEnd, trimInfo = gap
                        start, end = trimInfo.start, trimInfo.end
                        # -- If there isn't a single edge, then we need to populate an extender
                        numNodes = 0
                        for source, target in self.gapGraph.graph.edges(contigEnd):
                            if "Contig" not in self.gapGraph.graph.edge[source][target]['evidence'][0]:
                                numNodes += 1
                                #Ensure we don't make redundancies
                                l = [source, target]; l.sort(); source, target = l;
                                gapName = "%s_%s" % (source, target)
                                logging.debug("Writing %s to %s" % (usedRead, gapName))
                                numReads += 1
                                outputQueue[gapName].append(inputReads[usedRead].toString(start, end))
                        if numNodes == 0:
                            logging.debug("%s only has extending evidence" % (contigEnd))
                            logging.debug(self.gapGraph.graph.edges(contigEnd))
                            #self.__gapOutputs__(contigEnd, inputReads[usedRead].toString(start, end))
                            numReads += 1
                            outputQueue[contigEnd].append(inputReads[usedRead].toString(start, end))
                            
                    #we have an edge, so just write to the gap
                    elif len(gap) == 3:
                        source, target, trimInfo = gap
                        l = [source, target]; l.sort(); source, target = l;
                        gapName = "%s_%s" % (source, target)
                        start, end = trimInfo.start, trimInfo.end
                        #self.__gapOutputs__(gapName, inputReads[usedRead].toString(start, end))
                        numReads += 1
                        outputQueue[gapName].append(inputReads[usedRead].toString(start, end))
                
                if len(outputQueue.keys()) >= MAXGAPHOLD:
                    logging.info("Flushing Output Queue of %d gaps %d reads" % \
                                  (len(outputQueue.keys()), numReads))
                    self.flushQueue(outputQueue)
                    logging.info("Finshed Flush")
                    numReads = 0
                    del(outputQueue); 
                    outputQueue = defaultdict(list)
            
            logging.info("Parsed %d Reads" % (parsed))
        
        self.flushQueue(outputQueue)
        
    def flushQueue(self, outputQueue):
        """
        flushes the entire queue - dumps a file at a time
        closes up the file when it's done
        """
        for gapName in outputQueue:
            try:
                outFile = self.gapOutputs[gapName]
                #exists, it's a string
                outFile = open(outFile,'a')
                self.gapOutputs[gapName] = outFile
            except KeyError:
                #returns a file handler
                outFile = self.openGapOut(gapName)
            
            #No newline necessary
            outFile.write("".join(outputQueue[gapName]))
            
            outFile.close()
            self.gapOutputs[gapName] = outFile.name
        
    def openGapOut(self, gapName):
        """
        prepares the gap's output
        """
        #Create output file
        basedir = os.path.join(self.supportFolder, gapName.replace('/','.'))
        try:
            os.mkdir(basedir)
        except OSError as e:
            if e.errno == 17:
                logging.warning("%s already exists... overwriting results" % basedir)
            elif e.errno == 13:
                logging.critical("%s cannot be created... insufficent file permissions" % basedir)
                exit(13)
        
        fh = open(os.path.join(basedir,"input.fastq"),'w')
        self.gapOutputs[gapName] = fh   
        
        #Need to extract Seed Contigs
        #i = gapName.split('_')
        for flankName in gapName.split('_'):
            logging.debug("flankName %s" % flankName)
            sequence = self.reference[flankName[:10]]
            sequenceLength = len(sequence.seq)
            
            #contig end
            if flankName.count('.') == 0:
                if flankName[-1] == '5':
                    seq = sequence.getSeq(flankName, 0, FLANKAMT)
                elif flankName[-1] == '3':
                    s = max(sequenceLength-FLANKAMT, 0)
                    seq = sequence.getSeq(flankName, s, sequenceLength)
            #captured gap
            elif flankName.count('.') == 1:
                ref, part = flankName.split('.')
                part, end = part.split('e')
                if end == '5':
                    #-- This is what will messup during scaff incorp to cap gap -- ref\d{7}_-1_0
                    gap = "%s_%s_%s" % (ref, int(part)-1, part)
                    start = self.gapInfo[gap].end
                    end = min(start+FLANKAMT, sequenceLength)
                    seq = sequence.getSeq(flankName, start, end)
                    
                if end == '3':
                    gap = "%s_%s_%s" % (ref, part, int(part)+1)
                    end = self.gapInfo[gap].start
                    start = max(end - FLANKAMT, 0)#Prevent - start index
                    seq = sequence.getSeq(flankName, start, end)
            seq = self.cleanSeq(seq)
            fh.write( seq )
        
        return fh
          
    def cleanSeq(self, seq):
        seq = seq.replace('M','C')
        seq = seq.replace('R','A')
        seq = seq.replace('W','T')
        seq = seq.replace('S','G')
        seq = seq.replace('Y','C')
        seq = seq.replace('K','T')
        seq = seq.replace('V','G')
        seq = seq.replace('H','A')
        return seq
        
    def run(self):
        #Merge GML files
        logging.info("Opening GML Files")
        self.__loadGMLFiles__()
        #Figure out Support levels of evidence that we'll allow
        logging.info("Loading Reference Sequence")
        self.__loadReference__()
        #Figure out an efficient way of opening inputJobs and extracting necessary reads
        logging.info("Extracting Reads")
        self.__extractReads__()
        self.gapGraph.saveBML(os.path.join(self.supportFolder,"masterSupport.bml"))
        logging.info("Finished")

    def selectiveMergeFastaQual(self, fastaFn, qualFn):
        """
        Merges a fasta/qual but is selective about what gets the full processing
        selection is a dictionary of { readNames : True } for lookup
        I'm trying to speed up Extraction as much as possible
        """
        splRE = re.compile("\s+")
        fasta = FastaFile(fastaFn)
        qual = QualFile(qualFn, convert=False)
        logging.info("Selective loading from %d reads" % (len(fasta.values())))
        ret = {}
        for key in fasta:
            try:
                exists = self.readSupport[key]
                q = "".join([chr(int(x)+33) for x in splRE.split(qual[key])])
                ret[key] = FastqEntry(key, fasta[key], q)
            except KeyError:
                #doesn't exist
                continue
        return ret
    
if __name__ == '__main__':
    me = Extraction(sys.argv)
    me.run()
