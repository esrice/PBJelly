import re, sys, os, bisect, logging, copy
from collections import defaultdict
from StringIO import StringIO
from string import maketrans

from networkx import draw, write_gml, write_graphml, Graph

"""
An object for loading Fasta entries from a file into a
dictionary {SequenceName:sequence}
"""
def wrap(string, width=100):
    return os.linesep.join( \
        [ string[i:i+width] for i in xrange(0,len(string),width) ] )

def qwrap(lst, width = 40):
    buffer = []
    for i in xrange(0,len(lst),width):
        buffer.append( " ".join(map(str,lst[i:i+width])))
    return "\n".join(buffer)

def enum(**enums):
    return type('Enum', (), enums)

class FastaFile(dict):
    
    def __init__(self, fileName):
        super(dict)
        self.fileHandler = open(fileName,'r')
        self.__parse()
        self.fileHandler.close()
    
    def __parse(self):
        for line in self.fileHandler.readlines():
            if line.startswith('>'):
                curName = line.strip()[1:]
                self[curName] = StringIO()
                continue
            self[curName].write(line.strip())
    
        for key in self:
            self[key] = self[key].getvalue()
    
    def toString(self):
        buffer = []
        for key in self:
            buffer.append( ">" + key + "\n" + wrap(self[key]).strip())
        return "\n".join(buffer)

splRE = re.compile("\s+")
class QualFile(dict):
    
    def __init__(self, fileName, convert=True):
        super(dict)
        self.fileHandler = open(fileName,'r')
        self.convert = convert
        self.__parse()
        self.fileHandler.close()
    
    def __parse(self):
        splRE = re.compile("\s+")
        for line in self.fileHandler.readlines():
            if line.startswith('>'):
                curName = line.strip()[1:]
                if self.convert:
                    self[curName] = []
                else:
                    self[curName] = StringIO()
                continue
            
            if self.convert:
                self[curName].extend(map(int, splRE.split(line.strip())))
            else:
                self[curName].write(line.strip()+" ")
        
        if not self.convert:
            for key in self:
                self[key] = self[key].getvalue().strip()
    
    def __parse2(self):
        for line in self.fileHandler.readlines():
            if line.startswith('>'):
                curName = line.strip()[1:]
                continue
            self[curName].write(line.strip() + " ")
        
        for key in self:
            self[key] = map(int, re.split("\s+", self[key].getvalue().strip()))
    
    def toString(self):
        buffer = []
        for key in self:
            buffer.append( ">" + key + "\n" + qwrap(self[key]).strip())
        return "\n".join(buffer)
        
    def valToString(self,key):
        return qwrap(self[key]).strip()

def mergeFastaQual(fastaFile, qualFile):
    """
    opens a fasta and qual file and then merges them to make a fastq dict
    """
    fasta = FastaFile(fastaFile)
    qual = QualFile(qualFile)
    ret = {}
    for key in fasta:
        ret[key] = FastqEntry(key, fasta[key], qual[key])
    return ret
    
class FastqFile(dict):
    
    def __init__(self, fileName):
        self.fileName = fileName
        fh = open(fileName,'r')
        while True:
            name = fh.readline().strip()[1:]
            if name == "": break
            #seq grab
            line = fh.readline().strip()
            seq = StringIO()
            
            while not line.startswith('+'):#Assuming no name...
                seq.write(line)
                line = fh.readline().strip()
            seq = seq.getvalue()
            seqLen = len(seq)
    
            qual = ""
            curLen = 0
    
            while curLen != len(seq):
                line = fh.readline().strip()
                if line == "":
                    sys.stderr.write("Bad Fastq File: Last attempted entry = %s\n" % (name))
                    exit(10)
                curLen += len(line)
                qual += line
            
            self[name] = FastqEntry(name, seq, qual)
    
    def toString(self):
        buffer = []
        for key in self:
            buffer.append(self[key].toString())
        return "".join(buffer)
            
class FastqEntry():
    
    def __init__(self, name, seq, qual):
        self.name = name
        self.seq = str(seq)
        if type(qual) == list:
            self.qual = "".join(map(lambda x: chr(x+33), qual))
        elif type(qual) == str:
            self.qual = qual
        elif type(qual) == unicode:
            self.qual = str(qual)
        else:
            raise TypeError("FastqEntry.qual needs to be list, str or unicode. Not %s" % (type(qual)))
        
    def reverseCompliment(self):
        """
        reverse compliments this sequence
        """
        self.seq = self.seq.translate(revComp)[::-1]
        self.qual = self.qual[::-1]
       
    def lowerCase(self):
        """
        run tolower on the sequence
        """
        self.seq = self.seq.tolower()

    def upperCase(self):
        """
        """
        self.seq = self.seq.toupper()

    def getSeq(self,name, start=0, end = None):
        """
        Helper Method to rename the seq
        """
        return "@%s\n%s\n+\n%s\n" % (name, self.seq[start:end], self.qual[start:end])
        
    def translateQual(self):
        """
        Returns the expanded, non encoded qual as a list
        """
        ret = []
        for i in self.qual:
            ret.append(ord(i)-33)
        return ret
        
    def subSeq(self, start, end):
        """
        Return a new subsequence of this fastq entry
        """
        seq = self.seq[start:end]
        qual = self.qual[start:end]
        return FastqEntry(self.name, seq, qual)
    
    def toString(self, start=0, end = None):
        """
        Trim before returing entry as a string
        """
        if start != 0 or end is not None:
            name = self.name +  ("##%d#%d##" % (start, end))
            return "@%s\n%s\n+\n%s\n" % (name, self.seq[start:end], self.qual[start:end])
        #quicker?
        return "@%s\n%s\n+\n%s\n" % (self.name, self.seq, self.qual)
        
    def __str__(self):
        return "@%s\n%s\n+\n%s\n" % (self.name, self.seq, self.qual)

"""
An object for handling gapInfo files, which are in bed format:

targetName \t start \t end \t gapName
"""

class GapInfoFile(dict):
    
    def __init__(self, fileName):
        super(dict)
        self.fileHandler = open(fileName,'r')
        #inorder lift of gaps
        for line in self.fileHandler.readlines():
            curGap = Gap(*line.strip().split('\t'))
            self[curGap.name] = curGap
        
        self.fileHandler.close()
        
    def getSortedGaps(self):
        """
        Returns dictionary of gaps partitioned by scaffold and sorted by
        location
        """
        ret = defaultdict(list)
        for key in self:
            if self[key].start == 'na':
                continue
            #print str(self[key])
            #print ret[self[key].scaffold]
            bisect.insort(ret[self[key].scaffoldId], self[key])
            #print ret[self[key].scaffold]
            #raw_input()
            
        return dict(ret)
        
    def getScaffFromIndex(self, key):
        """
        Looks through file for scaffold name
        based on scaffold index
        """
        for x in self.keys():
            if x.startswith(key):
                return self[x].scaffold
        
        raise KeyError(key)

class Gap():
    #ends flags
    BEGIN = 1
    END = 2
    def __init__(self, scaff, start, end, name, endsFlag=0):
        #scaff, start, end, name = line.strip().split(delim)
        ref, lcontig, rcontig = name.split('_')
        self.name = name
        self.scaffold = scaff
        self.scaffoldId = ref
        self.leftContig = self.scaffoldId+"."+lcontig
        self.rightContig = self.scaffoldId+"."+rcontig
        
        if start == 'na':
            self.start = 'na'
            self.end = 'na'
            self.length = 0
        else:
            self.start = int(start)
            self.end = int(end)
            self.length = self.end - self.start
        self.endsFlag = int(endsFlag)
       
    def __str__(self):
        return "\t".join([self.scaffold, str(self.start), str(self.end), self.name, str(self.endsFlag)])

    def __lt__(self, other):
        if type(other) == type(self):
            return self.start - other.start < 0
        elif type(other) == int:
            return self.start - other < 0
        else:
            raise AttributeError
    
    def __gt__(self, other):
        if type(other) == type(self):
            return self.start - other.start > 0
        elif type(other) == int:
            return self.start - other > 0
        else:
            raise AttributeError
    
class GapGraph():
    """
    Holds Nodes and Edges of Contig Ends and their Connections.
    """
    def __init__(self, gapInfo=None):
        self.gapInfo = gapInfo
        self.graph = Graph()
        if gapInfo is not None:
            self.__createGraph()
        
    def __createGraph(self):
        """
        Create a graph with 
            Nodes = Contig Ends (5 and 3) and an arbitrary gap name
            Edges = inScaffolding A super heavy PointTo 
                    And reads that map across
            Taking these notes to testing/alignTestNotes.txt

        Create a graph with nodes == contigNames, Edges == 5' or 3'
        """
        for key in self.gapInfo:
            gap = self.gapInfo[key]
            
            if gap.endsFlag == (Gap.BEGIN + Gap.END):
                if gap.start == 'na': # singleton
                    prevNode = gap.scaffoldId + "e5"
                    nextNode = gap.scaffoldId + "e3"
                    self.graph.add_node(prevNode, extenders=[])
                    self.graph.add_node(nextNode, extenders=[])
                    self.graph.add_edge(prevNode, nextNode,  evidence=["Contig"])
                else:#one gap weirdo
                    startNode = gap.scaffoldId + "e5"
                    prevNode = gap.leftContig + "e3"
                    nextNode = gap.rightContig + "e5"
                    endNode = gap.scaffoldId + "e3"
                    
                    self.graph.add_node(startNode, extenders=[])
                    self.graph.add_node(prevNode,  extenders=[])
                    self.graph.add_node(nextNode,  extenders=[])
                    self.graph.add_node(endNode,   extenders=[])

                    self.graph.add_edge(startNode, prevNode, evidence=["Contig"])
                    self.graph.add_edge(prevNode,  nextNode, evidence=["Scaffold"])
                    self.graph.add_edge(nextNode,  endNode,  evidence=["Contig"])
                    
                continue
                
            prevNode = gap.leftContig + "e3"
            if gap.endsFlag & Gap.BEGIN:#is first gap - get that first contig   
                startNode = gap.scaffoldId + "e5"
                self.graph.add_node(startNode, extenders=[])
                self.graph.add_node(prevNode,  extenders=[])
                self.graph.add_edge(startNode, prevNode, evidence=["Contig"])
                
            nextNode = gap.rightContig + "e5"
            if gap.endsFlag & Gap.END:#is last gap
                endNode = gap.scaffoldId + "e3"
            else:
                endNode = gap.rightContig + "e3"
            
            self.graph.add_node(nextNode, extenders=[])
            self.graph.add_node(endNode,  extenders=[])
            
            self.graph.add_edge(prevNode, nextNode, evidence=["Scaffold"])
            self.graph.add_edge(nextNode, endNode,  evidence=["Contig"])
            
    def saveBML(self, fileName="briefGraph.bml"):
        """
        Only saves the evidence and extenders for the graph.
        This is for programatically building a graph quickly.
        evidence  start   end name,name,name
        extend  start   name,name,name
        """
        fout = open(fileName,'w')
        newGraph = copy.deepcopy(self.graph)
        for node in newGraph.nodes_iter():
            if len(newGraph.node[node]['extenders']) > 0:
                if newGraph.node[node].has_key('extenders'):
                    fout.write("extend\t%s\t%s\n" % (node, \
                        "::".join(newGraph.node[node]['extenders'])))
        for edgeA, edgeB in newGraph.edges_iter():
            ev = filter(lambda x: x != "Scaffold" and x != "Contig", \
                    newGraph[edgeA][edgeB]['evidence'])
            if len(ev) > 0:
                fout.write("evidence\t%s\t%s\t%s\n" % (edgeA, edgeB, "::".join(ev)))
        
        fout.close()
        
    def loadBML(self, fileName):
        """
        load the evidence and extenders from the BML file
        """
        fh = open(fileName,'r')
        for line in fh.readlines():
            data = line.strip().split('\t')
            if data[0] == 'extend':
                self.add_extend(data[1].replace('/','.'), data[2].split('::'))
            if data[0] == 'evidence':
                self.add_evidence(data[1].replace('/','.'), data[2].replace('/','.'), data[3].split('::'))

        fh.close()
        
    def saveGraph(self, fileName="graph.gml", spanOnly=False): 
        """
        Save a gml (default fileName is graph.gml). 
        if spanOnly, we won't write extenders evidence
        """
        newGraph = copy.deepcopy(self.graph)
        for node in newGraph.nodes_iter():
            if spanOnly:
                newGraph.node[node]['extenders'] = ""
            else:
                newGraph.node[node]['extenders'] = ":".join(newGraph.node[node]['extenders'])
                
        for edge in newGraph.edges_iter():
            newGraph[edge[0]][edge[1]]['evidence'] = ":".join(newGraph.get_edge_data(*edge)['evidence'])
        
        write_gml(newGraph, fileName)
      
    def add_extend(self, nodeName, extName):
        """
        Takes a string or a list
        """
        logging.debug("%s extends %s" % (extName, nodeName))
        if not isinstance(extName, list):
            extName = [str(extName)]
        
        logging.debug(nodeName)
        if nodeName not in self.graph.node:
            self.graph.add_node(nodeName, extenders=extName)
        else:
            self.graph.node[nodeName]['extenders'].extend(extName)
         
    def add_evidence(self, source, target, readName):
        """
        Add linkage support to the graph
        for span, we add left and then right
        same for contains
        readName can be a string or a list
        """
        if not isinstance(readName, list):
            readName = [str(readName)]
        
        if source == target:
            logging.warning(("Read %s self-extends Node %s " \
                             "Possible Evidence of Tandem Repeat on Singleton")\
                             % (readName, source))
        
            self.add_extend(source, readName)
            return
        
        logging.debug("%s gives evidence %s -> %s" % (readName, source, target))
        
        if source not in self.graph.node:
            self.graph.add_node(source, extenders=[])
        if target not in self.graph.node:
            self.graph.add_node(target, extenders=[])
        
        try:
            l = [source, target]; l.sort(); source, target = l
            self.graph.edge[source][target]['evidence'].extend(readName)
        except KeyError:#New edge
            self.graph.add_edge(source, target, evidence=readName)
     

"""
An object for handling m4 format mapped reads information efficently.

m4 format:

qname tname score pctsimilarity qstrand qstart qend qseqlength 
tstrand tstart tend tseqlength [mapqv ncells npaths]
"""

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
        self.mapqv          = int(data[12])
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
                                  self.tseqlength, self.mapqv]))
    
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
"""
Parses M5 alignment object

"""

revComp = maketrans("ATCGNatcgn","TAGCNtagcn")

class M5File(list):
    
    def __init__(self, file):
        super(list)
        if type(file) != type(StringIO()):
            file = open(file,'r')
        self.fileHandler = file
        self.__parse()
        self.fileHandler.close()
    
    def __parse(self):
        """
        Returns a list of all the M5 lines
        """
        for line in self.fileHandler.readlines():
            try:
                self.append(M5Line(line.strip()))
            except TypeError, err:
                sys.stderr.write("BadM5Line! \n%s\n%s\n" % (line, str(err)) )
                sys.exit(1)
            
class M5Line():
    
    def __init__(self, line):
        data = re.split("\s+",line)

        self.qname      = data[0]
        self.qseqlength = int(data[1])
        self.qstart     = int(data[2])
        self.qend       = int(data[3])
        self.qstrand    = data[4]
        self.tname      = data[5]
        self.tseqlength = int(data[6])
        self.tstart     = int(data[7])
        self.tend       = int(data[8])
        self.tstrand    = data[9]
        self.score      = int(data[10])
        self.nMatch     = int(data[11])
        self.nMismatch  = int(data[12])
        self.nInsert    = int(data[13])
        self.nDelete    = int(data[14])
        self.mapqv      = int(data[15])
        self.querySeq   = data[16]
        self.compSeq    = data[17]
        self.targetSeq  = data[18]
        
        self.queryPctAligned = (self.qend - self.qstart)/float(self.qseqlength)
        self.pctsimilarity = self.nMatch / float(self.qend - self.qstart)
        #"""newBlasr
        if self.tstrand == '-':
            self.negStrand = True
            #translating to + strand.
            self.targetSeq = self.targetSeq.translate(revComp)[::-1]
            self.querySeq = self.querySeq.translate(revComp)[::-1]
            self.compSeq = self.compSeq[::-1]
            self.tstrand = '1'
        else:
            self.negStrand = False
            self.tstrand = '0'
        
        #M5 is now always in + strand orientation
        self.qstrand = '0' if self.qstrand == '+' else '1'        
        self.flag = 0
        self.trim = False
    
    def toBed(self):
        """
        Returns the M5Line to a bed file
        """
        if self.negStrand:
            strand = "-"
            chromStart = str(self.tstart - (self.qseqlength - self.qend))
            chromEnd = str(self.tend + self.qstart)
        else:
            strand = "+"
            chromStart = str(self.tstart - self.qstart)
            chromEnd = str(self.tend + (self.qseqlength - self.qend))

        chrom = self.tname
        name = self.qname
        score = str(self.score)
        thickStart = str(self.tstart)
        thickEnd = str(self.tend)
        itemRgb = str(int(self.pctsimilarity))
        return "\t".join([chrom, chromStart, chromEnd, name, score, strand, \
                          thickStart, thickEnd, itemRgb])

    
    def __str__(self):
        """
        Undo changes 
        """
        #"""newBlasr
        if self.negStrand:
            targetSeq = self.targetSeq.translate(revComp)[::-1]
            querySeq = self.querySeq.translate(revComp)[::-1]
            compSeq = self.compSeq[::-1]
            tstrand = '-'
        else:
            targetSeq = self.targetSeq
            querySeq = self.querySeq
            compSeq = self.compSeq
            tstrand = '+'
        
        qstrand = '+' if self.qstrand == '0' else '-'        
        return " ".join(map(str, [self.qname, self.qseqlength, self.qstart, \
                                  self.qend, qstrand, self.tname, \
                                  self.tseqlength, self.tstart, self.tend, \
                                  tstrand, self.score, self.nMatch, \
                                  self.nMismatch, self.nInsert, self.nDelete, \
                                  self.mapqv, \
                                  querySeq, compSeq, targetSeq]))

class LiftOverTable():
    """
    TODO: 
        Make the entry take care of updating it's own coordinates
        instead of relying on the calling code to calculate shifts.
    """
    def __init__(self, fn=None):
        #Quick Look on a per entry basis
        self.hash = {}
        #For outputting all scaffolds later
        #And keeping contigs contained within 
        #scaffolds for updates
        self.scaffoldRoots = {}
        #For adding
        self.curRoot = None
        if fn is not None:
            self.__parse(fn)

    def __parse(self, fn):
        fh = open(fn,'r')
        head =  fh.readline()
        for line in fh.readlines():
            entry = LiftOverEntry(*line.strip().split('\t'))
            self.addEntry(entry)
        fh.close()

    def addEntry(self, entry):
        if not self.scaffoldRoots.has_key(entry.scaffold) or self.scaffoldRoots[entry.scaffold] is None:
            self.scaffoldRoots[entry.scaffold] = entry
            self.curRoot = entry
        else:
            self.curRoot.next = entry
            entry.prev = self.curRoot
            self.curRoot = entry
        
        key = entry.scaffold+str(entry.oStart)
        self.hash[key] = entry
            
    def getEntry(self, scaffold, oStart):
        return self.hash[scaffold+str(oStart)]
    
    def removeEntry(self, entry):
        if not self.scaffoldRoots.has_key(entry.scaffold):
            raise KeyError("Scaffold %s Not Found" % entry.scaffold)
        
        if entry.prev is not None:
            entry.prev.next = entry.next
        
        if entry.next is not None:
            entry.next.prev = entry.prev
        #I need a better hashing function
        #key = entry.scaffold+str(entry.oStart)
        #del(self.hash[key])
        
        #if this is the only guy in the entire scaffolding
        #Get rid of him
        if self.scaffoldRoots[entry.scaffold].next is None:
            del(self.scaffoldRoots[entry.scaffold])
    
    def updateScaffold(self, entry, shift):
        """
        Takes the key(entry) and changes all
        downstream coordinates as applicable
        """
        #entry.nStart += startShift
        #The previous end changes. 
        #With the start's shift
        #entry.prev.nEnd += startShift
        
        while entry.next is not None:
            entry = entry.next
            if type(entry.nStart) == str:
                continue
            entry.nStart += shift
            entry.nEnd += shift

    def insertEntry(self, existingEntry, newEntry, after=True):
        """
        Inserts a new LiftOverEntry into the Table
        #Note! Changes made to the scaffold before entry is inserted
        are not applied to this entry - so get everything inserted
        before you do this
        """
        key = newEntry.scaffold+str(newEntry.oStart)
        self.hash[key] = newEntry
        
        if after == True:
            newEntry.prev = existingEntry 
            newEntry.next = existingEntry.next 
            if newEntry.next is not None:
                newEntry.next.prev = newEntry
            existingEntry.next = newEntry
        else:
            newEntry.next = existingEntry
            newEntry.prev = existingEntry.prev
            if newEntry.prev is not None:
                newEntry.prev.next = newEntry
            existingEntry.prev = newEntry
            
        
    def __str__(self):
        """
        returns the table as a string
        """
        ret = ""
        for entry in self:
            ret += str(entry) + "\n"
        
        return ret

    def __iter__(self):
        for key in self.scaffoldRoots.keys():
            root = self.scaffoldRoots[key]
            while root.next is not None:
                yield root
                root = root.next
            yield root
        
class LiftOverEntry():
    def __init__(self, scaffold, oStart, oEnd, nStart, nEnd, gType, \
                 prev = None, next = None):
        self.scaffold = scaffold
        
        if oStart == 'na':
            self.oStart = oStart
            self.oEnd = oEnd
        else:
            self.oStart = int(oStart)
            self.oEnd = int(oEnd)

        if nStart == 'na':
            self.nStart = nStart
            self.nEnd = nEnd
        else:
            self.nStart = int(nStart)
            self.nEnd = int(nEnd)
        
        self.next = next
        self.gType = gType
        self.prev = prev
        
    def getNext(self, gType):
        """
        return the next feature that is gType
        """
        if self.next is None:
            return None
        if self.next.gType == gType:
            return self.next
        else:
            return self.next.getNext(gType)

    def getPrev(self, gType):
        """
        return the next feature that is gType
        """
        if self.prev is None:
            return None
        if self.prev.gType == gType:
            return self.prev
        else:
            return self.prev.getPrev(gType)
            
    def __str__(self):
        return "\t".join([self.scaffold, str(self.oStart), str(self.oEnd), \
                          str(self.nStart), str(self.nEnd), self.gType])



