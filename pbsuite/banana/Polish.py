#!/usr/bin/env python
import sys, re, random, argparse, textwrap, logging
import operator
from collections import defaultdict, namedtuple

from pbsuite.utils.FileHandlers import FastaFile, M5File, FastqFile
from pbsuite.utils.CommandRunner import exe
from pbsuite.utils.setupLogging import setupLogging

USAGE = """\
Takes reads.fastq and ref.fasta 
maps with blasr
creates consensus
"""

def blasr(query, target, nproc=1, bestn=1, outName="map.m5"):
    """
    runs blasr
    """
    r,o,e = exe("blasr %s %s -bestn %d -affineAlign -m 5 -nproc %d -out %s" \
                % (query, target, bestn, nproc, outName))
    

def realign(alignment):
    """
    realigns target, query so that every alignment should 
    have the comparable characteristics in the same sequence context
    regardless of differences due to the error rate and alignment 
    variation
    realignment happens inplace
    """
    def inner(align):
        target = list(align.targetSeq)
        query = list(align.querySeq)
        compSeq = list(align.compSeq)
        for pos in range(len(align.targetSeq)):
            if align.targetSeq[pos] == '-':
                i = re.match('-+[ATCGatcg]?', align.targetSeq[pos:])
                if align.targetSeq[pos+i.end()-1] == align.querySeq[pos]:
                    target[pos] = align.targetSeq[pos+i.end()-1]
                    compSeq[pos] = '|'
                    target[pos+i.end()-1] = '-'
                    compSeq[pos+i.end()-1] = '*'
                    align.targetSeq = "".join(target)
                    align.compSeq = "".join(compSeq)
    inner(alignment)
    alignment.targetSeq, alignment.querySeq = alignment.querySeq, alignment.targetSeq
    inner(alignment)
    alignment.targetSeq, alignment.querySeq = alignment.querySeq, alignment.targetSeq

def offSetSeqs(offset, align):
    q = [" "]*offset; q.extend(list(align.querySeq))
    c = [" "]*offset; c.extend(list(align.compSeq))
    t = [" "]*offset; t.extend(list(align.targetSeq))
    return q,c,t

def printSeqs(seqs):
    for s in seqs:
        q,c,t = s
        logging.debug( "".join(q))
        logging.debug( "".join(c))
        logging.debug( "".join(t))
    
def insert(base, pos, ch):
    base[0].insert(pos, ch)
    base[1].insert(pos, ch)
    base[2].insert(pos, ch)

def consensus(aligns):
    """
    expands alignment based on query, and then majority rules consensus the seqs
    """
    
    if len(aligns) > 500:#hard limit
        keep = []
        scores = map(lambda x: x.score, aligns)
        scores.sort()
        minS = scores[499]
        aligns = filter(lambda x: x.score <= minS, aligns)
        
    seqs = []
    for i in aligns:
        realign(i)
        seqs.append(offSetSeqs(i.tstart, i))
    
    #logging.debug("#Original Seqs (%d)" % (len(seqs)))
    #printSeqs(seqs)
    
    i = 0 #<-- target relative position
    remain = len(seqs) # number of sequences remaining
    while remain > 1:
        ins = False # Is there an insertion here
        for base in seqs:
            if i == len(base[1]):
                remain -= 1# kill it
            if  i < len(base[1]) and base[2][i] == '-':
                ins = True
        if ins: # apply insertion across all non-ins bases
            for base in seqs:
                if  i < len(base[1]):
                    if base[2][i] != '-' and base[2][i] != " ":
                        insert(base, i ,'_')
                    elif base[2][i] == " ":
                        insert(base, i, ' ')
        i += 1
    
    #logging.debug( "#Expanded Seqs" )
    #printSeqs(seqs)
    
    #majority vote consensus
    out = []
    contribBases = 0
    fillBases = 0
    if len(seqs) == 0:
        logging.info("no sequences")
    for p in range(max(map(lambda x: len(x[0]), seqs))):
        cnt = defaultdict(int)
        #Count it
        for s in seqs:
            if p < len(s[0]):
                cnt[s[0][p]] += 1
        cnt[" "] = 0
        contribBases += sum(cnt.values())
        #Maximum count
        nuc = max(cnt.iteritems(), key=operator.itemgetter(1))[1]
        #Who all has maximum count
        n = []
        for j in cnt:
            if cnt[j] == nuc:
                n.append(j)
        #get random one
        n = random.sample(n, 1)[0]
        if n not in ["-","_"]:
            fillBases += 1
        out.append(n)
    
    consen = "".join(out)
    logging.debug("# expanded consensus (%d nuc votes) <%d fill bases>" % (contribBases, fillBases))
    #logging.debug(consen)
    consen = consen.replace('_','').replace('-','').replace(' ','')

    results = namedtuple("polish_results", "contribSeqs contribBases fillBases sequence")
    return results(len(aligns), contribBases, fillBases, consen)

def parseArgs():
    parser = argparse.ArgumentParser(description=USAGE, \
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("reads", metavar="reads", type=str, \
                        help="Input reads .fasta or .fastq")
    
    parser.add_argument("-t", "--target", type=str, \
                        help="Target sequence name")
    parser.add_argument("-T", "--Target", type=str, \
                        help="Fasta file containing target sequence")
    parser.add_argument("-s", "--super", dest="super", action="store_true",\
                        help="Treat each read as the target once")
    
    parser.add_argument("-m", "--maxtail", type = int, default=sys.maxint, \
                        help="Max number of bases allowed to be in tail (inf)")
    
    parser.add_argument("-n", "--nproc", dest="nproc", default=1, type=int,\
                        help="Number of processors to use with blasr (1)")
    
    parser.add_argument("-o", "--outname", dest="outname", default="polish.out", \
                        type=str, \
                        help="Base name for output files (polish.out)")
    parser.add_argument("--debug", action="store_true")

    args = parser.parse_args()
    setupLogging(args.debug)
    
    #I don't think this is exhaustive
    if (args.target is not None and args.Target is not None) \
       or (args.super and (args.target is not None or args.Target is not None)):
        print "Error! only specify one of --super or --target or --Target"
        exit(1)
    
    return args

class NullDevice():
    def write(self, s):
            pass

if __name__ == '__main__':
    args = parseArgs()
    
    alignFile = args.outname+".m5"
    consensusFile = args.outname+".fasta"
    

    #extract the read I'm looking for    
    if args.target is not None:#Name
        tempOut = open("temp.fasta",'w')
        fasta = FastaFile(args.reads)
        tempOut.write(">%s\n%s\n" % (args.target, fasta[args.target]))
        tempOut.write
        blasr(args.reads, tempOut.name, nproc=args.nproc, outName=alignFile)
        
        aligns = M5File(alignFile)   
        fout = open(consensusFile, 'w')
        results = consensus(aligns)
        fout.write(">pbjpolish_%d_vote_%d_len\n" % (results.contribBases,\
                                     results.fillBases, results.sequence))
        #fout.write(">\n%s\n" % consensus(aligns))
    
        fout.close()    
    elif args.Target is not None:#File
        blasr(args.reads, args.Target, nproc=args.nproc, outName=alignFile)
        
        aligns = M5File(alignFile)   
        fout = open(consensusFile, 'w')
        results = consensus(aligns)
        fout.write(">pbjpolish_%d_vote_%d_len\n%s\n" % (results.contribBases,\
                                     results.fillBases, results.sequence))
        #fout.write(">%s\n%s\n" % consensus(aligns))
        fout.close()   
    elif args.super:#All
        tempfile = open("temp.fasta",'w')
        if args.reads.endswith(".fasta"):
            seqs = FastaFile(args.reads)
            #temp flie
            for s in seqs:
                tempfile.write(">%s\n%s\n" % (s, seqs[s]))
        elif args.reads.endswith(".fastq"):
            seqs = FastqFile(args.reads)
            #temp file
            for s in seqs:
                tempfile.write(">%s\n%s\n" % (s, seqs[s].seq))
        blasr(args.reads, tempfile.name, nproc=args.nproc, bestn=len(seqs), outName=alignFile)
        aligns = M5File(alignFile)   
        groups = defaultdict(list)
        for a in aligns:
            groups[a.tname].append(a)
        fout = open(consensusFile, 'w')
        for g in groups:
            results = consensus(aligns)
            fout.write(">pbjpolish_%d_vote_%d_len\n" % (results.contribBases,\
                                     results.fillBases, results.sequence))
        fout.close()
    
