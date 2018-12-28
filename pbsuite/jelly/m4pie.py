#!/usr/bin/env python
import os, re, argparse, logging, tempfile, sys
from collections import defaultdict

from pbsuite.utils.CommandRunner import exe
from pbsuite.utils.FileHandlers import revComp, M4File, FastaFile, FastqFile, M5File
from pbsuite.utils.setupLogging import setupLogging


USAGE="""\
Extracts softclip bases from aligned reads and remaps them to the provided reference. Produces a unified bam with reads containing updated information about tail-mapping.

WARNING! -- Input.bam should be produced without -noSplitSubreads in blasr (before blasr was fixed in smrtanalysis 2.1)
"""
    
def extractTails(aligns, reads, outFq, minLength=100):
    """
    0x1  -- template has multiple segments in sequencing
    0x40 -- first segment in template
    0x80 -- last segment in template

    Tail names will get _[pe][01].*:\d+
    on the end to hold metadata of:
    _    -- A delimieter
    [01] -- Strand of primary hit
    [pe] -- either p for prolog or e for epilog (strand ignorant)
    .*   -- refName of primary hit
    \d+  -- refPos of primary hit
    prolog, epilog, and the position of its primary (chr:pos)

    returns - nreads processed, ntails found, nreads with two ended tails
    """
    fout = open(outFq,'w')
    nreads      = 0
    ntails      = 0
    nmultitails = 0
    for r in reads.keys(): # protecting for spaces
        reads[r.split(' ')[0]] = reads[r]
        
    for read in aligns:
        nreads += 1
        pTail = read.qstart
        eTail = read.qseqlength - read.qend
        mateplace = read.tname
        
        strand = 1 if read.qstrand == '1' else 0
        
        hasTail = False
        len   = read.qseqlength
        if pTail >= minLength:
            hasTail = True
            ntails += 1
            seq   = reads[read.qname][:pTail]
            shift = 0
            fout.write(">%s_%d%s%d\n%s\n" % (read.qname, \
                       shift, 'p', len, seq))
             
        if eTail >= minLength:
            if hasTail:
                nmultitails += 1
            ntails += 1
            seq   = reads[read.qname][-eTail:]
            shift = read.qend
            fout.write(">%s_%d%s%d\n%s\n" % (read.qname, \
                       shift, 'e', len, seq))
        
    fout.close()
    return nreads, ntails, nmultitails
    
def mapTails(fq, ref, nproc=1, out="tailmap.sam", useSa=True):
    """
    automatically search for .sa
    """
    if os.path.exists(ref+".sa") and useSa:
        sa = "-sa " + ref + ".sa"
    else:
        sa = ""
    cmd = ("blasr %s %s %s -nproc %d -m 4 -bestn 1 -nCandidates 20 -out %s"
           " -minPctIdentity 75 -sdpTupleSize 6 -noSplitSubreads") \
           % (fq, ref, sa, nproc, out)
    
    logging.debug(cmd)
    r,o,e = exe(cmd)
    if r != 0:
        logging.error("blasr mapping failed!")
        logging.error("RETCODE %d" % (r))
        logging.error("STDOUT %s" % (str(o)))
        logging.error("STDERR %s" % (str(e)))
        logging.error("Exiting")
        exit(r)
    
    logging.info(str([r, o, e]))

def uniteTails(origAligns, tailMapFn, outMap="multi.m4", inplace=False):
    """
    Put the tails and original reads into a single m4.
    Add tags uniting the pieces

    every read comprises upto three pieces
        X->Y->Z
    or
        prolog->primary->epilog
    
    each piece has 3 tags added: (R) ref - (P) pos - (S) strand
    
    prolog and eplog will only point to the primary and the primary will point to both
    """
    datGrab = re.compile("^(?P<rn>.*)_(?P<shift>\d+)(?P<log>[pe])(?P<length>\d+)$")
    
    aligns = M4File(tailMapFn)
    mode = 'a' if inplace else 'w'
    aout = open(outMap, mode)
    
    nmapped = 0
    for read in aligns:
        nmapped += 1
        data = datGrab.search(read.qname).groupdict()
        read.qname = data["rn"]
        read.qseqlength = data["length"]
        read.qstart += int(data["shift"])
        read.qend += int(data["shift"])
        
        aout.write(str(read)+'\n')
    
    #consolidate information about the primary hits
    if not inplace:
        aout.write("\n".join([str(x) for x in origAligns]))
    
    aout.close()
    return nmapped

def parseArgs(argv):
    parser = argparse.ArgumentParser(description=USAGE, \
            formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("m4", metavar="M4", type=str, \
                        help="M4 containing mapped reads' alignments")
    parser.add_argument("reads", metavar="READS", type=str,\
                        help="Fasta/Fastq containing reads' sequence")
    parser.add_argument("ref", metavar="REFERENCE", type=str,\
                        help="REFERENCE to map tails to")
    parser.add_argument("-t", "--minTail", type=int, default=100,\
                        help="Minimum tail length to attempt remapping (100)")
    parser.add_argument("-n", "--nproc", type=int, default=1,\
                        help="Number of processors to use (1)")
    parser.add_argument("-o", "--output", type=str, default=None, \
                        help="Output Name (M4.tails.m4)")
    parser.add_argument("-i", "--inplace", action="store_true", \
                        help="Append the results to the input m4 file. Overrules --output")
    parser.add_argument("--noSa", action="store_true", \
                        help="Don't use reference's sa")
    parser.add_argument("--temp", type=str, default=tempfile.gettempdir(),
                        help="Where to save temporary files")
    parser.add_argument("--debug", action="store_true")
    
    args = parser.parse_args(argv)
    if args.inplace:
        args.output = args.m4
    elif args.output is None:
        args.output = args.m4[:-3] + ".tails.m4"
    
    setupLogging(args.debug)
    return args
    
def run(argv):
    print argv
    args = parseArgs(argv)
    if args.m4.endswith("m5"):
        aligns = M5File(args.m4)
    else:
        aligns = M4File(args.m4)
    if args.reads.endswith("fasta"):
        reads = FastaFile(args.reads)
    elif args.reads.endswith("fastq"):
        temp = FastqFile(args.reads)
        reads = {}
        for i in temp:
            reads[i] = temp[i].seq
        del(temp)
    else:
        logging.error("Expected Fasta or Fastq for READS (%s)" % args.reads)
        exit(1)
    
    logging.info("Extracting tails")
    tailfastq = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, dir=args.temp)
    tailfastq.close(); tailfastq = tailfastq.name
    logging.debug("Tail read tmp file %s " % (tailfastq))
    r, t, m = extractTails(aligns, reads, outFq=tailfastq, minLength=args.minTail)
    
    logging.info("Parsed %d reads" % (r))
    logging.info("Found %d tails" % (t))
    logging.info("%d reads had double tails" % (m))
    if t == 0:
        logging.info("No tails -- Exiting")
        exit(0)
    
    logging.info("Mapping Tails")
    tailmap = tempfile.NamedTemporaryFile(suffix=".m4", delete=False, dir=args.temp)
    tailmap.close(); tailmap = tailmap.name
    logging.debug("Read map tmp file %s " % (tailmap))
    mapTails(tailfastq, args.ref, nproc=args.nproc, out=tailmap, useSa=args.noSa)
    
    logging.info("Consolidating alignments")
    logging.debug("Final file %s " % (args.output))
    n = uniteTails(aligns, tailmap, args.output, args.inplace)
    logging.info("%d tails mapped" % (n))
    
if __name__ == '__main__':
    run(sys.argv[1:])
