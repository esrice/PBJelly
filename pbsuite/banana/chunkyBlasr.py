#!/usr/env/bin python

import sys, os, glob, subprocess, logging
from string import Template
from optparse import OptionParser
    

#Edit the template below to customize the submission for your cluster.
clusterTemplate = Template("echo '${CMD}' | msub -N \"${JOBNAME}\" -q analysis -o ${STDOUT} -e ${STDERR} -l nodes=1:ppn=8,mem=48000mb")


#These are your default blasr parameters. Adjust at will.
parameters = "-maxScore -1000 -bestn 24 -maxLCPLength 16 -nCandidates 24 -noSplitSubreads"


command = Template("blasr ${FAS} ${REF} ${SA} -m 4 -out ${OUT} -start ${START} -stride ${STRIDE} ${EXTRAPARAMS}")

USAGE="""%prog <reads.fasta> --output <outName> [ <options> ]
Creates All vs All alignment of reads.fasta
This script builds and submits cluster commands that will perform blasr mapping.
Protip - Don't play around with this script. It will submit many jobs very quickly.
         Be cautious and careful when thinking about running this.
Number of jobs created =
    (nChunks^2 + nChunks)
    -------------------- * stride
             2 """

def _exe(cmd):
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdoutVal, stderrVal =  proc.communicate()
    retCode = proc.returncode
    return retCode,stdoutVal,stderrVal
    
class ChunkyBlasr():
    
    def __init__(self, args):
        self.parseArgs(args)
        self.initLog()
    
    def parseArgs(self, argvs):
        
        parser = OptionParser(USAGE) 
        parser.add_option("--output", dest="output", \
                help=("Output results basepath - creates 'refChunks' "
                    "and 'aligns' subdirs"), default=None)
        parser.add_option("--stride", dest="stride", type="int", \
                help="See -stride and -start in blasr -h [DEFAULT=1]", default = 1)
        parser.add_option("--nChunks", dest="nChunks", type="int", default = 1, \
                help="Chunk reads into Number of Pieces [DEFAULT=1]")
        
        #parser.add_option("--chunkSize", dest="chunkSize", type="int", \
                #help=("Chunk Reference into Pieces of Size (in megabytes)"
                    #"overrides --nChunks"), default = None)
        
        parser.add_option("-f", "--filter", type="int", default=500,
                help="Filter reads < FILTER bp in length [DEFAULT=500]")
        parser.add_option("--sa", action="store_true",
                help="Create .sa indices for the chunks")
        parser.add_option("--skipSplit", action="store_true",
                help=("Skip Splitting Reference (only skipSplit if you've "\
                    " split in the past and are redoing your alignments)"))
        parser.add_option("-p", "--params", type="str", \
                help=("Parameters to pass to blasr. Surround string of " 
                    "params with \"'s"), default=parameters)
        parser.add_option("--debug", action="store_true",
                help="Verbose logging output")
        
        self.opts, args = parser.parse_args(argvs) 
        args = args[1:]#..stupid program..
        if len(args) != 1:
            parser.error("Error! Expceted 1 argument - the input reads")
        
        self.reads = os.path.abspath(args[0])
        
        
        if self.opts.output is None or not os.path.isdir(self.opts.output):
            parser.error("Error! Must specify output directory")
        
        self.opts.output = os.path.abspath(self.opts.output)
    
    def initLog(self):
        logLevel = logging.DEBUG if self.opts.debug else logging.INFO
        logFormat = "%(asctime)s [%(levelname)s] %(message)s"
        logging.basicConfig( stream=sys.stderr, level=logLevel, format=logFormat )
        logging.info("Running %s" % " ".join(sys.argv) )
    
    
    def setupPaths(self):
        try:
            self.refDir = os.path.join(self.opts.output, "refChunks")
            os.mkdir(self.refDir)
        except OSError:
            logging.warning("Reference Directory %s already exists" % self.refDir)
        
        try:
            self.alignDir = os.path.join(self.opts.output, "aligns")
            os.mkdir(self.alignDir)
        except OSError:
            logging.warning("Alignment Directory %s already exists" % self.alignDir)
        
        self.refBase = os.path.basename(self.reads)
        
    def chunkReference(self):
        """
        Chunks a readgroup into x pieces
        """
        nChunks = self.opts.nChunks
        
        logging.info("Chunking Reference into %d pieces" % (nChunks))
        r, o, e = _exe("fastasplit --fasta %s --output %s --chunk %d " % \
                    (self.reads, self.refDir, nChunks))
        logging.debug("RETCODE - %d\nSTDOUT - %s\nSTDERR - %s" %(r,str(o),str(e)))
        if r != 0:
            logging.error("Chunking Failed")
            exit(r)
        
        logging.info("Finished Chunking Reference")
        chunks = []
        logging.info("Filtering Reads < %d")
        for chunk in glob.glob(os.path.join(self.refDir, self.refBase+"*_chunk_*")):
            logging.info("Filtering out reads < 500bp from %s" % chunk)
            r, o, e = _exe(("fastalength {0} | "
                            "awk -F\  '{{if ($1 < 500) print $2}}' | "
                            "fastaremove {0} stdin > {0}.fasta").format(chunk))
            logging.debug("RETCODE - %d\nSTDOUT - %s\nSTDERR - %s" %(r,str(o),str(e)))
            if r != 0:
                logging.error("Filtering %s Failed" % (chunk))
                exit(r)
            
            logging.info("Removing chunk %s." % chunk)
            os.remove(chunk)
            #Filtering renamed it
            chunk = chunk+".fasta"
            
            #I'm worried that calling sawriter this way is breaking index
            if self.opts.sa:
                logging.info("Indexing %s" % (chunk))
                r,o,e = _exe("sawriter {0}.sa {0}".format(chunk))
                logging.debug("RETCODE - %d\nSTDOUT - %s\nSTDERR - %s" %(r,str(o),str(e)))
                if r != 0:
                    logging.error("Indexing %s Failed" % (chunk))
                    exit(r)
                    
                logging.info("Finished Indexing %s" % (chunk))
            
            chunks.append(chunk)
        
        return chunks
    
    def run(self):
            
        self.setupPaths()
        if self.opts.skipSplit:
            logging.info("Skipping Splitting")
            chunks = glob.glob(os.path.join(self.refDir, self.refBase+"*_chunk_*.fasta"))
        else:
            chunks = self.chunkReference()
        
        jobsSubmitted = 0
        for c1, ref in enumerate(chunks):
            for c2, reads in enumerate(chunks[c1:]):
                c2 += c1
                for i in range(int(self.opts.stride)):
                    chunkName = "chunk_%d.chunk_%d.stride_%d.m4" % (c1, c2, i)
                    myOutFile = os.path.join(self.alignDir, chunkName)
                    if os.path.exists(ref+".sa"):
                        saIdx = "-sa " + ref + ".sa"
                    else:
                        saIdx = ""
                    myParams = {"REF":ref, "SA": saIdx, "FAS": reads, \
                                "OUT":myOutFile, "START":i, "STRIDE":self.opts.stride, \
                                "EXTRAPARAMS":self.opts.params}
                    
                    #Build the stuff for the cluster Command
                    myCommand = {"CMD":command.substitute(myParams), \
                                "JOBNAME":chunkName, \
                                "STDOUT":myOutFile+".out", \
                                "STDERR":myOutFile+".err" } 
                    
                    #Submit the cluster command
                    logging.info("Submitting Chunk %d vs Chunk %d Stride %d" % (c1, c2, i))
                    runCmd = clusterTemplate.substitute(myCommand)
                    logging.debug("CMD: " + runCmd)
                    logging.debug(str(_exe(runCmd)))
                    jobsSubmitted += 1
        logging.info("Finished Running chunkyBlasr.py -- %d Jobs Submitted" % \
                                                                (jobsSubmitted))

if __name__ == '__main__':
    me = ChunkyBlasr(sys.argv)
    me.run()
