#!/usr/bin/env python
import sys, os, tempfile, logging, subprocess, shutil, re, signal
from optparse import OptionParser
from StringIO import StringIO

from pbsuite.utils.setupLogging import *
from pbsuite.utils.FileHandlers import *

import pbpy.io.AmosBank as amos

USAGE = """ Usage: %prog [<input.fasta> <input.qual> | <input.fastq>] [--options]

This is a wrapper around the de novo assembly process
used by PacificBioscience's SMRTAnalysis software.
"""
"""
Need more argument passing transparcency to wrapped PacBio scripts

>>> ord('O')
79 -- Exit Code for No Overlaps
>>> ord('L')
76 -- Exit Code for No Layout
>>> ord('C') -- May Not Every Actually Happen
67 -- Exit Code for No Contigs
"""
class FastQ():
    #Replacement for namedtuples
    def __init__(self, name, seq, qual):
        self.name = name
        self.seq = seq
        self.qual = qual
    
    def __cmp__(self, o):
        return len(self.seq) - len(o.seq)

def _exe( cmd , wait=True):
        log = StringIO()
        log.write("Running %s\n" % cmd)
        
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, \
                                stderr=subprocess.STDOUT, close_fds=True)
        if not wait:
            return proc#you'll have to communicate with it
        stdoutVal, stderrVal =  proc.communicate()
        retCode = proc.returncode
        
        
        log.write("Return Code = %s\n" % retCode)
        log.write("STDOUT - %s\n" % stdoutVal)
        log.write("STDERR - %s\n" % stderrVal)
        log.flush()
        return log.getvalue()

class OLCAssembly:
    def __init__(self):
        self._parseOptions()
        
    def _parseOptions( self ):
        parser = OptionParser( usage=USAGE )
        
        parser.add_option("--debug", action="store_true", help="Increases verbosity of logging" )
        parser.add_option("--nproc", type="int", help="Number of processes to use." )
        parser.add_option("-o", "--outName", type="string", help="Name of the output fasta and qual files (Don't include the extension)", default="out")
        parser.add_option("--fqOut", action="store_true", help="Create a .fastq output file")
        parser.add_option("--rename", type="string", help="Gives the ouput contigs more descriptive names")
        parser.add_option("--minSubreads", type="int", help="Minimum number of subreads required to attempt assembly")
        parser.add_option("--workDir", type="string", help="Directory to build the bank an everything in.")
        parser.add_option("--workTmp", type="string", help="Work in a temporary directory")
        parser.add_option("--threshold", type="int", help="Threshold when determining overlaps")
        parser.add_option("--transmax", type="int", help="Max links of transitivity")
        parser.add_option("-e", type="str", help="Alignment Error% e.g. 0.15 = 15%")
        
        parser.set_defaults(debug=False, nproc=1, outName="out", rename=None, minSubreads=2, \
            filtering=False, workTmp=None, threshold=800, transmax=1, e="0.15")
        
        self.options, args = parser.parse_args(sys.argv)
        setupLogging(self.options.debug)
        logging.warning("This program doesn't work with SMRTAnalysis v2.1 and on")
        logging.info("Reading Input Reads")
        if len(args) == 2:
            self.fastqFile = args[1]
            if not self.fastqFile.endswith(".fastq"):
                parser.error("Expected a Fastq File or Fasta/Qual")
            self.fastqSeq = FastqFile(self.fastqFile)
        elif len(args) == 3:
            fasta = args[1]
            qual = args[2]
            if not fasta.endswith(".fasta"):
                parser.error("Expected First Argument To End With .fasta")
            if not qual.endswith(".qual"):
                parser.error("Expected Second Argument To End With .qual")
            self.fastqSeq = mergeFastaQual(fasta, qual)
            self.fastaFile = fasta
            self.qualFile = qual
        else:
            parser.error("Expected <input.fastq> or <input.fasta> <input.qual> Arguments!")
        
        self.options.outName = os.path.abspath(self.options.outName)
        
        if self.options.workTmp is not None:
            self.options.workDir = tempfile.mkdtemp(dir=self.options.workTmp)
        if self.options.workDir is not None:
            os.chdir(self.options.workDir)
        
    def loadSequence(self):
        """
        populates what used to be pooledSubreads
        """
        
        self.fasta = FastaFile(self.fastaFile)
        
        logging.info("Reading Qual")
        self.qual = QualFile(self.qualFile)
        
        logging.info("Creating FastQ")
        self.sequences = []
        
        for key in self.fasta.keys():
            self.sequences.append(FastQ(key, self.fasta[key], self.qualToPhred(self.qual[key]) ) )
        
        #Sort Because blas's inconsistency. -- They're fixing it. Take this out when the fix is released
        #It is fixed but isn't part of SMRTPortal, yet -- so keep this in for consistency.
        #Note: This doesn't actually fix the problem, it just makes for reproduceable results
        #self.sequences.sort()
        
        logging.info("Total %d Reads to Process" % len(self.sequences))
        if len(self.sequences) < self.options.minSubreads:
            logging.info("Not enough for assembly, skipping")
            if self.options.rename is not None:
                fastOut = open(self.options.outName+".fasta", 'w')
                qualOut = open(self.options.outName+".qual", 'w')
                for i,key in enumerate(self.fasta.keys()):
                    fastOut.write(">%s_%d\n%s\n" % (self.options.rename, i+1, self.fasta[key]))
                    qualOut.write(">%s_%s\n%s\n" % (self.options.rename, i+1, " ".join(map(str,self.qual[key]))))
                fastOut.close()
                qualOut.close()
            else:
                logging.debug(_exe("cp %s %s" % (self.fastaFile, self.options.outName+".fasta")))
                logging.debug(_exe("cp %s %s" % (self.qualFile, self.options.outName+".qual")))
            exit(0)
    
    def setup(self):
        """
        Makes a Temporary Directory with the files we're going to need. Returns the directory's name
        Make Temporary Directory to work in.
        Make temp.fasta & temp.fastq for pipeline
        """
        fq = open("inputReads.fastq",'w')
        fa = open("reference.fasta",'w')
        for key in self.fastqSeq.keys():
            fq.write(str(self.fastqSeq[key]))
            fa.write(">%s\n%s\n" % (key, self.fastqSeq[key].seq))
        fa.close()
        fq.close()
        
    def mapReads(self):
        """
        Takes all of the pools and generates their fastq and alignments in their own folder. 
        """
        logging.info("Creating Overlap")
        log = _exe(("blasr inputReads.fastq reference.fasta -nproc %d -m 4 " 
                    "-out temp.rm4 -noSplitSubreads -useGuidedAlign -allowAdjacentIndels "#-minFrac 0.01 
                    "-nCandidates 20 -bestn 15 -minMatch 8 -maxLCPLength 16 ") % (self.options.nproc) )
        logging.debug(log)
        logging.info("Sorting the alignments")
        logging.debug(_exe("sort temp.rm4 > alignments.rm4"))

    def runAssembly(self):
        """
        Setup a bank in the given path with all of the filtered subreads
        """
        logging.info("Running Assembly")
        
        logging.info("Creating afg file")
        logging.debug(_exe("toAfg inputReads.fastq out.afg -noSplitSubreads"))
        
        logging.info("Creating out.bank")
        logging.debug(_exe("bank-transact -f -c -b out.bank -m out.afg"))
        #"""
        logging.info("Finding Overlaps")
        logging.debug(_exe(("transitiveOverlap.py --nproc=%d --assignIIDs "
                            "--threshold %d  --transmax %d --align zscore "
                            "alignments.rm4 out.bank --requireFullOverlap > "
                            "alignments.overlaps_realigned ") \
                            % (self.options.nproc, self.options.threshold, \
                               self.options.transmax)))
            

        #Sometimes one read is contained in another
        if os.path.getsize("alignments.overlaps_realigned") == -1:
            logging.info("Removing full overlap requirement and retrying")
            logging.debug(_exe(("transitiveOverlap.py --nproc=%d --assignIIDs "
                                "--threshold %d  --transmax %d --align zscore "
                                "alignments.rm4 out.bank > "
                                "alignments.overlaps_realigned ") 
                                % (self.options.nproc, self.options.threshold, \
                                self.options.transmax)))
           
        if os.path.getsize("alignments.overlaps_realigned") == 0:
            logging.error("No Overlaps. Exiting")
            _exe("touch %s && touch %s" % (self.options.outName+".fasta", self.options.outName+".qual"))
            exit(79)
        
        logging.info("Loading overlaps into out.bank")
        logging.debug(_exe("bank-transact -b out.bank -m alignments.overlaps_realigned"))

        logging.info("Constructing Layout")
        logging.debug(_exe("clusterOverlapIIDs.py --nproc=%d out.bank > __splits__" % (self.options.nproc) ))

        logging.info("Creating contigs")
        logging.debug(_exe("tigger -i 0_iids.split -b out.bank > layout.out"))

        if os.path.getsize("layout.out") == 0:
            logging.error("No Layout. Exiting")
            _exe("touch %s && touch %s" % (self.options.outName+".fasta", self.options.outName+".qual"))
            exit(76)
        
        logging.info("Loading layout into out.bank")
        logging.debug(_exe("bank-transact -b out.bank -m layout.out"))
        
        #"""
        #New Method...  Only works in Pre_Assembly correcting reads.
        """
        logging.info("Turning Blasr Alignments to Layout")
        logging.debug(_exe("align2layouts.py alignments.rm4 out.bank --overlapTolerance 100 > alignments.lay"))
        
        logging.info("Loading layout into out.bank")
        logging.debug(_exe("bank-transact -b out.bank -m alignments.lay"))
        
        #Here is where I can put multiprocessing on consensus building
        logging.info("Getting Layout IIDs")
        logging.debug(_exe("grep iid alignments.lay | cut -f2 -d: > layout.ids"))
       
        if os.path.getsize("layout.ids") == 0:
            logging.error("No Layout. Exiting")
            _exe("touch %s && touch %s" % (self.options.outName+".fasta", self.options.outName+".qual"))
            exit(76)
        """
        logging.info("Building Consensus")
        
        #if self.options.nproc > 1:
        #couldn't get this to work. tigger makes layouts that suck and id shit
            #splits = self.splitIIDs(self.options.nproc)
            #consensusOut = []
            #procWait = []
            #for pos,split in enumerate(splits):
                #cout = "consensus%d.out" % (pos)
                #procWait.append(_exe("make-consensus -A -L -i %s -b out.bank > %s" \
                                     #% (split, cout), wait=False))
                #consensusOut.append(cout)
            #for i in procWait:
                #o,e = i.communicate()
            #logging.info("Gathering Consensus")
            #_exe("awk '{ if (/^iid:[0-9]+$/) {print \"iid:\"++iid} else {print}; }' %s > consensus.out" \
                    #% (" ".join(consensusOut)))
        #else:
        logging.debug(_exe("make-consensus -e %s -b out.bank -A -L > consensus.out" % (self.options.e)))
        
        logging.info("Loading consensus into out.bank")
        logging.debug(_exe("bank-transact -b out.bank -m consensus.out"))
        eid = ""    
        if self.options.rename is not None:
            logging.info("Renaming Contig Sequences")
            self.renameContigs(self.options.rename) 
            eid = "-eid"
        
        logging.info("Dumping Contig Sequence")
        logging.debug(_exe("bank2fasta %s  -b out.bank -q %s.qual > %s.fasta" % (eid, self.options.outName, self.options.outName)))
        
        if self.options.fqOut:
            fout = open(self.options.outName + ".fastq", 'w')
            for entry in mergeFastaQual(self.options.outName+".fasta", self.options.outName+".qual").values():
                fout.write(entry.toString())
            fout.close()
            
    def splitIIDs(self, m):
        _exe("grep iid alignments.overlaps_realigned | cut -f2 -d\: > internal.ids")
        fh = open("internal.ids",'r')
        n = fh.readlines()
        fh.close()
        p = map(lambda x: list(), range(m))
        index = 0
        for item in n:
            p[index].append(item)
            if index < m-1:
                index += 1
            else:
                index = 0
        ret = []
        for pos, i in enumerate(filter(lambda x: len(x)>0, p)):
            spName = "splitCon%d.iids" % (pos)
            ret.append(spName)
            fout = open(spName, 'w')
            fout.write("".join(i))
            fout.close()
        return ret
    
    def renameContigs(self, name):
        """
        Renames contigs based on what reads compose it.
        TODO: Need to know the exception that's thrown when there are no contigOuts
        """
        logging.info("Renaming Contigs")
        getContigName = re.compile("(\d+)_contig")
        contigNames = {}
        bank = amos.AmosBank("out.bank")
        
        for x in bank.iAlignmentHits():
            read = int(x.query_id)
            read = bank.getRead(read).getEID()
            contig = int(getContigName.match(x.target_id).groups()[0])
            contigNames[contig] = read.split('/')[0]
        
        contigBankStream = amos.BankStream_AMOS.BankStream_t("CTG")
        contigBankStream.open(bank.bankPath)
        
        for index,iid in enumerate(contigNames.keys()):
            contig = amos.Contig_AMOS.Contig_t()
            contigBankStream.fetch(iid, contig)
            newName = "%s_%d" % (name, index + 1) 
            contig.setEID(newName)
            contigBankStream.replace(iid, contig)
         
        contigBankStream.close()
    
    def run(self):
        """
        Given an input.fasta and input.fastq - run the Assembly Pipeline
        """
        self.setup()
        self.mapReads()
        self.runAssembly()
        
        logging.info("Finished Assembly")
    
    def makeFastq(self, fastq):
        """
        Given a fastq named tuple, return fastq format string
        """
        return "@%s\n%s\n+%s\n%s\n" % (fastq.name, fastq.seq, fastq.name, fastq.qual)
    
    def qualToPhred(self, quals):
        """
        Given list of Integers that represent base phred scores, return the phred encoding
        """
        return "".join(map(lambda x: chr(x+33), quals))
    
    def makeFasta(self, fastq):
        """
        Given a fastq named tuple, return fasta format string
        """
        return ">%s\n%s\n" % (fastq.name, fastq.seq)
    
    
if __name__ == '__main__':
    me = OLCAssembly()
    me.run()
