#!/usr/bin/env python
import os
import sys
import shutil
import random
import argparse
import tempfile
import multiprocessing
from collections import namedtuple

import pysam

from pbsuite.utils.BedIO import *
from pbsuite.utils.FileHandlers import *
from pbsuite.utils.setupLogging import *
from pbsuite.utils.CommandRunner import *

USAGE = """\
Takes a list of putative SVs and a set of BAMs builds
as many of your sites as possible.
Remap Files are stored in 
<--output>.sam
<--output>.tails.sam
"""

CRITICAL = """\
Need to separate arguments
Need to make it explicit that Celera is ran blandly
Only the first bam's insert size is considered
"""

class Consumer(multiprocessing.Process):
    
    def __init__(self, task_queue, result_queue, args):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.args = args

    def run(self):
        try: 
            proc_name = self.name
            #Open the bams
            nBams = [] # nonTrim Bams
            tBams = [] # trim Bams
            for i in self.args.bam:
                nBams.append(pysam.Samfile(i))
            for i in self.args.pacBam:
                tBams.append(pysam.Samfile(i))
            
            while True:
                next_task = self.task_queue.get()
                if next_task is None:
                    # Poison pill means shutdown
                    logging.info('Thread %s: Exiting\n' % proc_name)
                    self.task_queue.task_done()
                    break
                try:
                    answer = next_task(nBams, tBams)
                except Exception as e:
                    logging.error("Exception raised in task %s" % (str(e)))
                    self.task_queue.task_done()
                    self.result_queue.put("Failure - UNK - %s" % str(e))
                    logging.info("fail in groupid=%s" % next_task.data.name)
                    continue
                self.task_queue.task_done()
                self.result_queue.put(answer)
            return
        except Exception as e:
            logging.error("Consumer %s Died\nERROR: %s" % (self.name, e))
            return
            
class Assembler(object):
    
    def __init__(self, data, args):
        """
        Args is a namedtuple of every parameter needed
        required args are buffer, tmpDir, and timeout
        """
        #buffer, tmpDir, timeout):
        self.data = data
        self.args = args
        self.buffer = args.buffer
        self.tmpDir = args.temp
        self.timeout = args.timeout
    
    def fetchReads(self, bam, chrom, start, end, trim=False):
        """
        Trim is for pacbio reads that have tails that were attempted to 
        be remapped -- helps reduce redundancy
        Also, I'm going to read tails that extend beyond my boundaries
        """
        logging.info("fetching %s from %s:%d-%d" % (bam.filename, chrom, start, end))
        ret = {}
        
        for id, read in enumerate(bam.fetch(reference=chrom, start=start, end=end)):
            name = read.qname 
            name += " DIRECTION: rev" if read.is_reverse else " DIRECTION: fwd"
            seq, qual = read.seq, read.qual 
            
            if trim and not read.is_unmapped:
                regTrim = 0
                
                trimS = None
                trimE = None
                if start > read.pos:
                    for queryPos, targetPos in read.aligned_pairs:
                        if trimS is None and targetPos is not None and targetPos >= start:
                            trimS = queryPos
                if end < read.aend:
                    for queryPos, targetPos in read.aligned_pairs[::-1]:
                        if trimE is None and targetPos is not None and targetPos <= end:
                            trimE = queryPos
                
                if trimS is not None:
                    upS = read.cigar[0][1] if read.cigar[0][0] == 4 else 0
                    trimS = max(0, trimS) + upS
                else:
                    trimS = 0
                
                if trimE is not None:
                    dnS = read.cigar[-1][1] if read.cigar[-1][0] == 4 else 0
                    trimE = min(len(read.seq), trimE)  - dnS
                else:
                    trimE = len(read.seq)
                
                seq = seq[trimS:trimE]
                qual = qual[trimS:trimE]
            
            ret[name + seq[:10]] = [name, seq, toQual(qual)]
        logging.info("%d reads retreived" % len(ret))
        return ret
    
    def cleanupTmp(self):
        for i in self.myTmpFiles:
            if os.path.exists(i):
                if os.path.isfile(i):
                    os.remove(i)
                else:
                    shutil.rmtree(i)
    
    def __str__(self):
        return self.result

class PhrapAssembler(Assembler):
    
    def __init__(self, data, args):
        #buffer, tmpDir, timeout, *args, **kwargs):
        Assembler.__init__(self, data, args)
    
    def __assemble(self, reads):
        """
        writes temp files
        assembles
        reads results
        clears temp files
        returns results as a string
        Calls the assembler
        """
        self.myTmpFiles = []
        #Temporary Files
        fout = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, dir=self.tmpDir)
        self.myTmpFiles.append(fout.name)
        qout = open(fout.name + '.qual', 'w')
        self.myTmpFiles.append(fout.name + '.qual')
        
        for name, seq, qual in reads:
            fout.write(">%s\n%s\n" % (name, seq))
            qout.write(">%s\n%s\n" % (name, qual))
        
        fout.close()
        qout.close()

        r, o, e = exe("phrap %s -minmatch 6 -minscore 20" % (fout.name),\
                      timeout=self.timeout)
           
        self.myTmpFiles.extend([fout.name + ".contigs",  fout.name + ".contigs.qual", \
                           fout.name + ".problems", fout.name + ".problems.qual", \
                           fout.name + ".log",      fout.name + ".singlets"])
        if r == 214:
            super(PhrapAssembler, self).cleanupTmp()
            return "Failure - Assembly Timeout " + self.data.name
         
        results = mergeFastaQual(fout.name + ".contigs", fout.name + ".contigs.qual")
        
        #Try to push the problems through, too
        if os.stat(fout.name + '.problems').st_size != 0:
            pfile = fout.name + ".problems"
            r, o, e = exe("phrap %s -minmatch 6 -minscore 20" % (pfile), \
                        timeout=self.timeout)
                    
            self.myTmpFiles.extend([pfile + ".contigs",  pfile + ".contigs.qual", \
                            pfile + ".problems", pfile + ".problems.qual", \
                            pfile + ".log",      pfile + ".singlets"])
            if r == 214:
                super(PhrapAssembler, self).cleanupTmp()
                return "Failure - Assembly Timeout " + self.data.name
            
            results.update(mergeFastaQual(fout.name + ".problems.contigs", fout.name + ".problems.contigs.qual"))
        
        #save to file
        fout = tempfile.NamedTemporaryFile(prefix = "asm" + self.data.name, \
                                    suffix=".fastq", delete=False, dir=self.tmpDir)
        for key in results:
            fout.write("@group" + self.data.name + "_" + key + "\n" + \
                        results[key].seq + '\n+\n' + \
                        results[key].qual + '\n')
        fout.close()
        self.results = fout.name
        
        #clean up
        super(PhrapAssembler, self).cleanupTmp()
        
        return self.results
    
    def __call__(self, nBams, tBams):
        #Fetch,
        logging.info("asm task groupid=%s start" % (self.data.name))
        reads = {}
        chrom = self.data.chrom
        start = self.data.start - self.buffer
        start = max(0, start)
        end = self.data.end + self.buffer
        
        for bam in nBams:
            if self.data.start + self.buffer >= self.data.end - self.buffer:
                reads.update(super(PhrapAssembler, self).fetchReads(bam, chrom, start, end))
            else:
                reads.update(super(PhrapAssembler, self).fetchReads(bam, chrom, \
                             max(0, self.data.start - self.buffer), self.data.start + self.buffer))
                reads.update(super(PhrapAssembler, self).fetchReads(bam, chrom, \
                             max(0, self.data.end - self.buffer), self.data.end + self.buffer))
                
        if len(reads) > self.args.maxreads:
            logging.info("Downsampling %s" % (self.data.name))
            nreads = {}
            for i in random.sample(reads.keys(), self.args.maxreads):
                nreads[i] = reads[i]
            reads = nreads
        
        for bam in tBams:
            if self.data.start + self.buffer >= self.data.end - self.buffer:
                reads.update(super(PhrapAssembler, self).fetchReads(bam, chrom, start, end, trim=True))
            else:
                reads.update(super(PhrapAssembler, self).fetchReads(bam, chrom, \
                        max(0, self.data.start - self.buffer), self.data.start + \
                        self.buffer, trim=True))
                reads.update(super(PhrapAssembler, self).fetchReads(bam, chrom, \
                             max(0, self.data.end - self.buffer), self.data.end +
                             self.buffer, trim=True))
             
        reads = reads.values() 
        totReads = len(reads)
        
        #Assemble
        logging.info("assembling %d reads" % (len(reads)))
        self.result = self.__assemble(reads)
        logging.info("asm task groupid=%s finish" % (self.data.name))
        return self.result
        
        #I can eventually save memory by just returning the output file's name - 
        #Though to do that I 'd need to edit the file with the reads' groupname
        
class MiniaAssembler(Assembler):
    
    def __init__(self, data, args):
        #buffer, tmpDir, timeout, *args, **kwargs):
        Assembler.__init__(self, data, args)
            
    def __assemble(self, reads):
        """
        writes temp files
        assembles
        reads results
        clears temp files
        returns results as a string
        Calls the assembler
        """
        self.myTmpFiles = []
        #Temporary Files
        fout = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, dir=self.tmpDir)
        self.myTmpFiles.append(fout.name)
        
        for name, seq, qual in reads:
            fout.write(">%s\n%s\n" % (name, seq))
        
        fout.close()
        resultOut = tempfile.NamedTemporaryFile(prefix="minia_", delete=False, dir=self.tmpDir)
        
        estSize = self.buffer * 2
        if self.data.rest[0] != 'DEL':
            estSize += int(self.data.rest[1])
        
        r, o, e = exe("minia {reads} 20 3 {estSize} {outPrefix}".format(reads=fout.name, \
                      estSize=estSize, outPrefix=resultOut.name ), \
                      timeout=self.timeout)
        
        logging.debug("RET - %d\nOUT - %s\nERR- %s" % (r, o, e))
        
        self.myTmpFiles.extend([resultOut.name + ".contigs.fa",  resultOut.name + ".debloom", \
                           resultOut.name + ".debloom2", resultOut.name + ".false_positive_kmers", \
                           resultOut.name + ".reads_binary",      resultOut.name + ".solid_kmers_binary"])
        
        if r == 214:
            super(MiniaAssembler, self).cleanupTmp()
            return "Failure - Assembly Timeout " + self.data.name
         
        fasta = FastaFile(resultOut.name + ".contigs.fa")
        
        results = {}
        for key in fasta:
            results[key] = FastqEntry(key, fasta[key], '?' * len(fasta[key]))
        
        #save to file
        fout = tempfile.NamedTemporaryFile(prefix = "asm" + self.data.name, \
                                    suffix=".fastq", delete=False, dir=self.tmpDir)
        for key in results:
            fout.write("@group" + self.data.name + "_" + key + "\n" + \
                        results[key].seq + '\n+\n' + \
                        results[key].qual + '\n')
            
        fout.close()
        self.results = fout.name
        
        #clean up
        super(MiniaAssembler, self).cleanupTmp()
        
        return self.results

    def __call__(self, nBams, tBams):
        #Fetch,
        logging.info("asm task groupid=%s start" % (self.data.name))
        reads = []
        chrom = self.data.chrom
        start = self.data.start - self.buffer
        start = max(0, start)
        #hope that fetching beyonde 3' boundary is okay
        end = self.data.end + self.buffer
        
        for bam in nBams:
            reads.extend(self.fetchReads(bam, chrom, start, end))
        if len(tBams) > 1:
            logging.warning("Minia isn't built to handle PacBio Reads.. results may be unoptimal")
        for bam in tBams:
            reads.extend(super(MiniaAssembler, self).fetchReads(bam, chrom, start, end, trim=True))
        
        if len(reads) > self.maxreads:
            return "Failure - Too Many Reads (%d) %s" % (len(reads), self.data.name)
        
        #Assemble
        logging.info("assembling %d reads" % (len(reads)))
        self.result = self.__assemble(reads)
        logging.info("asm task groupid=%s finish" % (self.data.name))
        return self.result
        
        #I can eventually save memory by just returning the output file's name - 
        #Though to do that I 'd need to edit the file with the reads' groupname
        
class SpadesAssembler(Assembler):
    
    def __init__(self, data, args):
        #buffer, tmpDir, timeout, threads, maxreads):
        super(SpadesAssembler, self).__init__(data, args)
        self.threads = args.nproc
    
    def __fetchPEReads(self, bam, chrom, start, end):
        def write_read(read):
            """
            Write read to open FASTQ file.
            """
            COMPLEMENT = {'A': 'T',
              'T': 'A',
              'C': 'G',
              'G': 'C',
              'N': 'N'}

            if read.is_reverse:
                ret = [read.qname + ("/%d" % (int(not read.is_read1) + 1)),
                       "".join((COMPLEMENT[b] for b in read.seq[::-1])),
                       read.qual[::-1]]
            else:
                ret = [read.qname + ("/%d" % (int(not read.is_read1) + 1)),
                       read.seq,
                       read.qual]
            return ret
        
        logging.info("fetching %s from %s:%d-%d" % (bam.filename, chrom, start, end))
        name = read_left = read_right = None
        
        sync_pairs = False
        for read in bam.fetch(reference=chrom, start = start, end = end):
            if name is not None and read.qname != name:
                if read_left and (not sync_pairs or read_right):
                    self.leftReads.append(write_read(read_left))
                if read_right and (not sync_pairs or read_left):
                    self.rightReads.append(write_read(read_right))
                read_left = read_right = None
            name = read.qname
            if read.is_read1:
                read_left = read
            else:
                read_right = read
        
        if read_left and (not sync_pairs or read_right):
            self.leftReads.append(write_read(read_left))
        if read_right and (not sync_pairs or read_left):
            self.rightReads.append(write_read(read_right))
        
    def __assemble(self):
        """
        writes temp files
        assembles
        reads results
        clears temp files
        returns results as a string
        Calls the assembler
        """
        self.myTmpFiles = []
        #Temporary Files
        fout = tempfile.NamedTemporaryFile(prefix="spades_pe1", suffix=".fastq", delete=False, dir=self.tmpDir)
        self.myTmpFiles.append(fout.name)
        for name, seq, qual in self.leftReads:
            fout.write("@%s\n%s\n+\n%s\n" % (name, seq, qual))
        fout.close()
        
        fout2 = tempfile.NamedTemporaryFile(prefix="spades_pe2", suffix=".fastq", delete=False, dir=self.tmpDir)
        self.myTmpFiles.append(fout2.name)
        for name, seq, qual in self.rightReads:
            fout2.write("@%s\n%s\n+\n%s\n" % (name, seq, qual))
        fout2.close()
        
        foutp = tempfile.NamedTemporaryFile(prefix="spades_pb", suffix=".fastq", delete=False, dir=self.tmpDir)
        self.myTmpFiles.append(foutp.name)
        for name, seq, qual in self.pbReads:
            foutp.write("@%s\n%s\n+\n%s\n" % (name, seq, qual))
        foutp.close()
        
        #working here
        resultOut = tempfile.mkdtemp(prefix="spades", dir=self.tmpDir)
        
        estSize = self.buffer * 2
        if self.data.rest[0] != 'DEL':
            estSize += int(self.data.rest[1])
        
        #r, o, e = exe("dipspades.py -1 {pe1} -2 {pe2} --pacbio {pacbio} -o {output} "\
        r, o, e = exe("spades.py -1 {pe1} -2 {pe2} --pacbio {pacbio} -o {output} "\
                      .format(pe1=fout.name, pe2=fout2.name, pacbio=foutp.name, output=resultOut), \
                      timeout=self.timeout)
                    
        logging.debug("RET - %d\nOUT - %s\nERR- %s" % (r, o, e))
        #just the output dir, maybe?
        self.myTmpFiles.append(resultOut)
        if r == 214:
            super(SpadesAssembler, self).cleanupTmp()
            return "Failure - Assembly Timeout " + self.data.name
         
        outFsta = os.path.join(resultOut, "dipspades", "consensus_contigs.fasta")
        fasta = FastaFile(outFsta)
        
        results = {}
        for key in fasta:
            results[key] = FastqEntry(key, fasta[key], '?' * len(fasta[key]))
        
        #save to file
        fout = tempfile.NamedTemporaryFile(prefix = "asm" + self.data.name, \
                                    suffix=".fastq", delete=False, dir=self.tmpDir)
        for key in results:
            fout.write("@group" + self.data.name + "_" + key + "\n" + \
                        results[key].seq + '\n+\n' + \
                        results[key].qual + '\n')
            
        fout.close()
        self.results = fout.name
        
        #clean up
        super(SpadesAssembler, self).cleanupTmp()
        
        return self.results

    def __call__(self, nBams, tBams):
        #Fetch,
        logging.info("asm task groupid=%s start" % (self.data.name))
        reads = []
        chrom = self.data.chrom
        start = self.data.start - self.buffer
        start = max(0, start)
        #hope that fetching beyonde 3' boundary is okay
        end = self.data.end + self.buffer
        
        self.leftReads = []
        self.rightReads = []
        self.pbReads = []
        
        for bam in nBams:
            self.__fetchPEReads(bam, chrom, start, end)
        
        for bam in tBams:
            self.pbReads.extend(super(SpadesAssembler, self).fetchReads(bam, chrom, start, end, trim=True))
            #self.pbReads.extend(Assembler.fetchReads(self, bam, chrom, start, end, trim=True))
        
        #Assemble
        totReads = len(self.leftReads) + len(self.rightReads) + len(self.pbReads)
        if totReads > self.args.maxreads:
            return "Failure - Too Many Reads (%d) %s" % (totReads, self.data.name)
        logging.info("assembling %d reads" % (totReads))
        self.result = self.__assemble()
        logging.info("asm task groupid=%s finish" % (self.data.name))
        return self.result
        
        #I can eventually save memory by just returning the output file's name - 
        #Though to do that I 'd need to edit the file with the reads' groupname
        
class CeleraAssembler(Assembler):
    
    def __init__(self, data, buffer, tmpDir, timeout, *args, **kwargs):
        Assembler.__init__(self, data, buffer, tmpDir, timeout)
            
    def __assemble(self, reads):
        """
        writes temp files
        assembles
        reads results
        clears temp files
        returns results as a string
        Calls the assembler
        """
        self.myTmpFiles = []
        #want to make a directory for temp files?
        
        """
        Celera procedure (current version does NOT use specs.. might be a bad idea)
        
        runCA -d directory -p prefix -s specfile <option=value> ... <input-files> ...

        input-files need to be output into FRG
        """
        #Temporary Files
        fout = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, dir=self.tmpDir)
        self.myTmpFiles.append(fout.name)
        
        for name, seq, qual in reads:
            fout.write(">%s\n%s\n" % (name, seq))
        
        fout.close()
        resultOut = tempfile.NamedTemporaryFile(prefix="minia_", delete=False, dir=self.tmpDir)
        
        estSize = self.buffer * 2
        if self.data.rest[0] != 'DEL':
            estSize += int(self.data.rest[1])
        
        r, o, e = exe("minia {reads} 20 3 {estSize} {outPrefix}".format(reads=fout.name, \
                      estSize=estSize, outPrefix=resultOut.name ), \
                      timeout=self.timeout)
        
        logging.debug("RET - %d\nOUT - %s\nERR- %s" % (r, o, e))
        
        self.myTmpFiles.extend([resultOut.name + ".contigs.fa",  resultOut.name + ".debloom", \
                           resultOut.name + ".debloom2", resultOut.name + ".false_positive_kmers", \
                           resultOut.name + ".reads_binary",      resultOut.name + ".solid_kmers_binary"])
        
        if r == 214:
            super(MiniaAssembler, self).cleanupTmp()
            return "Failure - Assembly Timeout " + self.data.name
         
        fasta = FastaFile(resultOut.name + ".contigs.fa")
        
        results = {}
        for key in fasta:
            results[key] = FastqEntry(key, fasta[key], '?' * len(fasta[key]))
        
        #save to file
        fout = tempfile.NamedTemporaryFile(prefix = "asm" + self.data.name, \
                                    suffix=".fastq", delete=False, dir=self.tmpDir)
        for key in results:
            fout.write("@group" + self.data.name + "_" + key + "\n" + \
                        results[key].seq + '\n+\n' + \
                        results[key].qual + '\n')
            
        fout.close()
        self.results = fout.name
        
        #clean up
        super(CeleraAssembler, self).cleanupTmp()
        
        return self.results

    def __call__(self, nBams, tBams):
        #Fetch,
        logging.info("asm task groupid=%s start" % (self.data.name))
        reads = []
        chrom = self.data.chrom
        start = max(0, self.data.start - self.buffer)
        #hope that fetching beyond 3' boundary is okay
        end = self.data.end + self.buffer
        
        for bam in nBams:
            #I actually need to convert these into a FRG file..
            #make temp
            for name, seq, qual, in super(PhrapAssembler, self).fetchReads(bam, chrom, start, end):
                #write read to illumina file
                pass
        
        for bam in tBams:
            #I actually need to convert these into a FRG..
            #And to I need to trim anymore?
            for name, seq, qual, in super(PhrapAssembler, self).fetchReads(bam, chrom, start, end, trim=True):
                #write to pacbio file
                pass
        
        if totReads > self.maxreads:
            return "Failure - Too Many Reads (%d) %s" % (totReads, self.data.name)
        #fastqToCA -insertsize {muins} {stdins} -librayname ill -technology illumina -nonrandom
        #fastqToCA -librayname pac -technology pacbio-raw -nonrandom
        #Assemble
        logging.info("assembling %d reads" % (len(reads)))
        self.result = self.__assemble(illFn, pacFn)
        
        logging.info("asm task groupid=%s finish" % (self.data.name))
        return self.result
        
        #I can eventually save memory by just returning the output file's name - 
        #Though to do that I 'd need to edit the file with the reads' groupname
 
def insertDist(bam):
    """
    Samples reads to get mean insert size and standard deviation
    """
    num_samp = 1000000
    counter = 0
    skip = 5000000
    skip_counter = 0
    ins_list = []
    for read in bam.fetch():
        if skip_counter < skip:
            skip_counter += 1
            continue
        if read.is_proper_pair and not read.is_reverse and not read.is_secondary:
            ins_list.append(read.tlen)
            counter += 1
        if counter == num_samp:
            break
    mean = sum(ins_list)/float(len(ins_list))
    v = 0
    for i in ins_list:
        v += (i-mean)**2
    variance = v/float(len(ins_list))
    stdev = variance**(0.5)
    return (mean, stdev)

def toQual(input):
    if type(input) == str:
        return " ".join([str(ord(x)-33) for x in input])
    elif type(input) == list and type(input[0]) == int:
        return "".join([chr(x+33) for x in input])
    raise TypeError("Expected string or list of ints for toQual")
    
def parseArgs(argv):
    parser = argparse.ArgumentParser(description=USAGE, \
                formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("putative", metavar="BED", type=str, \
                        help="Bed of regions to assemble")
    parser.add_argument("-b", "--bam", type=str, nargs="*", \
                        help="Input Bam (NonTrim)")
    parser.add_argument("-p", "--pacBam", type=str, nargs="*", \
                        help="PacBio Bam")
    parser.add_argument("-a", "--assembler", type=str, default='phrap', choices=["phrap", "minia", "spades"],
                        help="Assembly program to use (%(default)s)")
    parser.add_argument("-B", "--buffer", type=int, default=1000, \
                        help="Amount of buffer sequence around the variant to use (%(default)s)")
    parser.add_argument("-n", "--nproc", type=int, default=1, \
                        help="Number of processors to use (%(default)s)")
    parser.add_argument("-o", "--output", default="asm.fastq",\
                        help="Where to write the resultant assemblies (%(default)s)")
    parser.add_argument("-r", "--reference", default=None, \
                        help="Reference to map to (optional if --noRemap)")
    parser.add_argument("--noRemap", action="store_false", \
                        help="Do not remap assembly")
    parser.add_argument("--noSplitMap", action="store_false", \
                        help="Do not map tails from remapped assembly (off if --noRemap)")
    parser.add_argument("--timeout", type=int, default=30, \
                        help="Timeout assembly after N minutes (%(default)s)")
    parser.add_argument("--maxspan", type=int, default=100000, \
                        help="Maximum Span of SV to attempt assembling (%(default)s)")
    parser.add_argument("--maxreads", type=int, default=500, \
                        help="Maximum number of Illumina reads used to attempt assembling (%(default)s)")
    parser.add_argument("--temp", type=str, default=tempfile.gettempdir(),
                            help="Where to save temporary files")
    parser.add_argument("--start", type=int, default=0,
                        help="Index of the first variant to begin assembling. (%(default)s)")
    parser.add_argument("--stride", type=int, default=1,
                        help="Assemble one every N reads (%(default)s)")
    parser.add_argument("--debug", action="store_true",\
                        help="Verbose Logging")

    #parser.add_argument("--insertsize", type=int, default=None, \
                        #help=("Celera - insert size for PE Illumina reads (auto_detect)"))
    #parser.add_argument("--insertstd", type=float, default=None, \
                        #help=("Celera - insert std for PE Illumina reads (auto_detect)"))
    
    args = parser.parse_args(argv)
    setupLogging(args.debug)
    
    # Parameter checks
    if args.bam is None and args.pacBam is None:
        logging.error("Expected at least one BAM argument")
        exit(1)
    
    if not args.output.endswith(".fastq"):
        logging.error("Output needs to end with .fastq")
        exit(1)
    
    if not os.path.exists(args.putative):
        logging.error("Input {inp} does not exist".format(inp=args.putative))
        exit(1)
    
    if args.noRemap and args.reference == None:
        logging.error("Cannot remap without --reference")
        exit(1)
    
    if args.reference and not os.path.exists(args.reference):
        logging.error("Reference {ref} does not exist".format(ref=args.reference))
        exit(1)
    
    if args.bam is None:
        args.bam = []
        if args.insertsize is None and args.bam is not None:
            j = pysam.Samfile(args.bam[0])
            mu,std = insertDist(j)
            j.close()
            args.insertsize = mu
            args.insertstd = std if args.insertstd is None else args.insertstd
    
    if args.pacBam is None:
        args.pacBam = []
       
    return args


def run(argv):
    args = parseArgs(argv)
    # Establish communication queues
    tasks = multiprocessing.JoinableQueue()
    results = multiprocessing.Queue()
    
    # Start consumers
    num_consumers = args.nproc
    consumers = [ Consumer(tasks, results, args)
                    for i in xrange(num_consumers) ]
    for w in consumers:
        w.start()
    
    # Enqueue jobs
    num_jobs = 0

    #open the putative file
    bed = BedFile.fromFile(args.putative)
    
    entryPos = args.start
    
    ### Set the Assembler
    if args.assembler == 'phrap':
        MyAssembler = PhrapAssembler
    elif args.assembler == 'minia':
        MyAssembler = MiniaAssembler
    elif args.assembler == 'spades':
        MyAssembler = SpadesAssembler
    elif self.assembler == 'dipspades':
        raise NotImplementedError("Working on it")
    else:
        logging.error("%s - Not a valid assembler. Seek --help" % (assembler))
        exit(1)
    
    eout = open(args.output + ".err", 'w')
    while True:
        if entryPos >= len(bed):
            break
        
        entry = bed[entryPos]
        
        if abs(entry.start - entry.end) > args.maxspan:
            eout.write("Too Big %s %d\n" % (entry.name, abs(entry.start-entry.end)))
        else:
            tasks.put(MyAssembler(entry, args))
            num_jobs += 1
        
        #stride over
        #for i in xrange(args.stride):
            #entryPos += 1
        entryPos += args.stride
    
    
    # Add a poison pill for each consumer
    for i in xrange(num_consumers):
        tasks.put(None)
    
    logging.info("%d Tasks" % (num_jobs))
    # Wait for all of the tasks to finish -- I might not need to
    #logging.info("Joining")
    #tasks.join()
        
    fout = open(args.output, 'w')
    #tracking who failed
    asmFails = 0
    # Consolidate results
    logging.info("Consolidating")
    while num_jobs:
        result = results.get()
        num_jobs -= 1
        
        if result.startswith("Failure "):
            logging.error(result)
            eout.write(result + '\n')
            asmFails += 1
            continue
        
        fh = open(result,'r')
        fout.write(fh.read())
        fh.close()
        
        os.remove(result)
    
    eout.close()
    logging.error("%d assemblies failed" % (asmFails))
    fout.close()
    #for each, I'm going to need to 
    #Grabbing my reads
    if args.noRemap:
        logging.info("PIE mapping") 
        r, o, e = exe(("Honey.py pie --temp {tempDir} --nproc {np} "
                       "--minTail 50 {reads} {ref}" \
                       .format(tempDir=args.temp, reads=args.output, \
                               ref=args.reference, np=args.nproc)))
        if r != 0:
            logging.error("Honey pie quit (%d)\n%s\n%s" % (r, o, str(e)))
            exit(r)
        logging.info("Honey Log:\n" + o)

    logging.info("Finished")



if __name__ == '__main__':
    args = parseArgs(sys.argv[1:])
    takeMassivePhrap(args)
