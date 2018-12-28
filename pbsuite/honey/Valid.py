import sys, os, shutil, logging, argparse
from collections import namedtuple

from pbsuite.utils.CommandRunner import exe
from pbsuite.utils.FileHandlers import FastaFile, mergeFastaQual, FastqFile

USAGE = """\
    This procedure will validate a list of sites in a .vcf by
    A) samtools viewing the region of interest into an individual bam
    B) Creating a .fastq from that bam
    C) Assembling the .fastq result until the haplotype is satisfied or fails
    D) Compare the results to the call made in the .vcf  <-- not yet implemented
    E) Create an IGB url for quick visualization         <-- not yet implemented

None of this works anymore :-D
"""

reference = "/stornext/snfs0/next-gen/pacbio/data/references/human_g1k_v37/sequence/human_g1k_v37.fasta"

EXITCODES = {"OSError" : 1,\
             "ExeError": 2 }

USAGEPLUS = """
Notes

Main Commands:
    samtools view aligned_reads.bam ${region} | fix output

    awk '{print "@"  $1  "\n"  $10  "\n+\n"  $11}' > input.fastq

    OLCAssembly.py input.fastq -nproc 4 

    mkdir level2 && cd level2
    OLCAssembly.py ../out.fasta ../out.qual -nproc 4

    blasr out.fasta ~/pacbio/data/references/human_g1k_v37/sequence/human_g1k_v37.fasta -sa ~/pacbio/data/references/human_g1k_v37/sequence/human_g1k_v37.fasta.sa -nproc 4 -out remap.sam -sam -bestn 1
"""

VCFEntry = namedtuple("VCFEntry", "region sample haplotype")

def exeLog( func ):
    """
    Decorator for logging exe'd commands
    noFail=True will not allow the failure of this exe to end program's execution
    """
    def inner(*args, **kwargs):
        logging.info("Executing %s %s %s" % (func.__name__, args, kwargs))
        r,o,e = func(*args, **kwargs)
        if r != 0:
            logging.error("Problem executing %s" % (func.__name__))
            logging.error("\tRETCOD %d" % (r))
            logging.error('\tSTDOUT' + o.strip())
            logging.error('\tSTDERR' + str(e))
            if not noFail:
                exit(r)
        else:
            logging.info("%s exit success %s" % (func.__name__, r))
    inner.__doc__ = func.__doc__
    return inner

def exeLog_noFail( func ):
    """
    Decorator for logging exe'd commands
    noFail=True will not allow the failure of this exe to end program's execution
    """
    def inner(*args, **kwargs):
        logging.info("Executing %s %s %s" % (func.__name__, args, kwargs))
        r,o,e = func(*args, **kwargs)
        if r != 0:
            logging.error("Problem executing %s" % (func.__name__))
            logging.error("\tRETCOD %d" % (r))
            logging.error('\tSTDOUT' + o.strip())
            logging.error('\tSTDERR' + str(e))
            logging.error("Continuting")
        else:
            logging.info("%s exit success %s" % (func.__name__, r))
    return inner

def iterVCF(fn):
    """
    yield each vcfEntry
    """
    fh = open(fn,'r')
    for line in fh.readlines():
        if line.startswith("##"):
            continue
        if line.startswith("#"):
            h = line.strip()[1:].split('\t')
            header = {}
            for pos,i in enumerate(h):
                header[i] = pos
            continue
        
        data = line.strip().split()
        pos = int(data[header["POS"]])
        region = "{0}:{1}-{1}".format(data[header["CHROM"]], pos)
        
        #info grab
        info = data[header["INFO"]].split(';')
        sample = None
        haplotype = None
        for i in info:
            if i.startswith("validation_sample"):
                sample = i.split('=')[1]
                haplotype = data[header[sample]]
                break
        
        if sample is not None:
            yield VCFEntry(region, sample, haplotype)

def fixPBSam(fn):
    """
    The reference in the PB pipeline has fasta entry names that are too
    verbose. This will change the SQ names in the .bam so that they're
    usable for downstream processes
    """
    fh = open(fn,'r')
    output = []
    for line in fh.readlines():
        if line.startswith("@"):
            if line.startswith("@SQ"):
                #split fix
                data = line.split("\t")
                data[1] = data[1].split(' ')[0]
                line = "\t".join(data)
                output.append(line)
            else:
                output.append(line)
            continue
        data = line.split("\t")
        data[2] = data[2].split(' ')[0]
        line = "\t".join(data)
        output.append(line)
    fh.close()
    fout = open(fn,'w')
    fout.write("".join(output))
    fout.close()

def setupSite( entry ):
    """
    Makes the directory we'll be working in
    """
    if not os.path.exists(entry.sample):
        try:
            os.mkdir(entry.sample)
        except OSError:
            logging.error("Can't make directory %s" % entry.sample)
            exit(EXITCODES["OSError"])
            
    outPath = os.path.join(entry.sample, entry.region.replace(':','_'))
    if not os.path.exists(outPath):
        try:
            os.mkdir(outPath)
        except OSError:
            logging.error("Can't make directory %s" % outPath)
            exit(EXITCODES["OSError"])
    return outPath

def outputFastq( fastq, fn ):
    """
    writes the fastqs to fn
    """
    fout = open(fn, 'w')
    for i in fastq.values():
        fout.write(i.toString())
    fout.close()
    
@exeLog
def grabReads( inputBam, entry, outFn ):
    """
    Gets all of the reads for a region and puts them into outFn
    """
    return exe("samtools view -h %s %s > %s" % (inputBam, entry.region, outFn))

@exeLog
def bam2sam( fn, outName):
    """
    Turns a bam to a sam
    """
    return exe("samtools view -h %s > %s " % (fn, outName))
    
@exeLog
def sam2bam( fn ):
    """
    Creates BAM from SAM (only setup for hg19 -- see global variable reference)
    """
    name = fn[:-4]
    return exe(("samtools view -bt {0} {1} | samtools sort - {2}.sort && "
                "mv {2}.sort.bam {2}.bam && "
                "samtools index {2}.bam").format(reference, fn, name))

@exeLog
def samToFastq( inSam, outFq ):
    """
    Creates input.fastq from SAM file
    """
    return exe(('grep -v "^@" %s | '
                'awk \'{print "@"  $1  "\\n"  $10 "\\n+\\n"  $11}\' '
                '> %s') % (inSam, outFq))
@exeLog
def remapReads( reads, outName):
    """
    remaps reads to the provided reference (only setup for hg19 -- see 
    global variable reference)
    """
    return exe("blasr {0} {1} -sa {1}.sa -nproc 4 -out {2} -sam -bestn 1"\
               .format(reads, reference, outName))
    
@exeLog
def pileup( bam ):
    """
    create a pileup from the bam
    """
    return exe("samtools mpileup -f {0} {1} > {1}.plup".format(reference, bam))
    
def iterAssemble(entry, myReads):
    """
    """
    @exeLog_noFail
    def assemble(inputFq, workDir):
        return exe("OLCAssembly.py %s --nproc 4 --fqOut --workDir %s" % (inputFq, workDir))
    
    level = 0
    curReads = myReads
    
    ref,alt = entry.haplotype.split('/')
    isHet = ref != alt
    
    while True:
        logging.info("Running assembly level %d" % (level))
        workDir = os.path.join(os.path.dirname(myReads), "level%d" % (level))
        
        #Potential problem
        try:
            os.mkdir(workDir)
        except OSError:
            pass
            
        assemble(curReads, workDir)
        outName = os.path.join(workDir,"out.fastq")
        if not os.path.exists(outName):
            logging.error("Assembly iteration didn't return consensus")
            logging.error("Manual checking is required(?)")
            logging.error("Returning the best answer we have")
            return curReads
            
        output = FastqFile(outName)
        if len(output) == 0:
            logging.error("Couldn't assemble contigs after %d levels" % (level))
            logging.error("Returning the best answer we have")
            return curReads
        
        if not isHet and len(output) == 1 or isHet and len(output) == 2:
            logging.info("Made consensus after %d levels" % (level))
            return outName
        elif isHet and len(output) == 1:
            logging.warning("One consensus sequence created for het at level %d" % (level))
            logging.warning("Manual checking is required(?)")
            logging.error("Returning the best answer we have")
            return outName
        
        level += 1
        curReads = outName
    
def setupLogging(debug=False):
    logLevel = logging.DEBUG if debug else logging.INFO
    logFormat = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig( stream=sys.stderr, level=logLevel, format=logFormat )
    logging.info("Running %s" % " ".join(sys.argv))

def parseArgs():
    parser = argparse.ArgumentParser(description=USAGE, \
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bam", metavar="bam", type=str, \
                        help="BAM containing mapped reads")
    parser.add_argument("vcf", metavar="vcf", type=str, \
                        help="VCF containing sites to be validated")
    parser.add_argument("-s", "--formatsam", dest="formatsam", action="store_true",\
                        help="Reformat the input BAM's sequence names in-place")
    parser.add_argument("-m", "--maxdepth", dest="maxdepth", default=sys.maxint, \
                        help="Downsample to a maximum depth per site")
    parser.add_argument("-d", "--debug", dest="debug", action="store_true",\
                        help="Print verbose logging")
    args = parser.parse_args()
    
    return args
    
if __name__ == '__main__':
    args = parseArgs()
    setupLogging(args.debug)
    inputBam = args.bam
    inputVCF = args.vcf
    
    if args.formatsam:
        inputSam = inputBam[:-4] + ".sam"
        bam2sam(inputBam, inputSam)
        fixPBSam(inputSam)
        sam2bam(inputSam)
    
    for entry in iterVCF(sys.argv[2]):
        logging.info("Creating folders for site %s" % str(entry))
        outDir = setupSite(entry)
        logging.info("Extracting region of interest")
        myBam = os.path.join(outDir, "region.sam")
        grabReads(inputBam, entry, myBam)
        myReads = os.path.join(outDir, "region.fastq")
        samToFastq(myBam, myReads)
        sam2bam(myBam)
        
        logging.info("Assembling consensus for %s haplotype" % (entry.haplotype))
        fastqOut = iterAssemble(entry, myReads)
        logging.info("Remapping consensus")
        consensus = os.path.join(outDir, "consensus.fastq")
        shutil.copyfile(fastqOut, os.path.join(outDir, "consensus.fastq"))
        blasrOut = os.path.join(outDir, "remap.sam")
        remapReads(fastqOut, blasrOut)
        sam2bam(blasrOut)
        logging.info("Creating Pileup")
        pileup(blasrOut[:-4]+'.bam')
        logging.info("Finished %s" % str(entry))
    
    
