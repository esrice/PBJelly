import sys, re, os
from collections import OrderedDict
from pbsuite.utils.setupLogging import *

#setupLogging(debug=True)

HONTEMPLATE = os.path.join(os.path.dirname(__file__),"vcfTemplate.vcf")

class VCFFile():
    """
    Holds all the entries.. no indexing, yet
    """

    def __init__(self, filename, template=None, mode='r'):
        """
        Either [r]ead , [w]rite, or [a]ppend to a vcf file
        """
        #Open or write a file:
        self.filename= filename
        self.mode = mode
        self.filehandler = open(filename, mode)
        
        #INIT the metadata
        if template is not None:
            self.HEAD = template.HEAD
            self.FILTER = template.FILTER
            self.INFO = template.INFO
            self.ALT = template.ALT
            self.FORMAT = template.FORMAT
            self.SAMPLES = template.SAMPLES
        else:
            self.HEAD = OrderedDict()
            self.FILTER = OrderedDict()
            self.INFO = OrderedDict()
            self.ALT = OrderedDict()
            self.FORMAT = OrderedDict()
            self.SAMPLES = OrderedDict()
        
        if mode in ['r', 'a']:
            if template is not None:
                logging.warning("Overwriting any META data template has with file's")
            self.__parseFile()
        else:
            self.__entries = []
    
    #Need an iterator
    def __iter__(self):
        return self.__entries.__iter__()
    
    def __parseFile(self):
        self.__entries = []
        self.filehandler.seek(0)
        for line in self.filehandler.readlines():
            if line.startswith("##"):
                data = META.fromString(line)
                if data.id is None:
                    self.HEAD[data.key] = data.description
                elif data.key == "FILTER":
                    self.FILTER[data.id] = data
                elif data.key == "INFO":
                    self.INFO[data.id] = data
                elif data.key == "ALT":
                    self.ALT[data.id] = data
                elif data.key == "FORMAT":
                    self.FORMAT[data.id] = data
                else:
                    logging.error("Unrecognized META key %s" % data.key)
            elif line.startswith("#"):#Header
                for pos, name in enumerate(line.strip().split('\t')[9:]):
                    self.SAMPLES[name] = pos
            else:
                self.__parseEntry(line)
    
    def __parseEntry(self, line):
        try:    
            entry = VCFEntry.fromString(line)
            self.addEntry(entry)
        except ValueError as e:
            logging.error("Couldn't parse line %s" % line.strip())
            logging.error(e)
            exit(1)
    
    def addEntry(self, vcfentry):
        ##Validation and meta appending...
        #meta appending and validation
        vcfentry.fileValidate(self)
        self.__entries.append(vcfentry)

    def write(self, outName):
        """
        Output to new file
        """
        fout = open(outName, 'w')
        for i in self.HEAD:
            fout.write(str(self.HEAD[i])+'\n')
        for i in self.FILTER:
            fout.write(str(self.FILTER[i])+'\n')
        for i in self.INFO:
            fout.write(str(self.INFO[i])+'\n')
        for i in self.ALT:
            fout.write(str(self.ALT[i])+'\n')
        for i in self.FORMAT:
            fout.write(str(self.FORMAT[i])+'\n')
        fout.write(VCFEntry.header + "\t" + "\t".join(self.SAMPLES.keys()) + '\n')
        for entry in self.__entries:
            fout.write(str(entry) + '\n')
        


METARE = re.compile(('##\s*(?P<key>\w+)\s*=\s*\<\s*ID\s*=\s*(?P<id>\w+)\s*'
                     ',\s*(Number\s*=\s*(?P<number>[0-9.])\s*,'
                     '\s*Type\s*=\s*(?P<type>\w+)\s*,)?'
                     '\s*Description\s*=\s*"(?P<description>.*)"\s*\>\s*'),\
                     re.IGNORECASE)

class META():
    """
    Key=<ID=x,Number=y,Type=z,Description="everything has one of these">
     or
    Key=some text description
    Holds the rules on how to parse the various META types
    If id == None, we assume some simple text passed in through
    description
    """
    def __init__(self, key, id=None, number=None, type=None, description=None):
        self.key = key
        self.id = id
        self.number = number
        self.type = type
        self.description = description
        
    @classmethod
    def fromString(cls, line):
        match = METARE.match(line)
        if match is None:
            try:
                key,val = line.lstrip("#").strip().split('=')
                key = key.strip()                
                val = val.strip()
                return cls(key=key, description=val)
            except Exception as e:
                logging.error("Error parsing %s\n" % line)
                logging.error(str(e)+'\n')
                exit(1)
        
        g = match.groupdict()
        if g['number'] == '':
            g['number'] = None
        return cls(**g)
   
    def convert(self, data):
        """
        Based on the meta tag's rules, convert this data
        """
        if self.type is None:
            return data
        
        if self.number == '1':
            items = [data]
        else:
            items = data.split(',')

        if self.type.lower() == 'integer':
            convert = int
        elif self.type.lower() == 'float':
            convert = float
        elif self.type.lower() == 'flag':
            convert = bool
        elif self.type.lower() == 'string':
            convert = str
        else:
            convert = lambda x: x
        
        newItems = []
        for i in items:
            newItems.append(convert(i))
        
        return newItems
    
        
    def __str__(self):
        if self.id is None:
            return "##%s=%s" % (self.id, self.description)
        
        return '##{0}=<ID={1},Number={2},Type={3},Description="{4}">'\
               .format(self.key, self.id, self.number, self.type, self.description)

class VCFEntry():
    """
    CHROM   the chromosome.
    POS the genome coordinate of the first base in the variant. Within a chromosome, VCF
        records are sorted in order of increasing position.
    ID  a semicolon-separated list of marker identifiers.
    REF the reference allele expressed as a sequence of one or more A/C/G/T nucleotides
        (e.g. "A" or "AAC")
    ALT the alternate allele expressed as a sequence of one or more A/C/G/T nucleotides
        (e.g. "A" or "AAC"). If there is more than one alternate alleles, the field should 
        be a comma-separated list of alternate alleles.
    QUAL    probability that the ALT allele is incorrectly specified, expressed on the the
            phred scale (-10log10(probability)).
    FILTER  Either "PASS" or a semicolon-separated list of failed quality control filters.
    INFO    additional information (no white space, tabs, or semi-colons permitted).
    FORMAT  colon-separated list of data subfields reported for each sample. The format
            fields in the Example are explained below.
    
    Only CHROM and POS are required. Everything else will be blank if necessicary
    """
    header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    
    def __init__(self, CHROM, POS, ID=None, REF=None, ALT=None, QUAL=None, \
                 FILTER=None, INFO=None, FORMAT=None, SAMPLES=None):
        """
        #CHROM  str
        #POS    int
        #ID     str or None
        #REF    str or None
        #ALT    list or None
        #QUAL   float or None
        #FILTER list or None
        #INFO   dict or None
        #FORMAT list or None
        #SAMPLE dict or None
        """
        
        #Populate VCF Entry
        self.CHROM = str(CHROM)
        self.POS = int(POS)
        self.ID = str(ID) if ID is not None else '.'
        self.REF = str(REF) if REF is not None else '.'
        
        if ALT is not None:
            if type(ALT) == list:
                self.ALT = [re.sub('[<>]', '', x) for x in ALT if x != '']
            else:
                raise ValueError("Expected ALT to be list")
                
        else:
            self.ALT = []
        self.QUAL = float(QUAL) if QUAL not in [None, '.'] else '.'
        if FILTER is not None:
            if type(FILTER) == list:
                self.FILTER = FILTER
            else:
                raise ValueError("Invalid Value for FILTER")
                
        else:
            self.FILTER = []
        
        if INFO is not None:
            if type(INFO) == dict:
                self.INFO = INFO
            else:
                raise ValueError("Invalid Value for INFO")
        else:
            self.INFO = {}

        if FORMAT is not None:
            if type(FORMAT) == list:
                self.FORMAT = FORMAT
            else:
                raise ValueError("Invalid Value for INFO")
        else:
            self.FORMAT = []

        if SAMPLES is not None:
            self.SAMPLES = SAMPLES
            if type(SAMPLES) == dict:
                self.SAMPLES = SAMPLES
            else:
                raise ValueError("Invalid Value for SAMPLE")
        else:
            self.SAMPLES = {}
        
    @classmethod
    def fromString(cls, line):
        data = line.strip().split('\t')
        chrom = data[0]
        pos = int(data[1])
        id = data[2]
        ref = data[3]
        alt = data[4].split(',') if data[4] != '.' else []
        qual = float(data[5]) if data[5] != '.' else None
        filter = data[6].split(';') if data[6] != '.' else []
        info = {}
        for i in data[7].split(';'):
            k,v =  i.split('=')
            info[k] = v
        format = data[8].split(':') if data[8] != '.' else []
        samples = {}
        for idx, samp in enumerate(data[9:]):
            samples[idx] = samp.split(':')
        
        return cls(chrom, pos, id, ref, alt, qual, filter, info, format, samples)
    
    def fileValidate(self, vcffile):
        """
        validates entry as if it's from vcffile .. so all info must exist
        in header
        updates any blanks (e.g. indexed samples) with metadata info 
        makes a dict out of the samples information 
        This is super useful for when you're reading in a file
        If you're creating new entries, use VCFFile.addEntry()
        and it'll ensure that everything you're adding from that entry is
        already documented within the VCFFile's metadata
        """
        self.__myparent = vcffile
        
        for i in self.ALT:
            if i.upper() not in vcffile.ALT and len(re.sub(i.upper(), "[ATCGN.]", '')) != 0:
                raise ValueError("Unknown ALT %s in entry %s" % (i, str(self)))

        for i in self.FILTER:
            if i not in vcffile.FILTER and i != "PASS":
                raise ValueError("Unknown FILTER %s in entry %s" % (str(i), str(self)))
        
        for i in self.INFO:
            if i not in vcffile.INFO:
                raise ValueError("Unknown INFO %s in entry %s" % (i, str(self)))
            self.INFO[i] = vcffile.INFO[i].convert(self.INFO[i])
        
        for i in self.FORMAT:
            if i not in vcffile.FORMAT:
                raise ValueError("Unknown FORMAT %s in entry %s" % (i, str(self)))
        
        newSamp = {}       
        for key in self.SAMPLES:    
            if type(key) == int:
                try:
                    newKey = vcffile.SAMPLES.keys()[key]
                    newSamp[newKey] = self.SAMPLES[key]
                except IndexError:
                    raise ValueError("More SAMPLES than VCFFile has in entry %s" % (str(self)))
            elif key not in vcffile.SAMPLES:
                raise ValueError("Unknown SAMPLE %s" % (str(key)))
            else:
                newSamp[key] = self.SAMPLES[key]
        self.SAMPLES = newSamp
    
    def __str__(self):
        if self.__myparent is None:
            logging.warning("Writing entry without parent. Sample order may suffer")
            nALT = ",".join(self.ALT) if len(self.ALT) != 0 else "."
            nINFO = ";".join(["%s=%s" % (x,",".join([str(y) for y in self.INFO[x]])) for x in self.INFO])
            nFORMAT = ':'.join(self.FORMAT)
            nSAMPLES = []
            for s in self.SAMPLES:
                dat = []
                for i in nFORMAT:
                    dat.append(i)
                    
        else:
            nALT = []
            for i in self.ALT:
                if i in self.__myparent.ALT:
                    nALT.append("<%s>" % i)
                else:
                    nALT.append(i)
            nALT = ",".join(nALT)
            
            nINFO = []
            for i in self.__myparent.INFO:
                if i in self.INFO:
                    if self.__myparent.INFO[i].number == '1':
                        nINFO.append("%s=%s" % (i, self.INFO[i][0]))
                    else:
                        nINFO.append("%s=%s" % (i, ",".join([str(y) for y in self.INFO[i]])))
            nINFO = ";".join(nINFO)
            
            nFORMATd = []
            for i in self.__myparent.FORMAT:
                if i in self.FORMAT:
                    nFORMATd.append(i)
            
            nSAMPLES = []
            for i in self.__myparent.SAMPLES:
                if i in self.SAMPLES:
                    nSAMPLES.append(":".join([str(x) for x in self.SAMPLES[i]]))
                else:
                    nSAMPLES.append(':'.join(['.']*len(nFORMATd)))
            
            if len(nFORMATd) == 0:
                nFORMAT = '.'; nSAMPLES = ['.']
            else:
                nFORMAT = ':'.join(nFORMATd)
        
        return "\t".join([self.CHROM, \
                          str(self.POS),\
                          self.ID, \
                          self.REF, \
                          nALT, \
                          str(self.QUAL), \
                          ";".join(self.FILTER), \
                          nINFO, \
                          nFORMAT, \
                          "\t".join(nSAMPLES)])
                            
        #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
        #otherwise we'll format correctly
        #for i in self.__myparent:
        #Need . for unknown samples
        #

def test(filename):
    f = VCFFile(filename, 'r')
    for i in f:
        print str(i)
    return f

if __name__ == '__main__':
    test(sys.argv[1])
    test(sys.argv[1])
