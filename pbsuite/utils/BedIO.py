
"""
Bed Entry and Files objects
"""

class BedEntry:
    """
    Holds a Bed Entry with accessors to chrom, start, end, name
    all other fields are held in BedEntry.rest as a list
    """
    def __init__(self, chrom, start, end, name, *args):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.rest = list(args)
    
    def plainStr(self):
        """
        Doesn't write all the infomation, just the first 4 columns
        """
        return "%s\t%d\t%d\t%s" % (self.chrom, self.start, self.end, self.name)
    
    def __str__(self):
        extra = ""
        if len(self.rest) > 0:
            extra = "\t" + "\t".join([str(x) for x in self.rest])
        
        return self.plainStr() + extra
    
    def __lt__(self, other):
        if self.chrom != other.chrom:
            return cmp(self.chrom, other.chrom)
        return self.start < other.start
    
    def __gt__(self, other):
        if self.chrom != other.chrom:
            return cmp(self.chrom, other.chrom)
        return self.start > other.start
 
class BedFile(list):
    """
    Create a .bed file either from a file or from a set of entries
    """
    def __init__(self, fileName=None, entries=None):
        super(list)
        if entries is not None:
            self.extend(entries)
        self.fileName = fileName

    @classmethod
    def fromFile(cls, fileName):
        """
        Create a BedFile directly from file
        """
        fh = open(fileName, 'r')
        entries = []
        for line in fh.readlines():
            entries.append(BedEntry(*line.strip().split('\t')))
        fh.close()
        return cls(fileName, entries)
        
    def plainStr(self):
        """
        Doesn't write all the infomation, just the first 4 columns
        """
        ret = ""
        for i in self:
            ret += i.plainStr() + '\n'
        return ret
        
    def __str__(self):
        ret = ""
        for i in self:
            ret += str(i)+"\n"
        return ret
 
class BedPEEntry():
    """
    Same as a bed except the less than doesn't sort for chromosomes
    """
    def __init__(self, chrom1, start1, end1, chrom2, start2, end2, name, *args):
        self.chrom1 = chrom1
        self.start1 = start1
        self.end1 = end1
        self.chrom2 = chrom2
        self.start2 = start2
        self.end2 = end2
        self.name = name
        self.rest = args
    
    def plainStr(self):
        """
        Doesn't write all the infomation, just the first 7 columns
        """
        return "%s\t%d\t%d\t%s\t%d\t%d\t%s" % (self.chrom1, self.start1, self.end1, \
                                               self.chrom2, self.start2, self.end2, self.name)
    
    def __str__(self):
        extra = ""
        if len(self.rest) > 0:
            extra = "\t" + "\t".join([str(x) for x in self.rest])
        
        return self.plainStr() + extra
    
    def __lt__(self, other):
        #if self.chrom != other.chrom:
            #return cmp(self.chrom, other.chrom)
        return self.start < other.start
    
    def __gt__(self, other):
        #if self.chrom != other.chrom:
            #return cmp(self.chrom, other.chrom)
        return self.start > other.start
 
