import sys, random
from pbsuite.utils.FileHandlers import FastaFile, revComp, wrap

def getRandomSeq(length):
    return "".join([random.choice(['A', 'T', 'C', 'G']) for i in xrange(length)])
    
if __name__ == '__main__':
    fasta = FastaFile(sys.argv[1])
    key = fasta.keys()[0]
    ref = list(fasta[key])
    
    #800bp insertion in the sample (deletion in the reference) 
    ref[5000:5800] = ""
    #5000 Insertion

    #Inversion in the sample (inversion in the reference) tails
    ref[9000:12000] = list("".join(ref[10000:13000]).translate(revComp)[::-1])
    #9000-12000 - INversion
    
    #1kb deletion in sample (insert into the reference) tails
    seq = getRandomSeq(1000)
    ref[20000:20000] = list(seq)
    #20000-21000 -- Deletion 
    
    #100bp insertion in sample (deletion in the reference) spots
    ref[30000:30100] = ""
    #30000 - Insertion

    #200bp deletion in sample (insert into the reference) spots
    seq = getRandomSeq(200)
    ref[35000:35000] = list(seq)
    
    #35000 - 35200 -- Deletion
    print ">%s\n%s" % (key, wrap("".join(ref)))    

