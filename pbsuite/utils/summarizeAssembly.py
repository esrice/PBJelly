#!/usr/bin/env python
import sys, re, math
from optparse import OptionParser
from string import Template
from FileHandlers import FastaFile

USAGE="""%prog <file.fasta> [options]
Returns basic statistics (like N50s) about an assembly"""

def parseArgs():
    parser = OptionParser(USAGE)
    parser.add_option("-b", "--binsize", dest="binsize", type="int", default=0,\
                        help="Bin size for creating gap frequency data. (Default is to not print the frequency)")
    parser.add_option("-m", "--min", dest="min", type="int", default=25,\
                        help="Minimum gap size to be considered. DEFAULT=25")
    parser.add_option("-M", "--max", dest="max", type="str", default="",\
                        help="Maximum gap size to be considered. DEFAULT=inf")
    parser.add_option("-c", "--consolidate",dest="consolidate", type=int, default=25,\
                        help=("Concolidate gaps within XXbp as a single gap since this is more" \
                              "indicative of a single LowQuality Region than multiple gaps. DEFAULT=25 [0 means off]"))
    #Put this in. I think it would be nice.
    #parser.add_option("-o","--output", dest="output", type="str", default=None,
                        #help="File name to output a Bed File with chromosome, start, end of gaps features. Default=Off")
    opts, args = parser.parse_args()
    
    if len(args) != 1:
        parser.error("No Fasta Specified!")
    if opts.min < 1:
        parser.error("Minimum gap size must be at least 1.")
    return opts, args[0]
   
def getStats(seqLengths):
    data = {}

    seqLengths.sort(reverse=True)
    
    data["numSeqs"] = len(seqLengths)
    data["totalLength"] = sum(seqLengths)
    tl = data["totalLength"]
    n50_mark = data["totalLength"] * .5
    n90_mark = data["totalLength"] * .90
    n95_mark = data["totalLength"] * .95
    
    data["n50"] = None
    data["n90"] = None
    data["n95"] = None
    basesSeen = 0
    
    for n in seqLengths:
        basesSeen += n
        if data["n50"] is None and basesSeen > n50_mark:
            data["n50"] = n
        if data["n90"] is None and basesSeen > n90_mark:
            data["n90"] = n
        if data["n95"] is None and basesSeen > n95_mark:
            data["n95"] = n
            break
    #may not have gaps
    if data["numSeqs"] == 0:
        return data
    data["min"] = seqLengths[-1]
    data["FstQu"] = seqLengths[ int(math.floor(data["numSeqs"]*.75)) ]
    median = data["numSeqs"]*.50
    data["median"] = int( (seqLengths[ int(math.floor(median)) ] + \
                           seqLengths[ int(math.floor(median)) ]) / 2)
    data["mean"] = data["totalLength"]/data["numSeqs"]
    data["TrdQu"] = seqLengths[ int(math.floor(data["numSeqs"]*.25)) ] 
    data["max"] = seqLengths[0]

    return data
            
def printBins(seq, binsize):
    """
    Print Bin Sizes.
    """
    if binsize < 1:
        exit(0)

    seq.sort()
    
    BINSIZE = binsize
    bin_mark = BINSIZE
    bincount = 0
    i = 0
    while i < len(seq):
        if seq[i] <= bin_mark:
            bincount += 1
            i += 1
        else:
            if bincount != 0:
                print str(bin_mark-BINSIZE+1)+"bp : "+str(bin_mark)+"bp\t"+str(bincount)
            bincount = 0
            bin_mark += BINSIZE

    if bincount != 0:
        print str(bin_mark-BINSIZE+1)+"bp : "+str(bin_mark)+"bp\t"+str(bincount)

if __name__ == '__main__':
    opts, ref = parseArgs()
    
    reference = FastaFile(ref)
    
    gapLengths = []
    lowQualNs = 0
    contigLengths = []
    contigLengthsNoN = []
    scaffoldLengths = []
    scaffoldLengthsNoN = []
    gapRE = re.compile("[^Nn]([Nn]{%d,%s})[^Nn]" % (opts.min, opts.max))
    for entry in reference:
        seq = reference[entry]
        mySeqLen = len(seq)
        myGapLen = []
        gapCoords = []
        
        for gap in gapRE.finditer( seq ):
            #Finditer gives the full span of the match.
            #The first and last characters of the match are not N's
            #Therefore they are not part of the gap
            gapCoords.append([gap.start() + 1, gap.end() - 1])
        
        if len(gapCoords) == 0:
            contigLengths.append(len(seq))
            ns = seq.count('N') + seq.count('n')
            lowQualNs += ns
            contigLengthsNoN.append(len(seq) - ns)
            scaffoldLengths.append( mySeqLen )
            scaffoldLengthsNoN.append( mySeqLen )
            continue
        
        #Consolidate gaps that are too close
        i = 0
        while i < len(gapCoords)-1:
            if gapCoords[i+1][0] - gapCoords[i][1] < opts.consolidate:
                gapCoords[i+1][0] = gapCoords[i][0]
                del(gapCoords[i])
            else:
                i += 1
        
        contigLengths.append(gapCoords[0][0])
        ns = seq[:gapCoords[0][0]].count('N') + \
            seq[:gapCoords[0][0]].count('n')
        lowQualNs += ns
        contigLengthsNoN.append(gapCoords[0][0] - ns)
        myGapLen.append(gapCoords[0][1]-gapCoords[0][0])

        for i in range(1, len(gapCoords)):
            size = gapCoords[i][0] - gapCoords[i-1][1]
            contigLengths.append(size)
            ns = seq[gapCoords[i-1][1]:gapCoords[i][0]].count('N') \
                 + seq[gapCoords[i-1][1]:gapCoords[i][0]].count('n')
            lowQualNs += ns
            contigLengthsNoN.append(size - ns)
            myGapLen.append(gapCoords[i][1] - gapCoords[i][0])
        size = len(seq) - gapCoords[-1][1]
        contigLengths.append(size)
        ns = seq[gapCoords[-1][1]:].count('N') \
                 + seq[gapCoords[i-1][1]:].count('n')
        lowQualNs += ns

        contigLengthsNoN.append(size - ns)
            
        gapLengths.extend(myGapLen)
        scaffoldLengths.append( mySeqLen )
        scaffoldLengthsNoN.append( mySeqLen - sum(myGapLen) )
        
        #prevStart = 0 # previous contig start
        #contigLengths.append(gap.start() - prevStart - 1)
        #prevStart = gap.end() - 1
        #contigLengths.append(len(seq) - prevStart)
    
    
    scafStats = getStats(scaffoldLengths)
    scafStats2 = getStats(scaffoldLengthsNoN)
    contStats = getStats(contigLengths)
    contStats2 = getStats(contigLengthsNoN)
    gapStats = getStats(gapLengths)
    
    space = str(max([len(str(x)) for x in scafStats.values()])+2)
    
    report = ("#Seqs  | {numSeqs:%d,}\n"
              "Min    | {min:%d,}\n"
              "1st Qu.| {FstQu:%d,}\n" + \
              "Median | {median:%d,}\n" + \
              "Mean   | {mean:%d,}\n" + \
              "3rd Qu.| {TrdQu:%d,}\n" + \
              "Max    | {max:%d,}\n" + \
              "Total  | {totalLength:%d,}\n" + \
              "n50    | {n50:%d,}\n" + \
              "n90    | {n90:%d,}\n" + \
              "n95    | {n95:%d,}\n").replace("%d", str(space))
        
    reportDoub = ("#Seqs  | {numSeqs:%d,}\n" 
                  "Min    | {min:%d,} | {noNMin:%d,}\n" 
                  "1st Qu.| {FstQu:%d,} | {noN1q:%d,}\n" 
                  "Median | {median:%d,} | {noNmed:%d,}\n" 
                  "Mean   | {mean:%d,} | {noNmea:%d,}\n" 
                  "3rd Qu.| {TrdQu:%d,} | {noN3q:%d,}\n" 
                  "Max    | {max:%d,} | {noNmax:%d,}\n" 
                  "Total  | {totalLength:%d,} | {noNtot:%d,}\n" 
                  "n50    | {n50:%d,} | {noNn50:%d,}\n" 
                  "n90    | {n90:%d,} | {noNn90:%d,}\n" 
                  "n95    | {n95:%d,} | {noNn95:%d,}\n").replace("%d",str(space))
    
    if scafStats["numSeqs"] == 0:
        print "="*20
        print "No Scaffolding!"
        print "="*20
    else:
        scafStats["noNMin" ] = scafStats2["min"]
        scafStats["noN1q"]   = scafStats2["FstQu"]
        scafStats["noNmed"]  = scafStats2["median"]
        scafStats["noNmea"]  = scafStats2["mean"]
        scafStats["noN3q"]   = scafStats2["TrdQu"]
        scafStats["noNmax"]  = scafStats2["max"]
        scafStats["noNtot"]  = scafStats2["totalLength"]
        scafStats["noNn50"]  = scafStats2["n50"]
        scafStats["noNn90"]  = scafStats2["n90"]
        scafStats["noNn95"]  = scafStats2["n95"]
        print "="*20
        print "Scaffolds | withGaps | withoutGaps"
        print "="*20
        print reportDoub.format(**scafStats)
        print "="*20
    
    if contStats["numSeqs"] == 0:
        print "No Contigs! (or gaps betwen them)"
        print "="*20
    else:
        contStats["noNMin" ] = contStats2["min"]
        contStats["noN1q"]   = contStats2["FstQu"]
        contStats["noNmed"]  = contStats2["median"]
        contStats["noNmea"]  = contStats2["mean"]
        contStats["noN3q"]   = contStats2["TrdQu"]
        contStats["noNmax"]  = contStats2["max"]
        contStats["noNtot"]  = contStats2["totalLength"]
        contStats["noNn50"]  = contStats2["n50"]
        contStats["noNn90"]  = contStats2["n90"]
        contStats["noNn95"]  = contStats2["n95"]
        print "Contigs | withNs | withoutNs"
        print "="*20
        print reportDoub.format(**contStats)
        print "="*20
    
    if gapStats["numSeqs"] == 0:
        print "No Gaps!"
        print "="*20
    else:
        print "Gaps"
        print "="*20
        print report.format(**gapStats)
        print "="*20
    print "Non-gapped Ns Count: ", lowQualNs
    if opts.binsize != 0:
        printBins(gapLengths, opts.binsize)
