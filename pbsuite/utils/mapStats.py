#!/usr/bin/env python

"""
Parse a directory of m4 files and report mapping stats like:

Number of Reads
Number of Subreads
Mean Subread Length
Mean Alignment Length
Mean #Hits per Read
Mean Accuracy
Mean Score

"""

import sys, os, glob
from collections import defaultdict
from pbsuite.utils.FileHandlers import M4File
from pbsuite.utils.summarizeAssembly import *

def collectStats(stats, input):
    if input.endswith(".m4"):
        parser = M4File(input)
    elif input.endswith(".m5"):
        parser = M5File(input)
    
    reads = defaultdict(list)
    subreads = defaultdict(list)    
    
    nSubs = 0
    numBases = 0
    numAlignedBases = 0
    tailLengths = []
    
    for line in parser:
        name = line.qname.split('/')[1]
        reads[name].append(line)
        subreads[line.qname].append(line)
    
    reads = dict(reads)
    subreads = dict(subreads)
    
    for read in reads:
        readLength = {}
        alignedLength = 0
        for subr in reads[read]:
            readLength[subr.qname] = subr.qseqlength
            alignedLength += subr.qend - subr.qstart
        stats["read"].append(sum(readLength.values()))
        stats["aln_read"].append(alignedLength)
    
    for sub in subreads.keys():
        readLength = subreads[sub][0].qseqlength
        stats["subread"].append(readLength)
        
        alignedBases = 0
        pieHits = []
        for hit in subreads[sub]:
            alnLength = hit.qend - hit.qstart
            alignedBases += alnLength
            pieHits.append(alnLength)
            
            stats["hit"].append(alnLength)
            stats["hitsim"].append(hit.pctsimilarity)
            
        stats["aln_subread"].append(alignedBases)
        
        if len(subreads[sub]) > 1:
            stats["tail"].append(readLength)
            stats["aln_tail"].extend(pieHits)
        
        stats["unmappedTail"] += readLength - alignedBases

if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.stderr.write("Error! Invalid Number of Arguments\nMapStats.py <mappingDir | input.m4>\n Generates stats for every m4 alignemnt file in an inputDir or for a single m4 file\n")
        sys.exit(1)
    
    inputDir = sys.argv[1]
    stats = {"read": [],
             "aln_read": [],
             "subread": [],
             "aln_subread":[],
             "hit": [],
             "hitsim": [],
             "tail": [],
             "aln_tail": [],
             "unmappedTail": 0}
    
    if os.path.isfile(inputDir):
        if not inputDir.endswith('.m4'):
            sys.stderr.write("Error! Specified file is not .m4!")
            exit(1)
        collectStats(stats,inputDir)
    else:
        a = glob.glob(os.path.join(inputDir, "*.m4"))
        a.extend(glob.glob(os.path.join(inputDir, "*.m5")))
        for input in a:
            collectStats(stats, input)
    
    print "Read Stats:"
    read = getStats(stats["read"])
    alnread = getStats(stats["aln_read"])
    print "Number of Reads\t", read["numSeqs"]
    print "Number of Bases\t", read["totalLength"]
    print "Mean Read Length\t", read["mean"]
    print "N50 Read Length\t", read["n50"]
    print "Number of Aligned Bases\t", alnread["totalLength"]
    print "Mean Aligned Read Length\t", alnread["mean"]
    print "N50 Aligned Read Length\t", alnread["n50"]
    print "Unmapped Bases\t", stats["unmappedTail"]
    print 
    print "Subread Stats:"
    subread = getStats(stats["subread"])
    alnsubread = getStats(stats["aln_subread"])
    print "Number of Subreads\t", subread["numSeqs"]
    print "Total Subread Bases\t", subread["totalLength"]
    print "Mean Subread Length\t", subread["mean"]
    print "N50 Subread Length\t", subread["n50"]
    print "Total Subread Aligned Bases\t", alnsubread["totalLength"]
    print "Mean Subread Aligned Length\t", alnsubread["mean"]
    print "N50 Subread Aligned Length\t", alnsubread["n50"]
    print 
    print "Alignment Hit Stats:"
    hit = getStats(stats["hit"])
    hitsim = getStats(stats["hitsim"])
    print "Number of Alignments\t", hit["numSeqs"]
    print "Total Hit Aligned Bases\t", hit["totalLength"]
    print "Mean Hit Aligned Length\t", hit["mean"]
    print "N50 Hit Aligned Length\t", hit["n50"]
    print "Mean Hit Percent Similarity\t", hitsim["mean"]
    print 
    print "Tailed (Split) Hits Stats:"
    tailhit = getStats(stats["tail"])
    alntail = getStats(stats["aln_tail"])
    print "Number of Tailed Reads\t", tailhit["numSeqs"]
    print "Number of Tailed Bases\t", tailhit["totalLength"]
    print "Mean Tailed Read Length\t", tailhit["mean"]
    print "N50 Tailed Read Length\t", tailhit["n50"]
    print "Total Tailed Aligned Bases\t", alntail["totalLength"]
    print "Mean Tailed Alignment Length\t", alntail["mean"]
    print "N50 Tailed Alignment Length\t", alntail["n50"]
    print 
    
