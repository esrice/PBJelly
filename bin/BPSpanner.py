#!/usr/bin/env python
import sys, math, json
import pysam
from pbsuite.honey.Force import *

fh = open(sys.argv[1])
fh.readline(); fh.readline() #Header

bam = pysam.Samfile(sys.argv[2])

#400bp around bp must be spanned
def bpCheck(bam, chrom, point, BUFFER=200):
    nSpan = 0
    nCount = 0
    spans = {}
    tot = {}
    for read in bam.fetch(reference=chrom, start = max(0, point-BUFFER), end=point+BUFFER):
        #Get the shawties outta chere
        if read.flag & 0x40 or read.flag & 0x80:
            continue
        if read.pos >= point-BUFFER and read.aend <= point:
            #print "ever?"
            continue
        if read.pos >= point and read.aend <= point-BUFFER:
            #print "ever?"
            continue
        
        if read.qname in tot:
            continue
        
        if read.pos <= point-BUFFER and point+BUFFER <= read.aend:
            spans[read.qname] = 1
        tot[read.qname] = 1
        #elif read.pos < point and point < read.aend:
        #nCount += 1
    return len(spans), len(tot)

def log_choose(n, k):
    r = 0.0
    # swap for efficiency if k is more than half of n
    if k * 2 > n:
        k = n - k

    for  d in xrange(1,k+1):
        r += math.log(n, 10)
        r -= math.log(d, 10)
        n -= 1

    return r

# return the genotype and log10 p-value
def bayes_gt(ref, alt):
    # prob seeing an alt read with true genotype of of hom_ref, het, hom_alt respectively
    p_alt = [0.2, 0.4, 0.8]

    total = ref + alt
    
    lp_homref = log_choose(total, alt) + alt * math.log(p_alt[0], 10) + ref * math.log(1 - p_alt[0], 10)
    lp_het = log_choose(total, alt) + alt * math.log(p_alt[1], 10) + ref * math.log(1 - p_alt[1], 10)
    lp_homalt = log_choose(total, alt) + alt * math.log(p_alt[2], 10) + ref * math.log(1 - p_alt[2], 10)

    return (lp_homref, lp_het, lp_homalt)

gtlookup = {0:'0/0', 1:'0/1', 2:'1/1'}
print "id\tnReads\tnZMWs\tavgCov\tpctVar\tSpan1\tCoverage1\tPct1\tSpan2\tCoverage2\tPct2\tgt\tgt1\tgt2"
for line in fh.readlines():
    data = line.strip().split("\t")
    chr1 = data[2]
    bp1 = int(data[3])
    chr2 = data[5]
    bp2 = int(data[6])

    nReads = int(data[10])
    nZMWs = int(data[11])
    
    nSpan1, nCount1 = bpCheck(bam, chr1, bp1)
    nSpan2, nCount2 = bpCheck(bam, chr2, bp2)
    if nCount1 == 0 or nCount2 == 0:
        continue
    nRef = nSpan1 + nSpan2
    nAlt = (nCount1 + nCount2) - nRef
    
    p_bp1 = bayes_gt(nSpan1, nCount1 - nSpan1)
    p_bp2 = bayes_gt(nSpan2, nCount2 - nSpan2)
    p_bp = bayes_gt(nRef, nAlt)
    
    
    gt1 = gtlookup[p_bp1.index(max(p_bp1))]
    gt2 = gtlookup[p_bp2.index(max(p_bp2))]
    gt = gtlookup[p_bp.index(max(p_bp))]
    
    avgCov = (nCount1 + nCount2) / 2

    var = (float(nReads)/avgCov)
    p1 = (float(nSpan1)/nCount1)
    p2 = (float(nSpan2)/nCount2)

    threshold = 0.15
    #if the percent spans are really low, it's hom variant
    if p1 <= threshold and p2 <= threshold:
        #print 'HOM', data[keygenoType], data[keyisvalid], line,
        genoType = 'HOM'
    #if the percent spans is really close to 50%
    elif 0.5-threshold <= p1 <= 0.5+threshold and 0.5-threshold <= p2 <= 0.5+threshold:
        genoType = "HET"
    #if the fraction of variant reads is close to threshold
    elif abs(var - 0.5) < threshold:
        genoType = "HET*"
    #if the fraction of variant reads is closer to 0
    elif abs(var - 1) < threshold:
        genoType = "HOM*"
    else:#Complex genotype.. I don't know what's happening
        #print 'COX', data[keygenoType], data[keyisvalid], line,
        genoType = "COX"
    
    output = {}
    print "\t".join([str(x) for x in [data[0], \
                                      nReads, \
                                      nZMWs, \
                                      avgCov, \
                                      "%.3f" % var, \
                                      nSpan1, \
                                      nCount1, \
                                      "%.3f" % p1, \
                                      nSpan2, \
                                      nCount2, \
                                      "%.3f" % p2, \
                                      genoType, \
                                      gt, gt1, gt2]])
        
    sys.stdout.flush()
