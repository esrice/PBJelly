#!/usr/bin/env python

import sys, logging, argparse
import pysam
from pbsuite.honey.HSpots import *
from pbsuite.utils.setupLogging import setupLogging
USAGE = "Recall spots in a hon.h5 file"


args = parseArgs(sys.argv[1:], established=True)

setupLogging(args.debug)

f = h5py.File(args.hon,'a')
bam = pysam.Samfile(args.bam)
tsp = 0
makeKernals(args.binsize)
print "#CHROM\tOUTERSTART\tSTART\tINNERSTART\tINNEREND\tEND\tOUTEREND\tTYPE\tSIZE\tINFO"
for chrom in f.keys():
    logging.info("Calling %s" % (chrom))
    container = f[chrom]["data"]
    start = f[chrom].attrs["start"]
    spots = callHotSpots(container, start, args)
    logging.info("Filtering spots")
    fspot = 0
    for spot in spots:
        spot.chrom = chrom
        spot.offset(start)
                
        spot.estimateSize()
        #the sv spans too far
        if spot.size > args.spanMax:
            continue
                #the sv doesn't have adequate read support
        if supportingReadsFilter(bam, spot, args):
            continue
        
        spot.estimateSize()#get a better size estimate
        #Filter based on minimum size expectaions
        if spot.size < args.minIndelSize or spot.size > args.spanMax:
            continue
        
        fspot += 1
        print spot
    tsp += fspot
    logging.info("Found %d spots" % (fspot))
logging.info("Finished %d spots" % (tsp))
