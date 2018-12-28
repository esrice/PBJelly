#!/usr/bin/env python
import sys, os, argparse
from pbsuite.honey.HSpots import SpotResult
from pbsuite.utils import VCFIO

USAGE = "Turn .spots results into a .vcf"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=USAGE, \
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", metavar="SPOTS", \
                        help="Results to convert")
    
    args = parser.parse_args()
    
    template = VCFIO.VCFFile(os.path.join(os.path.dirname(VCFIO.__file__),"vcfTemplate.vcf"))
    vcfout = VCFIO.VCFFile(args.input+'.vcf', template, 'w')
    
    fh = open(args.input,'r')
       
    for line in fh.readlines():
        if line.startswith("##"):
            vcfout.HEAD["parameters"] = line.strip()[2:]
            continue
        elif line.startswith("#"):
            continue
        
        spot = SpotResult.parseLine(line)
        info = {"SVTYPE":spot.svtype, "END":str(spot.end)}

        filter = []
        if "noSpan" in spot.tags:
            filter.append("NoSpan")
            format = []
            sample = []#might be able to put some of hes in.. coverage? alt?
        else:
            format = ['GT', 'GQ', 'COV', 'ALT', 'ALNID', 'MQF'] 
            sample = [spot.tags["GT"], spot.tags["GQ"], int(spot.tags["coverage"]), \
                      spot.tags["strandCnt"], float(spot.tags["alnIdentityEstimate"]), \
                      int(spot.tags["mqfilt"])]
            if sample[0] == '0/0':
                filter.append("LowQual")
            else:
                filter.append("PASS")
        
        if spot.svtype == 'DEL':
            #this is fucked
            info["SVLEN"] = '-'+str(spot.size)
            ref = spot.tags["seq"] if "seq" in spot.tags else '.'
            alt = ['.']
        elif spot.svtype == 'INS':
            info["SVLEN"] = str(spot.size)
            ref = '.'
            alt = [spot.tags["seq"] if "seq" in spot.tags else '.']
        
        entry = VCFIO.VCFEntry(spot.chrom, spot.start, '.', ref, alt, '.', \
                               filter, info, format, {"SAMPLE":sample})
        vcfout.addEntry(entry)
    
       
    for i in vcfout:
        print i
            
    fh.close()

