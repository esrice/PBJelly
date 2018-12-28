#!/usr/bin/env python

import argparse, sys
from pbsuite.utils.setupLogging import *
from pbsuite.honey import bampie, TGraf, HSpots, Force, ComplexResolver, \
                         massivePhrap

STAGES = {"pie":   bampie.run, \
          "tails": TGraf.run, \
          "spots": HSpots.run, \
          "force": Force.run , \
          "asm":  massivePhrap.run, \
          "cpxres": ComplexResolver.run}

USAGE = """\
   Honey - genomic variant calling with long sequencing reads

   STAGE is one of
     pie        Extract and map soft-clipped Tails from a bam.
     tails      Cluster mapped tails to make break-points of larger events.
     spots      Find genomic variants within reads' spans.
     force      Given a BedFile of predicted variants, force search for matching
     asm        Assemble reads around a variant and remap contigs to a reference
     cpxres     Complex multi-break-point resolution (beta)
    
   See HoneyReadme.txt for documentation or --help for details\
"""

def parseArgs():
    parser = argparse.ArgumentParser(prog="Honey.py", description=USAGE, \
            formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument("-h", "--help", action="store_true")
    parser.add_argument("stage", metavar="STAGE", choices=STAGES.keys(), type=str, \
                        help="Stage to execute")
    parser.add_argument("options", metavar="OPTIONS", nargs=argparse.REMAINDER,\
                         help="Options to pass to the stage")

    args = parser.parse_args()
    
    sys.stderr.write("""
Please Cite: English, Adam C., William J. Salerno, Jeffery G.
             Reid. "PBHoney: identyfying genomic variants via
             long-read discordance and interrupted mapping."
             BMC Bioinformatics 2014, 15:180 (June 10, 2014).
             doi:10.1186/1471-2105-15-180\n\n""")
    
    STAGES[args.stage](args.options)
    
if __name__ == '__main__':
    parseArgs()
