#!/usr/bin/env python
import sys, re, os, logging, subprocess
from optparse import OptionParser

from pbsuite.utils.setupLogging import *
from pbsuite.utils.FileHandlers import FastaFile, QualFile, Gap, wrap, qwrap
from pbsuite.utils.CommandRunner import exe

USAGE = """%prog <inputScaffolding.fasta> [outputContigs.fasta]

Take the input scaffolding and split it into contigs around all [Nn]+ 
If an <inputScaffolding.qual> exists beside the <inputScaffolding.fasta>,
we'll create and incorporate it into PBJelly. Otherwise, it is not necessary
"""

refParser = re.compile("(.*)\|(ref\d{7})/?(\d+)?$")

class Setup():
    def __init__(self):
        self.parseArgs()
        setupLogging(self.opts.debug)
    
    def parseArgs(self):
        parser = OptionParser(USAGE)
        
        parser.add_option("-g", "--gapOutput", dest="gapOutput", \
            help="Create the table for gapInformation", default=None)
        parser.add_option("-i", "--index", dest="index", action="store_true", default=False, \
            help="Create the .sa index for faster blasr alignments")
        #parser.add_option("--noRename", dest="rename", action="store_true", default=False, \
            #help="Flag to indicate there will be no renaming of scaffolding.")
        parser.add_option("--debug",action="store_true", help="Increases verbosity of logging" )
        parser.add_option("--minGap", dest="minGap", type="int", default=25, \
            help="Minimum number of consecutive Ns to be considered a gap. default=25")
        parser.add_option("--maxGap", dest="maxGap", type="string", default="", \
            help="Maximum number of consecutive Ns to be considered a gap default=Inf")
        
        self.opts, args = parser.parse_args()
        if len(args) < 1:
            parser.error("Error! Incorrect number of arguments")
        if len(args) == 1:
            self.scaffInput = args[0]
            if not self.scaffInput.endswith(".fasta"):
                parser.error("Reference must end in extension .fasta! Please rename it.")
            if not os.path.isfile(self.scaffInput):
                parser.error("Error! Scaffold File is not a file / does not exist")
        else:
            parser.error("Error! Incorrect number of arguments")

        qn = self.scaffInput[:self.scaffInput.rindex('.fasta')]
        qualInputName = qn+".qual"
        if not os.path.isfile(qualInputName):
            self.qualInput = None
        else:
            self.qualInput = qualInputName
   
    def run(self):
        #Fasta Ref Output
        scaffTempName = self.scaffInput+".tempFasta"
        scaffOutput = open(scaffTempName, 'w')
        
        #Qual Ref Output
        if self.qualInput is not None:
            qualTempName= self.qualInput+".tempQual"
            qualOutput = open(qualTempName, 'w')
        
        #Gaps Output
        if self.opts.gapOutput is not None:
            gapTableOut = open(self.opts.gapOutput,'w')
        else:
            gapTableOut = False
        
        logging.info("Creating reference sequence index names and identifying gaps")
        
        refTemplate = "ref%07d"
        refId = 1
        
        #Read References
        reference = FastaFile(self.scaffInput)
        if self.qualInput is not None:
            qualReference = QualFile(self.qualInput)    
        
        for key in reference:
            
            scaffIndex = refTemplate % refId
            scaffName = key.replace(' ','_')
            
            refId += 1
            
            scaffName = scaffName + "|" + scaffIndex
            scaffOutput.write(">"+scaffName+"\n"+wrap(reference[key])+"\n")
            
            if self.qualInput is not None:
                qualOutput.write(">"+scaffName+"\n"+qwrap(qualReference[key])+"\n")
            
            gapCoords = []
            for gap in re.finditer("[^Nn]([Nn]{%d,%s})[^Nn]" % \
                    (self.opts.minGap, self.opts.maxGap), reference[key]):
                gapCoords.append([gap.start() + 1, gap.end() - 1])
            
            if len(gapCoords) == 0:#no Gaps
                gapTableOut.write("\t".join([scaffName, 'na', 'na', scaffIndex+"_0_0", '3'])+'\n')
                logging.debug("Scaffold %s is empty" % scaffName)
                continue
            
            #Consolidate gaps that are too close -- indicating LQ regions.
            i = 0
            while i < len(gapCoords)-1:
                if gapCoords[i+1][0] - gapCoords[i][1] < 25:
                    gapCoords[i+1][0] = gapCoords[i][0]
                    del(gapCoords[i])
                else:
                    i += 1
            
            prevEnd = 0#Contig Start Tracking
            idx = 0
            #Make the first gap
            prevEnd = gapCoords[0][1]
            gapCoords[0][1]-gapCoords[0][0]
            
            flag = Gap.BEGIN
            if len(gapCoords) == 1:
                flag += Gap.END
            if gapTableOut:
                gapTableOut.write("%s\t%i\t%i\t%s_%i_%i\t%d\n" \
                        % (scaffName, gapCoords[0][0], gapCoords[0][1], scaffIndex, idx, idx+1, flag))

            #Now Go Through the rest of the gaps
            for i in range(1, len(gapCoords)):
                idx += 1
                prevEnd = gapCoords[i][1]
                gapCoords[i][1]-gapCoords[i][0]
            
                if gapTableOut:
                    if i == len(gapCoords)-1:
                        flag = Gap.END
                    else:
                        flag = 0
                    gapTableOut.write("%s\t%i\t%i\t%s_%i_%i\t%d\n" \
                        % (scaffName, gapCoords[i][0], gapCoords[i][1], scaffIndex, idx, idx+1, flag))
            
        #Close shop
        scaffOutput.close()
        os.rename(self.scaffInput, self.scaffInput+".original")
        os.rename(scaffTempName, self.scaffInput)
        
        if self.qualInput is not None:
            qualOutput.close()
            os.rename(self.qualInput, self.qualInput+".original")
            os.rename(qualTempName, self.qualInput)
        
        if gapTableOut:
            gapTableOut.close()
        
        if self.opts.index:
            logging.info("Creating .sa indexes for references")
            r, o, e = exe("sawriter %s.sa %s" % (self.scaffInput, self.scaffInput))
            if r != 0:
                logging.error("sawriter returned %d" % r)
                logging.error("Ensure it's in your path")
                exit(1)
            logging.debug(str(o) + ' ' + str(e))
        
        logging.info("Finished!")

if __name__ == '__main__':
    me = Setup()
    me.run()
