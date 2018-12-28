#!/usr/bin/env python
"""
document the protocol
"""

import os, sys, subprocess, logging, time
from optparse import OptionParser
from string import Template
from xml.etree import ElementTree
from glob import glob

STAGES = ["setup", "mapping", "support", "extraction", "assembly", "output"]


from pbsuite.jelly import Stages
from pbsuite.utils.setupLogging import *
from pbsuite.utils.CommandRunner import * 


USAGE = """USAGE: Jelly.py <stage> <protocol.xml> [-x \"--options for stage\"]
    
    Jelly is the driver for each stage of the 
    reference genome upgrade. 
    
    <stage> is one of
        %s
    <protocol.xml> contains the information about the
    data and parameters Jelly will run. See README.txt
    or the documentation for the Protocol's format.""" % ("\n\t".join(STAGES))

class JellyProtocol():
    """
    Independent method to parse protocols
    """
    def __init__(self, fileName):
        self.protocolName = os.path.abspath(fileName)
        self.parseProtocol()
    
    def parseProtocol(self):
        """
        Parses a Jelly Protocol
        """
        try:
            file = ElementTree.parse(self.protocolName)
        except Exception:
            logging.error(("If you're actually sure the input Protocol is "\
                           "where you said it is, then the Protocol doesn't "\
                           "have valid XML Format. Use an XML validator."))
            sys.exit(1)
        
        root = file.getroot()
        
        refNode = root.find("reference")
        
        if refNode is None:
            logging.error("Protocol doesn't have <reference> element.")
            sys.exit(1)
        else:
            self.scaffoldName = refNode.text
            if not os.path.exists(self.scaffoldName):
                logging.error("Reference %s Does Not Exist!" % (self.scaffoldName))
                sys.exit(1)
            self.referenceNameBase = self.scaffoldName[:self.scaffoldName.rindex(".fasta")]
            self.scaffoldQualName = self.referenceNameBase+".qual"
            if not os.path.exists(self.scaffoldQualName):
                sys.stderr.write(("[WARNING] Couldn't find %s... assuming there " \
                                "isn't a quality file for the reference\n") % \
                                (self.scaffoldQualName))
                self.scaffoldQualName = None
                
            self.gapTable = self.referenceNameBase+".gapInfo.bed"
            
            self.reference = self.scaffoldName
             
            self.referenceIndex = self.reference+".sa"
            
        inputNode = root.find("input")
        if inputNode.attrib.has_key("baseDir"):
            self.baseDir = inputNode.attrib["baseDir"]
        else:
            self.baseDir = ""
        
        self.inputs = []
        for input in inputNode:
            
            #Do some checking on the inputs
            fullPath = os.path.join(self.baseDir, input.text)
            if not os.path.exists(fullPath):
                logging.error("Input %s doesn't exit! Exiting" % (fullPath))
                exit(0)
            if not (fullPath.lower().endswith('.fasta') or fullPath.lower().endswith('.fastq')):
                logging.error("Input %s needs to end with .fasta or .fastq! Exiting" % (fullPath))
                exit(0)
            
            self.inputs.append(fullPath)
        if len(self.inputs) == 0:
            logging.error("Protocol doesn't specifiy any input inputs!")
            sys.exit(1)

        
        outputNode = root.find("outputDir")
        if outputNode is None:
            logging.warning("Output directory not specified. Using pwd. Hope you're cool with that...")
            self.outDir = os.getcwd()
        else:
            self.outDir = outputNode.text
        
        #Command Runner is going to sort all this out.
        #And create an object to do the calls.
        self.runXML = root.find("cluster")
        
        blasrNode = root.find("blasr")
        if blasrNode is None:
            logging.warning("No blasr parameters!?")
            self.blasrParams = ""
        else:
            self.blasrParams = blasrNode.text
    

class JellyRunner():
    """
    Take a JellyProtocol and loads in the variables in a way that 
    JellyRunner can use it easily!
    """
    def __init__(self):
        """
        Given a protocol fn, load it up so we are ready to run. 
        """
        self.parseArgs()
        setupLogging(self.options.debug)
        sys.stderr.write("""
Please Cite: English, Adam C., Stephen Richards, Yi Han, Min Wang,
             Vanesa Vee, Jiaxin Qu, Xiang Qin, et al. "Mind the
             Gap: Upgrading Genomes with Pacific Biosciences RS
             Long-Read Sequencing Technology." PLoS ONE 7, no. 11
             (November 21, 2012): e47768.
             doi:10.1371/journal.pone.0047768.\n\n""")
        self.parseProtocol()
        
    def parseProtocol(self):
        self.protocol = JellyProtocol(self.protocolName)
        self.parseCluster(self.protocol.runXML)
        
    def parseCluster(self, xmlNode):
        if xmlNode is None:
            self.runCmd = CommandRunner()
        else:
            command = xmlNode.find("command")
            if command is None:
                logging.error(("You're trying to use a cluster " \
                               "template, but you didn't specify the " \
                               "template. Please read the documentation." \
                               "Exiting.\n"""))
                sys.exit(1)
            nJobs = xmlNode.find("nJobs")
            
            if nJobs is None or nJobs.text == '0':
                logging.warning(("Not Specifying number of jobs may " \
                                  "overload clusters."))
                nJobs = 0
            else:
                nJobs = int(nJobs.text)
            
            cmdTemplate = command.text
            self.runCmd = CommandRunner(cmdTemplate, nJobs)
                    
    def parseArgs(self):
        """
        Uses OptionParser to parse out input
        Jelly.py <stage> <protocol>
        """
        parser = OptionParser(USAGE)
        parser.remove_option("-h")
        parser.add_option("-h", "--help", action="store_true", default=False)
        
        parser.add_option("--debug",action="store_true",default=False)
        parser.add_option("-x", dest="extras", type="string", default="", 
                help="-x \"<options>\" are options to pass into the stage you're running")
        
        self.options, args = parser.parse_args()

        if self.options.help == True:
            if len(args) == 1:
                if args[0] in STAGES:
                    print exe(Stages.PRINT_HELPS[args[0]])[1]
                    sys.exit(0)
                #Else, this will drop down to the next parser.error
            else:
                print parser.format_help()
                sys.exit(0)
        if len(args) != 2 or args[0] not in STAGES:
            parser.error("Invalid Arguments. Expected one of\n'%s'" % "', '".join(STAGES))
            sys.exit(1)
        self.executeStage = args[0]
        self.protocolName = os.path.abspath(args[1])
        
        
    def run(self):
        logging.info("Executing Stage: %s" % self.executeStage)
        
        if self.options.debug:
            Stages.DEBUG = "--debug"
        
        #Setup before a stage and run
        if self.executeStage == "setup":
            wDir = os.path.dirname(self.protocol.scaffoldName)
            myCommands = [Stages.setup(self.protocol.scaffoldName, \
                        self.protocol.scaffoldQualName, self.protocol.gapTable, \
                        self.options.extras)]
        
        elif self.executeStage == "mapping":
            wDir = os.path.join(self.protocol.outDir, "mapping")
            try:
                os.mkdir(wDir)
            except OSError:
                logging.debug("%s already exists" % wDir)
            if not os.path.exists(wDir):
                logging.warning("%s was not created. Check write permissions!" % wDir)
                sys.exit(1)
            
            myCommands = Stages.mapping(self.protocol.inputs, wDir, \
                        self.protocol.reference, self.protocol.referenceIndex, \
                        self.protocol.blasrParams, self.options.extras)
        
        elif self.executeStage == "support":
            wDir = os.path.join(self.protocol.outDir, "support")
            try:
                os.mkdir(wDir)
            except OSError:
                logging.debug("%s already exists" % wDir)
            if not os.path.exists(wDir):
                logging.warning("%s was not created. Check write permissions!" % wDir)
                sys.exit(1)
       
            myCommands = Stages.support(self.protocol.outDir, self.protocol.gapTable, \
                    wDir, self.options.extras) 
        
        elif self.executeStage == "extraction":
            wDir = os.path.join(self.protocol.outDir, "assembly")
            try:
                os.mkdir(wDir)
            except OSError:
                logging.debug("%s already exists" % wDir)
            if not os.path.exists(wDir):
                logging.warning("%s was not created. Check write permissions!" % wDir)
                sys.exit(1)
            
            myCommands = [Stages.extraction(wDir, \
                                            self.protocolName, \
                                            self.options.extras)]
        
        elif self.executeStage == "assembly":
            wDir = os.path.join(self.protocol.outDir, "assembly")
            
            myCommands = Stages.assembly(wDir, self.protocol.gapTable, \
                                         self.options.extras)
        
        elif self.executeStage == "output":
            wDir = os.path.join(self.protocol.outDir, "assembly")
            
            myCommands = [Stages.collection(self.protocol.outDir, self.protocol, self.options.extras)]
        
        logging.debug("CommandRunner Returned: " + 
            str(self.runCmd(myCommands, wDir, self.executeStage )) )
        
        logging.info("Finished %s Stage: %s" % (self.runCmd.runType, self.executeStage))
        
if __name__ == '__main__':
    prog = JellyRunner()
    prog.run()
