#!/usr/bin/env python

from pbsuite.utils.FileHandlers import FastqFile, M5File
from pbsuite.utils.CommandRunner import exe

"""
This can be run inside of an assembly folder and create our filling sequence into
polish.out.fasta
"""
if __name__ == '__main__':
    input = FastqFile("input.fastq")

    fout = open("ref.fasta",'w')
    for i in input.values():
        if i.name.startswith("ref"):
            fout.write(">%s\n%s\n" % (i.name, i.seq))
    fout.close()
    
    print exe(("blasr input.fastq ref.fasta  -bestn 2 -m 5 -noSplitSubreads > out.m5"))
    print exe(("python /stornext/snfs5/next-gen/scratch/english/Jelly/"
               "DevJelly/branches/consensusDev/GetSubs.py out.m5 input.fastq"))
    print exe(("python /stornext/snfs5/next-gen/scratch/english/Jelly/"
               "DevJelly/branches/sv/pbjPolish.py "
               "reads.fastq seed.fasta -n 4 -l"))

