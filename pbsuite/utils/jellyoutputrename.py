from pbsuite.utils.FileHandlers import FastaFile
import json, re, sys

"""
## Arguments

1 - reference.fasta -- input reference created by Setup.py
2 - liftOverTable.json - created at end of Jelly run
3 - jelly.out.fasta -- new reference created by Jelly

"""
fasta = FastaFile(sys.argv[1])
nameLookup = {}
for entry in fasta:
    data = entry.split('|')
    refKey = data[-1]
    origName = "|".join(data[:-1])
    nameLookup[refKey] = origName

liftOver = json.load(open(sys.argv[2],'r'))
jellyFasta = FastaFile(sys.argv[3])
regex = re.compile("ref\d{7}")
for key in liftOver:
    myRefs = set()
    for id, strand, size in liftOver[key]:
        myRefs.update(regex.findall(id))
    newName = []
    for refId in myRefs:
        newName.append(nameLookup[refId])
    sys.stdout.write(">%s\n%s\n" % ("_".join(newName), jellyFasta[key]))

    

