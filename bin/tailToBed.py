#!/usr/bin/env python
import sys

if __name__ == '__main__':
    fh = open(sys.argv[1],'r')
    fh.readline()#args
    h = fh.readline()[1:].strip()
    header = {}
    for pos, item in enumerate(h.split('\t')):
        exec("%s=%d" % (item, pos))
    
    for line in fh.readlines():
        data = line.strip().split()
        if data[annot] == "TLOC":
            #I can't do tlocs, yet
            continue
        if data[annot] in ["INS","DEL"] and \
            abs(int(data[uBreak]) - int(data[dBreak])) < int(data[remainSeq]):
                size = int(data[remainSeq])
        elif data[annot] in ["INS", "DEL"]:
            size = abs(int(data[uBreak]) - int(data[dBreak]))
        elif data[annot] == "INV":
            #does it have all of the pieces to say it's fully an inversion
            nbp = len(data[evidence].split(';'))
            if nbp == 8:
                data[annot] = "INV"
            elif nbp <= 4:
                data[annot] = "MIS"
            else:
                data[annot] = "BKP"
            size = abs(int(data[uBreak]) - int(data[dBreak]))
        else:
            print "PROBLEM", line,
        
        print "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(data[uRef], \
                        data[uBreak], data[dBreak], \
                        data[id], data[annot], size)
    fh.close()




