import sys
from pbsuite.utils.FileHandlers import FastqFile, M5File
from pbsuite.jelly.Support import AlignmentConnector, SUPPORTFLAGS

"""
Need to do work here
"""
if __name__ == '__main__':
    connector = AlignmentConnector()
    aligns = connector.parseAlignments(M5File(sys.argv[1]))

    reads = FastqFile(sys.argv[2])

    bestScore = None
    best = None
    fout = open("reads.fastq",'w')
    spanCount = 0
    for readGroup in aligns:
        if readGroup[0].qname.startswith("ref"):
            continue
        if len(readGroup) == 2:
            r1, r2 = readGroup
            a = connector.extendsTarget(r1)
            b = connector.extendsTarget(r2)
            if a != SUPPORTFLAGS.none and b != SUPPORTFLAGS.none:
                spanCount += 1
                print r1.qname, "spans"
                
                rStart = min(r1.qend, r2.qend)
                rEnd = max(r1.qstart, r2.qstart)
                t = reads[r1.qname].subSeq(rStart, rEnd)
                fout.write(str(t))
                gout = open("seed%d.fasta" % spanCount, 'w')
                gout.write(">%s\n%s\n" % (t.name, t.seq))
                gout.close()
                if bestScore is None:
                    bestScore = r1.score + r2.score
                    best = reads[r1.qname].subSeq(rStart, rEnd)
                else:
                    if (r1.score + r2.score) < bestScore:
                        best = reads[r1.qname].subSeq(rStart, rEnd)
        else:
            a = readGroup[0]
            if a.tname.endswith('e5'):
                fout.write(str(reads[a.qname].subSeq(0, a.qstart)))
            elif a.tname.endswith('e3'):
                fout.write(str(reads[a.qname].subSeq(a.qend, a.qseqlength)))
    fout.close()
    print "%d spans" % spanCount
    fout = open("seed.fasta",'w')
    fout.write(">%s\n%s\n" % (best.name, best.seq))
    fout.close()
                
            
