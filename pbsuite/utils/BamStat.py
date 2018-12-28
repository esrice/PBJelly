#!/usr/bin/env python
import sys, json, re, logging, argparse

from pysam import Samfile
from pbsuite.utils.summarizeAssembly import getStats

def expandAlign(alignment):
    """
    Takes a pysam Alignment and creates 
    (reference, query) alignments
    For example:
        query     =  ATCGC-GT
        reference =  AT-GCGGA
        Where C inserted, G deleted, A->T Sub
    """
    seq = alignment.seq
    cigar = expandCigar(alignment.cigar)

    mdTag = None
    for i in alignment.tags:
        if i[0] == "MD":
            mdTag = expandMd(i[1])
    
    if mdTag is None:# and alignment.target:
        logging.debug(("Mapped read %s doesn't have MD tag. Mismatches will be 0") \
                         % (alignment.qname))
        mdTag = "-" * len(cigar)
    
    qPos = 0
    tPos = 0
    
    tSeq = []
    qSeq = []
    #Expanding query seq and filling in target seq
    try:
        for i in cigar:
            if i == 0:
                qSeq.append(seq[qPos])
                if mdTag[tPos] != '-':
                    tSeq.append(mdTag[tPos])
                else:
                    tSeq.append(seq[qPos])
                qPos += 1
                tPos += 1
            elif i == 1:
                qSeq.append(seq[qPos])
                tSeq.append('-')
                qPos += 1
            elif i == 2:
                qSeq.append('-')
                tSeq.append(mdTag[tPos])
                tPos += 1
    except IndexError:
        return None, None
    return (qSeq,tSeq)

def expandCigar(cigar):
    """
    Turns the abbreviated cigar into the full array
    
    0 = M
    1 = I
    2 = D
    """
    ret = []
    for t,s in cigar:
        ret.extend([t]*s)
    return ret

def expandMd(md):
    """
    Turns abbreviated MD into a full array
    """
    ret = []
    for i in re.findall("\d+|\^?[ATCGN]+", md):
        if i.startswith('^'):
            ret.extend(list(i[1:]))
        elif i[0] in ["A","T","C","G","N"]:
            ret.extend(list(i))
        else:
            ret.extend(['-']*int(i))
    return ret
    
def counter(query, reference):
    mat = ins = dels = sub = 0
    for q,r in zip(query, reference):
        if q == '-':
            dels += 1
        elif r == '-':
            ins += 1
        elif q != r:
            sub += 1
        else:
            mat += 1

    acc = mat/float(mat+ins+dels+sub)
    tot = len(filter(lambda x: x != '-', query))
    return acc, tot, ins, dels, sub
    
if __name__ == '__main__':
    if len(sys.argv) == 1:
        sam = Samfile("-", "r")
    elif len(sys.argv) == 2:
        try:
            sam = Samfile(sys.argv[1])
        except ValueError:
            sys.stderr.write("%s may not be a valid sam/bam file\n" % (sys.argv[1]))
            exit(0)
    else:
        sys.stderr.write("Expected 1 argument - the sam/bam file\n")
        exit(0)
    
    readLengths = []
    
    accuracy = 0
    insertions = 0
    deletions = 0
    subs = 0
    soft = 0
    tot = 0.0
    cnt = 0.0
    unmapped = 0
    
    for align in sam:
        if align.is_unmapped:
            unmapped += 1
            continue
            
        query, refer = expandAlign(align)
        if query is None:
            continue#does this happen?
        readLengths.append(len(align.seq))
        cnt += 1
        a,t,i,d,s = counter(query, refer)
        sc = 0
        if align.cigar[0][1] == 4:
            sc += align.cigar[0][0]
        if align.cigar[-1][1] == 4:
            sc += align.cigar[0][0]
            
        accuracy += a
        tot += t
        insertions += i
        deletions += d
        subs += s
        soft += sc
    
    errCnt = float(insertions + deletions + subs)
    stats = getStats(readLengths)
    space = str(max([len(str(x)) for x in stats.values()])+2)
    report = ("#Seqs  | {numSeqs:%d,}\n"
              "Min    | {min:%d,}\n"
              "1st Qu.| {FstQu:%d,}\n" + \
              "Median | {median:%d,}\n" + \
              "Mean   | {mean:%d,}\n" + \
              "3rd Qu.| {TrdQu:%d,}\n" + \
              "Max    | {max:%d,}\n" + \
              "Total  | {totalLength:%d,}\n" + \
              "n50    | {n50:%d,}\n" + \
              "n90    | {n90:%d,}\n" + \
              "n95    | {n95:%d,}\n").replace("%d", str(space))
     
    #print stats
    print "Read Stats"#, json.dumps(getStats(readLengths), indent=4)
    print report.format(**stats)
    print "Bases Counted %d" % tot
    print "Average Accuracy %.2f" % (accuracy/cnt)
    print "Total Unmapped %d" % (unmapped)
    print "Percent Unmapped %.2f" % (unmapped/cnt)
    print
    print "Total Insertions %d" % insertions
    print "Average Insertions per Read %.2f" % (insertions/cnt)
    print "Percentage of errors Insertions %.2f" % (insertions/errCnt)
    print
    print "Total Deletions %d" % deletions
    print "Average Deletions per Read %.2f" % (deletions/cnt)
    print "Percentage of errors Deletions %.2f" % (deletions/errCnt)
    print
    print "Total Substitutions %d" % subs
    print "Average Substitutions per Read %.2f" % (subs/cnt)
    print "Percentage of errors Substitutions %.2f" % (subs/errCnt)
    print
    print "Total SoftClipped %d" % soft
    print "Average SoftClipped per Read %.2f" % (soft/cnt)
