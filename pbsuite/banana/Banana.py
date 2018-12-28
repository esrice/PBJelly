import sys, json

from pbsuite.utils.FileHandlers import FastaFile, M5File
from pbsuite.utils.CommandRunner import exe
from pbsuite.jelly.Support import *

import networkx as nx

USAGE = "%prog - A Simple Local Assembler"

"""
Gameplan --- OLCAssembler for not so far spans

A) Blasr all vs all
B) Turn overlaps into a graph
    I'm going to weigh edges by similarity
C) Find easiest path through graph (build in find node to node path)
D) Turn these overlaps into a b.s. backbone
E) Map reads back over the backbone
F) Run Polish -- and I'm going to trim the ends since it's local
"""
def blasr(query, target, nproc = 1, outname = "out.m5"):
    """
    Simple overlapper
    """
    r,o,e = exe(("blasr %s %s -m 5 -bestn 200 -nCandidates 200 -minMatch 12 "
                 "-affineExtend 3 -nproc %d -noSplitSubreads -out %s -maxScore -1000") % \
                 (query, target, nproc, outname))

def m5ToOvlGraph(readNames, fileName):
    """
    Create the graph
    """
    connector = AlignmentConnector()
    alignments = M5File(fileName)
    graph = nx.Graph()       
    

    filt = []
    #get only the single best alignment between any two reads
    fdict = {}
    for align in alignments:
        if align.qname == align.tname:
            continue
        name = [align.qname, align.tname]
        name.sort()
        name = ":".join(name)
        if name in fdict:
            if align.score < fdict[name].score:
                fdict[name] = align
        else: 
            fdict[name] = align
    
    alignments = fdict.values()
    #make edges for all overlaps
    for align in alignments:
        if align.qname == align.tname:
            continue
        extend = connector.extendsTarget(align)
        align.support = extend
        if extend != SUPPORTFLAGS.none:
            graph.add_edge(align.qname, align.tname, data = align)
    
    return graph
    
def ovlSimplify(graph):
    """
    Find the most continuous path through and return the subgraphs
    I'll find the best score extending whichver end through all of the reads
    """
    def getStrongestEdge(node, end, used):
        """
        Gets the name of the strongest edge this node has and returns the neighboring node's name
        """
        best = 0
        bestN = None
        for edge in s.edges(node):
            align = s.edge[edge[0]][edge[1]]["data"]
            if align.tname != node or align.qname in used:
                continue
            if align.support in [end, SUPPORTFLAGS.span]:
                if align.score < best:
                    best = align.score
                    bestN = align
        return bestN
        
    subG = nx.connected_component_subgraphs(graph)
    paths = []
    for s in subG:
        for cnode in s:
            path = []
            used = []
            direction = SUPPORTFLAGS.left
            node = cnode
            while True:
                used.append(node)
                next = getStrongestEdge(node, direction, used)
                if next is None:    
                    break
                path.insert(0, next)
                if next.qstrand == '-':
                    if direction == SUPPORTFLAGS.left:
                        direction = SUPPORTFLAGS.right
                    elif direction == SUPPORTFLAGS.right:
                        direction = SUPPORTFLAGS.left
                node = next.qname
            node = cnode
            direction = SUPPORTFLAGS.right
            while True:
                used.append(node)
                next = getStrongestEdge(node, direction, used)
                if next is None:    
                    break
                path.append( next )
                if next.qstrand == '-':
                    if direction == SUPPORTFLAGS.left:
                        direction = SUPPORTFLAGS.right
                    elif direction == SUPPORTFLAGS.right:
                        direction = SUPPORTFLAGS.left
                node = next.qname
            paths.append(path)
    for p in paths:
        if len(p) == 0:
            continue
        for i in p:
            if not i.tname.startswith("ref"):
                print i.tname.split('/')[1],i.tstrand,"\t",
            else:
                print i.tname,i.tstrand,"\t",
        print
    
if __name__ == '__main__':
    reads = sys.argv[1]
    fasta = FastaFile(reads)
    #blasr(reads, reads, 4)
    ovl = m5ToOvlGraph(fasta.keys(), "out.m5")
    ovlSimplify(ovl)
    nx.write_gml(ovl, "ovl.gml")
