import sys


from difflib import SequenceMatcher

USAGE= "%prog A small assembler that quickly guesses at consensus of small, completely overlapping sequences"""

ss = ["CGCGTACGCG",
     ["
a = "dsa jld lal"
b = "dsajld kll"
c = "dsc jle kal"
d = "dsd jlekal"

ss = [a,b,c,d]

s = SequenceMatcher()

for i in range(len(ss)):
    x = ss[i]
    s.set_seq1(x)
    for j in range(i+1,len(ss)):    
        y = ss[j]
        s.set_seq2(y)   
        print
        print s.ratio()
        print s.get_matching_blocks()
