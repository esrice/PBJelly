OLCAssembly - A shortcut to calling pacbio's Allora Assembler

Polish - A quick consensus caller

PBJNovo - DeNovo Assembler of PacBio only reads

run order
1) chunkyBlasr.py
	creates aligns/chunk*.m4
2) MakeReciprocal.py
	Creates B->A since chunkyBlasr.py saves times and only makes A->B
3) sortChunks.sh
	Put the reads on order for easy parsing
4) MakeOverlapTable.py
	Creates your Overlap Table
todo:
	need layouts and then to include polish
