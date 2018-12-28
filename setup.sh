#/bin/bash

#If you use a virtual env - source it here
source activate py27

#This is the path where you've install the suite.
export SWEETPATH=/home/anscgenomics/esrice/software/PBSuite_15.8.24
#for python modules
export PYTHONPATH=$PYTHONPATH:$SWEETPATH
#for executables
export PATH=$PATH:$SWEETPATH/bin/
