#!/usr/bin/env python
import sys, h5py, argparse

from pbsuite.utils.setupLogging import *

USAGE = """\
Merge multple spots.h5 files.
"""

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=USAGE, \
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--output", required=True, \
                        help="Output spot.h5 file name")
    parser.add_argument("spots", metavar="SPOTH5", type=str, nargs="+",\
                        help="Spot.h5 files to merge")
    args = parser.parse_args()
    
    output = h5py.File(args.output, 'w')
    loadCols = True
    for fn in args.spots:
        print fn
        newData = h5py.File(fn,'r')
        
        if loadCols:
            output.attrs["columns"] = newData.attrs["columns"]
            loadCols = False
            
        for group in newData.keys():
            if group not in output.keys():
                myGroup = output.create_group(group)
                #could put a check here that attrs are the same...
                for attr in newData[group].attrs.keys():
                    output.attrs[attr] = newData[group].attrs[attr]
                data = myGroup.create_dataset("data", data=newData[group]["data"], compression="gzip")
                print "setup", group
            else:
                print "merge", group
                data = newData[group]["data"].value + output[group]["data"].value
                output[group]["data"].write_direct(source=data)
    
    output.close()
    """
    for all other args
        
        open other
        for group in other
            if group not in output
                create it
            merge groups output[group][data] = other[group][data]
    
    done
    """
