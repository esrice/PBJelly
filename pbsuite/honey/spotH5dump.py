#!/usr/bin/env python
import sys
import h5py

#view with chr:start-end
if __name__ == '__main__':
    h5 = h5py.File(sys.argv[1])
    print "chrom\tposition\t"+"\t".join(h5.attrs["columns"])
    for chrom in h5.keys():
        for pos,i in enumerate(h5["/%s/data" % chrom].value.transpose().tolist()):
            print chrom + "\t" + str(pos) + '\t' + "\t".join(map(str, i)) 
