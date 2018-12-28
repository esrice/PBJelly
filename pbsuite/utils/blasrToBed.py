#!/usr/bin/env python
import sys
from pbsuite.utils.FileHandlers import M4File, M5File

if __name__ == '__main__':
    try:
        fn = sys.argv[1]
    except:
        sys.stderr.write(("Error! Expected One Argument, " \
                          "an m4 or m5 alignment file\n"))
        exit(1)

    if fn.endswith('.m4'):
        file = M4File(sys.argv[1])
    elif fn.endswith('.m5'):
        file = M5File(sys.argv[1])
    else:
        print "Unrecognized File Type (expecting  .m4 or .m5)"
        exit(1)
    
    if len(sys.argv) == 3:
        out = open(sys.argv[2],'w')
    else:
        out = sys.stdout
    
    out.write("\n".join(map(lambda x: x.toBed(), file))+"\n")
