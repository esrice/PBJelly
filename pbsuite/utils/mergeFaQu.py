#!/usr/bin/env python
import sys

from pbsuite.utils.FileHandlers import mergeFastaQual

if __name__ == '__main__':
    f = mergeFastaQual(sys.argv[1], sys.argv[2])
    for i in f:
        print f[i].toString(),
