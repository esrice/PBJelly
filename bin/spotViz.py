#!/usr/bin/env python
import h5py, sys, numpy, random
from pbsuite.honey.HSpots import *
import matplotlib.pyplot as plt


#transposon
#start = 1976517 - buffer
#end   = 1977304 + buffer
#tandem duplication
#start = 1096192 - buffer
#end   = 1096825 + buffer
#tandem deletion
#start = 4294209 - buffer
#end   = 4294314 + buffer    
#p-element inversion
#appx 8192387
#start, end = 304417-buffer, 304460+buffer #tails
#start, end = 1456604 - buffer, 1456650 + buffer
#start, end = 3255896 - buffer, 3255932 + buffer

avgWindow = numpy.ones(50, dtype=numpy.float16) / 50

def makeTransformPlots(key, data, start, end, buffer, binsize, fignum, normalize=True):
    orig = numpy.convolve(data, avgWindow, "same") #smooth
    if normalize:
        norm = orig / cov
    else:
        norm = orig
        
    slopeT = numpy.convolve(norm, slopWindow, "same") #slope transform
    abs = numpy.abs(slopeT) # absolute value
    smo = numpy.convolve(abs, avgWindow, "same") #smooth
    space = numpy.convolve(smo * orig, avgWindow, "same") # Take back to absolute space and smooth again
    win = range(start-buffer, end+buffer)
    plt.figure(fignum)
    
    plt.subplot(711)
    plt.plot(win, data, 'b-')
    plt.axhline(1, color='k'); plt.axvline(start, color='k'); plt.axvline(end, color='k')
    
    plt.subplot(712)
    plt.plot(win, orig, 'b-')
    plt.axhline(1, color='k'); plt.axvline(start, color='k'); plt.axvline(end, color='k')
    
    plt.subplot(713)
    plt.plot(win, norm, 'b-')
    plt.axhline(1, color='k'); plt.axvline(start, color='k'); plt.axvline(end, color='k')
    
    plt.subplot(714)
    plt.plot(win, slopeT, 'b-')
    plt.axhline(1, color='k'); plt.axvline(start, color='k'); plt.axvline(end, color='k')
    
    plt.subplot(715)
    plt.plot(win, abs, 'b-')
    plt.axhline(1, color='k'); plt.axvline(start, color='k'); plt.axvline(end, color='k')
    
    plt.subplot(716)
    plt.plot(win, smo, 'b-')
    plt.axhline(1, color='k'); plt.axvline(start, color='k'); plt.axvline(end, color='k')
    
    plt.subplot(717)
    plt.plot(win, space,  'b-')
    plt.axhline(1, color='k'); plt.axvline(start, color='k'); plt.axvline(end, color='k')
    
    plt.show()
    plt.savefig("transform_%s.png" % key)

def makeLinePlots(dataOrig, start, end, offset):
    plt.figure()
    print start, end, offset
    print len(dataOrig)
    cov  = numpy.convolve(data[COV], avgWindow, "same") 
    print "cov", numpy.max(numpy.abs(data[COV][10:-10])), numpy.mean(cov)
    
    
    insR, imu, isd = preprocessSignal(data[INS], data[COV])
    insS = signalTransform(insR)
    
    delR, dmu, dsd = preprocessSignal(data[DEL], data[COV])
    delS = signalTransform(delR)
    
    win = range(start, end)
    
    covRg = plt.plot(win, cov,  'k-')
    misRg = plt.plot(win, misR * cov, 'r-')
    insRg = plt.plot(win, insR * cov, 'g-')
    delRg = plt.plot(win, delR * cov, 'b-')
    
    
    ticks = range(start, end, (end-start)/6)#[:-1]
    labels = range(offset, offset + (end-start)+1, (end-start)/6)
    print ticks, labels
    plt.xticks(ticks, ticks, horizontalalignment="left", rotation=17)
    plt.xlabel("position")
    plt.ylabel("rate")
    plt.legend([covRg, insRg, delRg], ["COV", "INS", "DEL"],)
    plt.suptitle("%d bp sv (%d - %d)" % (end - start, start, end))
    plt.show()
    plt.savefig("rates.png")

    plt.figure()
    
    insSg = plt.plot(win, insS, 'g-')
    delSg = plt.plot(win, delS, 'b-')
    
    ticks = range(start, end+1, (end-start)/6)#[:-1]
    #labels = range(offset, offset + (end-start)+1, (end-start)/6)
    plt.xticks(ticks, ticks, horizontalalignment="center", rotation=17)
    plt.xlabel("position")
    plt.ylabel("signal")
    plt.legend([insSg, delSg], ["INS", "DEL"])
    plt.suptitle("%d bp sv (%d - %d)" % (end - start, start, end))
    plt.show()
    plt.savefig("signals.png")

    
    
if __name__ == '__main__':
    h5    = h5py.File(sys.argv[1], 'r')
    cols  = h5.attrs["columns"]
    key   = sys.argv[2]
    offset = int(sys.argv[3])
    start = int(sys.argv[3]) - h5[key].attrs["start"]
    end   = int(sys.argv[4]) - h5[key].attrs["start"]
    
    #buffer = int((start-end) *.1)
    buffer = 100
    start  = max(0, start-buffer)
    end    = min(h5[key].attrs["end"], end+buffer)
    data   = h5[key]["data"][: , start:end]
    h5.close()
    makeKernals()
    makeLinePlots(data, start, end, offset)
    
