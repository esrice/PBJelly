#!/hgsc_software/python/python-2.6/bin/python
import sys, math

def getStats(seqLengths):
    data = {}

    seqLengths.sort(reverse=True)
    
    data["numItems"] = len(seqLengths)
    data["itemSum"] = sum(seqLengths)
    tl = data["itemSum"]
    n50_mark = data["itemSum"] * .5
    n90_mark = data["itemSum"] * .90
    n95_mark = data["itemSum"] * .95
    
    data["n50"] = None
    data["n50_gt_count"] = None
    data["n90"] = None
    data["n90_gt_count"] = None
    data["n95"] = None
    data["n95_gt_count"] = None
    basesSeen = 0
    
    for pos,n in enumerate(seqLengths):
        basesSeen += n
        if data["n50"] is None and basesSeen > n50_mark:
            data["n50"] = n
            data["n50_gt_count"] = pos
        if data["n90"] is None and basesSeen > n90_mark:
            data["n90"] = n
            data["n90_gt_count"] = pos
        if data["n95"] is None and basesSeen > n95_mark:
            data["n95"] = n
            data["n95_gt_count"] = pos
            break
    #may not have gaps
    if data["numItems"] == 0:
        return data
    data["min"] = seqLengths[-1]
    data["FstQu"] = seqLengths[ int(math.floor(data["numItems"]*.75)) ]
    median = data["numItems"]*.50
    data["median"] = int( (seqLengths[ int(math.floor(median)) ] + \
                           seqLengths[ int(math.floor(median)) ]) / 2)
    data["mean"] = data["itemSum"]/data["numItems"]
    data["TrdQu"] = seqLengths[ int(math.floor(data["numItems"]*.25)) ] 
    data["max"] = seqLengths[0]

    return data

def run(data):
    """
    list of numbers - can be  a string if you want
    """
    data = map(float, data)
    ret = getStats(data)
    
    outputOrder = ["itemSum",
                   "numItems",
                   "min",
                   "FstQu",
                   "mean",
                   "median",
                   "n50",
                   "n50_gt_count",
                   "TrdQu",
                   "n90",
                   "n90_gt_count",
                   "n95", 
                   "n95_gt_count",
                   "max"]

    for key in outputOrder:
        print "{0}\t{1:.2f}".format(key, ret[key])

if __name__ == '__main__':
    run(sys.stdin.read().strip().split('\n'))
    
