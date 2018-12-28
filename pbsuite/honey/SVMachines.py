import sys, os, glob, copy, numpy, logging
from ctypes import *
#from svmutil import *
#from svm import *

"""
This folder holds different svmachines.
This program will help load any particular model as well as give information on it
"""

#SVMachines.py
class Machine():
    """
    Holds information and easy access to models
    """
    def __init__(self, svmDir):
        """
        svmDir is a directory containing all of the training information
        from a libsvm created svm model
        """
        self.name = os.path.dirname(svmDir).split(os.path.sep)[-1]
        self.dataPath = svmDir
        self.modelFile = os.path.join(svmDir, "model")
        self.rangeFile = os.path.join(svmDir, "range")
        self.trainFile = os.path.join(svmDir, "train")
        self.scaleFiles = (os.path.join(svmDir, "scale"),\
                           os.path.join(svmDir, "scale.out"),\
                           os.path.join(svmDir, "scale.png"))
        if self.isModelDir():
            self.loadModel()
            self.loadRange()
    
    def isModelDir(self):
        """
        Ensure that all of the files exist and I'm working with a
        usable model -- currently i only do the former
        """
        return (os.path.exists(self.modelFile) and \
                os.path.exists(self.rangeFile))
                #os.path.exists(self.trainFile) and \
                #os.path.exists(self.scaleFiles[0]) and \
                #os.path.exists(self.scaleFiles[1]) and \
                #os.path.exists(self.scaleFiles[2]))
    
    def loadModel(self):
       """
       svmutils.load_svm_model or something -- I probably won't need more work
       should probably monkey patch over the svmutils stuff
       """
       self.model = svm_load_model(self.modelFile)

    def loadRange(self):
        """
        x
        -1 1
        1 0 1
        2 0 1.611111
        3 -1 1.611111
        becomes 
        myscale = {1: [0, 1], 2: [0, 1.611111], 3:[-1, 1.611111]}
        """
        ret = {}
        ret = {}
        with open(self.rangeFile,'r') as fh:
            t = fh.readline(); t = fh.readline()
            for line in fh.readlines():
                data = line.strip().split()
                ret[int(data[0])] = map(float, data[1:])
        self.ranges = ret
     
    def predict(self, x):
        """
        predict for each of the objects. 
        x is a list of either dictionaries or lists
            [{1:0, 2:0, ...}, ...]
            [[0, 0, ...], ...]
          
        reuturn 2-tuple of [ (p_label, p_val), ... ]
        you usually want to call normalizeInstances x beforehand
        """
        p_label, p_acc, p_val = svm_predict( [0] * len(x), x, self.model )
        return p_label, p_val
        
    def nppredict(self, x):
        """
        predicts for the numpy array x. -- this goes a little quicker than
        Model.predict
        array structure should be (n, N) where n is the number of 
        measurements per feature and N is the number of features
        usually want to call normalize beforehand
        --- This is mainly me cherry picking from libsvm/python/svmutil.py and svm.py
            Right now it will only support models created using easy.py (i.e. no prediction 
            probablity and stuff)
        
        """
        m = self.model
        
        svm_type = m.get_svm_type()
        is_prob_model = m.is_probability_model()
        nr_class = m.get_nr_class()
        
        nr_classifier = nr_class*(nr_class-1)//2

        ret = numpy.zeros( len(x[0]) )
        
        index_range = range(1, len(x)+1)
        nelem = len(index_range) + 1
        
        dec_values = (c_double * nr_classifier)()
        for pos in xrange(len(x[0])):
            #apply across axis?
            data = (svm_node * nelem)()
            data[-1].index = -1
            #for each index
            for idx, j in enumerate(index_range):
                data[idx].index = j 
                data[idx].value = x[idx][pos]
            label = libsvm.svm_predict_values(m, data, dec_values)
            ret[pos] = label
            #ret[pos][1] =  prob_estimates[:nr_class]
            
        return ret
        
    def normalize(self, instances):
        """
        This works on numpy.arrays
        normalizes in place 
        """
        lower = -1.
        upper =  1
        
        ranges = self.ranges
        def getVal(a):
            """
            Mathmatics for normalizing
            This is pretty slow.
            """
            ret = numpy.zeros(len(a))
            for pos,i in enumerate(a):
                if i == minimum:
                    ret[pos] = lower
                elif i == maximum:
                    ret[pos] = upper
                else:
                    ret[pos] = lower + (upper-lower) * (i - minimum) / (maximum - minimum)
            return ret
        
        for i in xrange(len(ranges.keys())):
            attribute = i + 1
            minimum, maximum = ranges[attribute]
            instances[i] = numpy.apply_along_axis(getVal, 0, instances[i])
        
    def normalizeInstances(self, instances) :
        """
        normalize the instances sent to the model 
        This works for now, but I'll need it to work on numpy.arrays
        --- If I can get it out of the dictionary, I can supply the sparse matrix
            to svm_predict(y, x, m)
        """
        lower = -1.
        upper =  1
        
        ranges = self.ranges

        normalized_instances = copy.deepcopy(instances)
        if ranges == None :
            ranges_dict = dict()
        max_attribute = max(normalized_instances[0].keys())
        for attribute in range(1, max_attribute + 1) :  # we iterate on the attributes
            column = []
            for instance in normalized_instances:
                if attribute not in instance.keys():
                    instance[attribute] = 0.
                column.append(instance[attribute])
            if ranges != None :
                minimum = ranges[attribute][0]
                maximum = ranges[attribute][1]
            else :
                minimum = min(column)
                maximum = max(column)
                ranges_dict[attribute] = [minimum, maximum]
            for i in range(len(column)) :
                if column[i] == minimum :
                    column[i] = lower
                elif column[i] == maximum :
                    column[i] = upper
                else :
                    column[i] = lower + (upper-lower) * (column[i] - minimum) / (maximum - minimum)
            # Copying normalized values in memory
            for elem, instance in zip(column, normalized_instances):
                instance[attribute] = elem
        
        if ranges == None :
            return normalized_instances, ranges_dict
        else :
            return normalized_instances
    
def loadMachine(name = None):
    """
    Loads machines.
    returns the appropriate Machine
    returns a list of machines if name is None
    raises errors about machines not being able to load
    will fail and sys.exit(1) if name.svm can't be loaded
    """
    availMachines = []
    if name is None:
        #load all machines
        svmDirs = glob.glob(os.path.join(os.path.dirname(__file__), "svm", "*.svm"))
    else:
        svmDirs = glob.glob(os.path.join(os.path.dirname(__file__), "svm", "%s.svm" % name))
    
    for machine in svmDirs:
        if not os.path.isdir(machine):
            continue
        m = Machine(machine)
        if not m.isModelDir():
            logging.error("SVMachine %s couldn't be loaded" % name)
            if name is None:
                sys.exit(1)
        availMachines.append(Machine(machine))
        
    if name is not None:
        availMachines = availMachines[0]
    
    return availMachines
    

if __name__ == '__main__':
    setupLogging()
    loadMachine()
    

