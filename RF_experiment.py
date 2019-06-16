# -*- coding: utf-8 -*-
"""
Created on Wed May 22 14:35:06 2019

@author: hunte
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 15:47:45 2018

@author: hunte
"""

from Common import CommonMethods
import sys
if len(sys.argv) > 1:
    trainFile = sys.argv[1]
    experimentOutputFile = sys.argv[2]
    modesIncluded = sys.argv[3]
    cutoff = int(sys.argv[4])
    iterations = int(sys.argv[5])

CommonMethods.experiment(trainFile, modesIncluded, cutoff, experimentOutputFile, iterations)