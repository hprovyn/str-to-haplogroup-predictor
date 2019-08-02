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
    panelHierarchyFile = sys.argv[3]
    policyFileStem = sys.argv[4]
    modelPickleFileStem = sys.argv[5]
    utilityUnderSpecificityError = int(sys.argv[6])
    utilityOverSpecificityError=int(sys.argv[7])
    utilityCompletelyWrongError=int(sys.argv[8])
    utilityCorrect=int(sys.argv[9])
    if utilityUnderSpecificityError == 0 and utilityOverSpecificityError == 0 and utilityCompletelyWrongError == 0 and utilityCorrect == 0:
        utilityWeights = None
    else:
        utilityWeights = [utilityUnderSpecificityError, utilityOverSpecificityError, utilityCompletelyWrongError, utilityCorrect]
    experimentMapFileStem = sys.argv[10]

if __name__ ==  '__main__':
    CommonMethods.experimentErrorPolicy(trainFile, experimentOutputFile, panelHierarchyFile, policyFileStem, modelPickleFileStem, utilityWeights, experimentMapFileStem)