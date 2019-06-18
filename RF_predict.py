# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 18:26:41 2019

@author: hunte
"""

from Common import CommonMethods
import sys
if len(sys.argv) > 1:
    trainFile = sys.argv[1]
    dataDir = sys.argv[2]
    panelHierarchyFile = sys.argv[3]
    policyFileStem = sys.argv[4]
    modelPickleFileStem = sys.argv[5]
    sampleId = sys.argv[6]

CommonMethods.predict(trainFile, dataDir, sampleId, panelHierarchyFile, policyFileStem, modelPickleFileStem)