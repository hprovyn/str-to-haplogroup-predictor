# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 18:26:41 2019

@author: hunte
"""

from Common import CommonMethods
import sys
if len(sys.argv) > 1:
    dataDir = sys.argv[1]
    panelHierarchyFile = sys.argv[2]
    policyFileStem = sys.argv[3]
    modelPickleFileStem = sys.argv[4]
    percentMissingSTRThreshold = int(sys.argv[5])/100
    sampleId = sys.argv[6]

CommonMethods.predict(dataDir, sampleId, panelHierarchyFile, policyFileStem, modelPickleFileStem, percentMissingSTRThreshold)