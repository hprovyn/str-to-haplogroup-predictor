# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 18:26:41 2019

@author: hunte
"""

from Common import CommonMethods
import sys
if len(sys.argv) > 1:
    panelHierarchyFile = sys.argv[1]
    policyFileStem = sys.argv[2]
    modelPickleFileStem = sys.argv[3]
    percentMissingSTRThreshold = float(sys.argv[4])/100
    strAlleleString = sys.argv[6]
    haplogroupClassConfigPath = sys.argv[5]

    
CommonMethods.predict(strAlleleString, panelHierarchyFile, policyFileStem, modelPickleFileStem, percentMissingSTRThreshold, haplogroupClassConfigPath)
