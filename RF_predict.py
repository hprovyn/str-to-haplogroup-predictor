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
    sampleId = sys.argv[3]

CommonMethods.predict(trainFile, dataDir, sampleId)