# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 16:47:58 2019

@author: hunte
"""

import sys
if len(sys.argv) > 1:
    generalModelMetadataFile = sys.argv[1]
    
with open(generalModelMetadataFile, "r") as r:
    print(r.readline())
r.close()