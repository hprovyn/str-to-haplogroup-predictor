# -*- coding: utf-8 -*-
"""
Created on Thu May 23 13:14:31 2019

@author: hunte
"""

theKits = {}
def readResultsFile(fil):
    with open(fil, 'r') as f:
        for line in f.readlines():
            splitLine = line.split(",")
            theid = splitLine[0]
            marker = splitLine[1]
            allele = splitLine[2]
            #print(theid, marker, allele)
            if theid not in theKits.keys():
                theKits[theid] = {"STRs": {}, "pos":[], "neg":[]}
            if allele[-1] == "+":
                theKits[theid]["pos"].append(marker)
            if allele[-1] == "-":
                theKits[theid]["neg"].append(marker)
            if allele[-1] != "+" and allele[-1] != "-":
                theKits[theid]["STRs"][marker] = allele

import os
def writeOutPosNegs(theid, outDir):
    thepath = os.path.join(outDir, theid)
    print('attempting to make dir', thepath)
    os.mkdir(thepath)
    posFile = os.path.join(outDir, theid, "pos")
    negFile = os.path.join(outDir, theid, "neg")
    strFile = os.path.join(outDir, theid, "str")
    
    with open(posFile, 'w') as f:
        f.write(",".join(theKits[theid]["pos"]))
    f.close()
    with open(negFile, 'w') as f:
        f.write(",".join(theKits[theid]["neg"]))
    f.close()
    with open(strFile, 'w') as f:
        kvs = []
        for k in theKits[theid]["STRs"]:
            kvs.append(k + "=" + theKits[theid]["STRs"][k])
        f.write(",".join(kvs))
    f.close()

import sys

if len(sys.argv) > 1:
    resultsFile = sys.argv[1]
    outDir = sys.argv[2]

readResultsFile(resultsFile)

import shutil
try:
    shutil.rmtree(outDir)
except:
    0
os.mkdir(outDir)

for i in theKits.keys():
    if len(theKits[i]["STRs"].keys()) > 12 and len(theKits[i]["pos"]) > 0:
        print(theKits[i])
        writeOutPosNegs(i, outDir)
    else:
        print(i,"ignored because not enough STRs or positive SNPs")