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
                theKits[theid] = {"str": {}, "pos":[], "neg":[]}
            if allele[-1] == "+":
                theKits[theid]["pos"].append(marker)
            if allele[-1] == "-":
                theKits[theid]["neg"].append(marker)
            if allele[-1] != "+" and allele[-1] != "-":
                theKits[theid]["str"][marker] = allele

def writeOutPosNegs(outFile):
        
    with open(outFile, 'w') as f:
        for theid in theKits:
            pos = theKits[theid]["pos"]
            neg = theKits[theid]["neg"]
            strs = theKits[theid]["str"]
            if hasEnoughInfoToProceed(pos, strs):
                for p in pos:
                    f.write("\t".join(["pos", theid, theid, p, "."]) + "\n")
                for n in neg:
                    f.write("\t".join(["neg", theid, theid, n, "."]) + "\n")
                for marker in strs:
                    f.write("\t".join(["str", theid, theid, marker, strs[marker]]) + "\n")
            else:
                print(theid,"ignored because not enough STRs or positive SNPs")
    f.close()
    
def hasEnoughInfoToProceed(pos, strs):
    if len(strs) > 15 and len(pos) > 0:
        return True
    else:
        return False
    
import sys

if len(sys.argv) > 1:
    resultsFile = sys.argv[1]
    outFile = sys.argv[2]

readResultsFile(resultsFile)
writeOutPosNegs(outFile)        