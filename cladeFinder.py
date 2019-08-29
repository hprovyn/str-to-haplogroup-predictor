# -*- coding: utf-8 -*-
"""
Created on Fri Sep  7 15:23:13 2018

@author: hunte
"""

"""
python cladeFinder.py YFull_YTree_v6.02__20180402.tree.json pos negatives output

create negatives
bcftools filter -O z -i '(GT=="1/1" && AA==ALT) || (GT=="0/0" && AA=REF)' chrY_cleaned_1_hg38.vcf | bcftools query -f '%ID,' > negatives

create positives
bcftools query -f '%ID,' chrY_derived_1_hg38.vcf.gz > pos

"""


import json
import sys

#TODO this may need updating
toIgnore = ["PF129", "Z2533", "S6868", "BY2285", "YP2229", "YP2250", "YP2228", "YP2129", "YP1838", "YP1807", "YP1841", "YP1740", "Y17293", "PF6234",  "L132.2"]

if len(sys.argv) > 1:
    treeFile = sys.argv[1]
    tabixFilePath = sys.argv[2]
    outputFile = sys.argv[3]
    maxThreads = int(sys.argv[4])
 
def parseFile(file):
    fr = open(file, 'r')
    parsed = fr.readline().split(",")
    fr.close()
    return parsed

hierarchy = {}
childMap = {}
snps = {}

def parseTreeJSON(fil):
    thefile = open(fil)
    root = json.load(thefile)
    thefile.close()
    recurseTreeJson(root, hierarchy, snps)
    return (root["id"], hierarchy, snps)

def parseSNPsString(snpsString):
    thesnps = set([])
    for snps in snpsString.split(", "):
        for snp in snps.split("/"):
            thesnps.add(snp)
    return thesnps
            
def recurseTreeJson(node, hierarchy, snps):
    if "children" in node:
        childMap[node["id"]] = []
        for child in node["children"]:
            childMap[node["id"]].append(child["id"])
            hierarchy[child["id"]] = node["id"]
            snps[child["id"]] = parseSNPsString(child["snps"])
            recurseTreeJson(child, hierarchy, snps)
                
def getChildren(clade, childParents):
#    children = []
#    for child in childParents:
#        if childParents[child] == clade:
#            children.append(child)
#
    if clade in childMap:
        return childMap[clade]
    else:
        return []

def isInChildrenThisLevel(clade, positives, childParents):
    children = getChildren(clade, childParents)
    inChildren = []
    for child in children:
        if any(snp in positives for snp in snps[child]):
            inChildren.append(child)
    return inChildren

def recurseDownTreeUntilFirstHits(clade, positives, childParents):
    posChildrenThisLevel = isInChildrenThisLevel(clade, positives, childParents)
    for child in getChildren(clade, childParents):
        if child not in posChildrenThisLevel:
            childResult = recurseDownTreeUntilFirstHits(child, positives, childParents)
            for cres in childResult:
                posChildrenThisLevel.append(cres)
    return posChildrenThisLevel

def removeDuplicates(arr): 

    n = len(arr)
    # Return, if array is  
    # empty or contains 
    # a single element 
    if n == 0 or n == 1: 
        return n 
  
    temp = list(range(n)) 
  
    # Start traversing elements 
    j = 0; 
    for i in range(0, n-1): 
  
        # If current element is 
        # not equal to next 
        # element then store that 
        # current element 
        if arr[i] != arr[i+1]: 
            temp[j] = arr[i] 
            j += 1
  
    # Store the last element 
    # as whether it is unique 
    # or repeated, it hasn't 
    # stored previously 
    temp[j] = arr[n-1] 
    j += 1
      
    # Modify original array 
    for i in range(0, j): 
        arr[i] = temp[i] 
  
    return arr

def refineHitsRecursively(sequences, positives, childParents, solutions):
    for sequence in sequences:
        refinedResults = recurseDownTreeUntilFirstHits(sequence[-1], positives, childParents)
        if len(refinedResults) == 0:
            solutions.append(sequence)
        else:
            print(sequence, refinedResults)
            for refRes in refinedResults:
                #print(sequence, refRes)
                seqCopy = sequence[:]
                seqCopy.append(refRes)
                refineHitsRecursively([seqCopy], positives, childParents, solutions)               

def recurseDownTree(positives, childParents, solutions):
    sequences = recurseDownTreeUntilFirstHits("", positives, childParents)
    newSequences = []
    for sequence in sequences:
        newSequences.append([sequence])
    refineHitsRecursively(newSequences, positives, childParents, solutions)

def getTotalSequence(clade, hierarchy):
    sequence = [clade]
    thisClade = clade
    while thisClade in hierarchy:
        thisClade = hierarchy[thisClade]
        sequence.append(thisClade)
    return sequence[:-1]
    
def getScore(sequence, totalSequence):
    return float(len(sequence)) / len(totalSequence)

def printSolutions(solutions):
    for solution in solutions:
        print(" ".join(solution), getScore(solution))

def getConflicts(sequence, negatives, hierarchy):
    conflictingNegatives = []
    for hg in sequence:
        if any(snp in negatives for snp in snps[hg]):
            conflictingNegativeSnps = ""
            for snp in snps[hg]:
                if snp in negatives:
                    conflictingNegativeSnps += " " + snp
            conflictingNegatives.append(hg + " @" + conflictingNegativeSnps + ";")
    return conflictingNegatives

def getWarnings(sequence, negatives, hierarchy):
    messages = []
    conflicts = getConflicts(sequence, negatives, hierarchy)
    for conflict in conflicts:
        messages.append(" " + conflict)
    return messages

def getWarningsConf(conflicts):
    messages =[]
    for conflict in conflicts:
        messages.append(" " + conflict)
    return messages

def unconflictedPercent(sequence, conflicts):
    return (float(len(sequence)) - len(conflicts)) / len(sequence)

import numpy as np
def getPathScores(fullSequence, confirmed, negatives, positives, conflicts):
    scores = []
    for thing in fullSequence:
        if thing in confirmed:
            negs = len(negatives.intersection(set(snps[thing])))
            poses = len(positives.intersection(set(snps[thing])))
            toappend = 0
            if poses + negs > 0:
                toappend = float(poses) / (poses + negs)
            scores.append(toappend)
        else:
            if thing in conflicts:
                scores.append(0)
            else:
                scores.append(0.5)
    return scores

def isBasal(clade, negatives, positives, hierarchy):
    basal = False
    children = getChildren(clade, hierarchy)
    if len(children) > 0:
        basal = True
        for child in children:
            isNeg = len(negatives.intersection(set(snps[child]))) > 0
            isPos = len(positives.intersection(set(snps[child]))) > 0
            if isNeg and not isPos:
                basal = basal and True
    return basal        

def writeFile(scoredSolutions):
    w = open(outputFile,"w+")
    w.write("#Haplogroup\tPath\tPurity\tDepth\tScore\tConflicting Negative SNPs\n")
    for scoredSolution in scoredSolutions:
        w.write(scoredSolution[1] + "\t" + " > ".join(scoredSolution[0]) + "\t" + str(scoredSolution[2]) + "\t" + str(scoredSolution[3]) + "\t" + str(scoredSolution[4]) + "\t" + " ".join(scoredSolution[5]))
        w.write("\n")
    w.close()

def writeSampleCladesFile(cladeMap):
    w = open(outputFile,"w+")
    for sample in cladeMap:
        w.write(",".join([sample,cladeMap[sample]]) + "\n")
    
    
from operator import itemgetter
def getRankedSolutions(positives, negatives, hierarchy):
    solutions = []
    recurseDownTree(positives, hierarchy, solutions)
    scoredSolutions = []
    print(solutions)
    uniqueSolutions = removeDuplicates(solutions)
    print(uniqueSolutions)
    for solution in solutions:
        lastChainMoreNegThanPos = True
        removed = 0
        while lastChainMoreNegThanPos and removed < len(solution):
            totalSequence = getTotalSequence(solution[-1 - removed], hierarchy)
            totalSequence.reverse()
            conflicts = getConflicts(totalSequence, negatives, hierarchy)
            scores = getPathScores(totalSequence, solution, negatives, positives, conflicts)
            clade = solution[-1 - removed]
            if isBasal(clade, negatives, positives, hierarchy):
                clade = clade + "*"
            scoredSolutions.append([totalSequence, clade, np.average(scores), np.sum(scores), np.average(scores) * np.sum(scores), getWarningsConf(conflicts)])
            removed = removed + 1
            if scores[-1] > 0.5:
                lastChainMoreNegThanPos = False
            #else:
                #print(totalSequence, scores[-1], scores)
            
        #print(scoredSolutions[-1])
    scoredSolutions = sorted(scoredSolutions, key=itemgetter(4), reverse=True)
    
    return scoredSolutions

import time
start = time.time()

a = parseTreeJSON(treeFile)
end = time.time()
print('load tree elapsed time,' + str(round((end-start)/60,2)) + " minutes\n")

hierarchy = a[1]
snps = a[2]

cladeMap = {}

import tabix

#TODO get unique column values tabix query?

tb = tabix.open(tabixFilePath)
idresults = tb.querys("id:1-9999999")
ids = []
for theid in idresults:
    ids.append(theid[1])

class ParallelCladeFind():
        
    def findClade(self, tb, sampleId, hierarchy, toIgnore, queue):
        posResults = tb.querys("pos:" + sampleId + "-" + sampleId)
        positives = []
        for p in posResults:
            positives.append(p[3])
        positives = set(positives)
        negResults = tb.querys("neg:" + sampleId + "-" + sampleId)
        negatives = []
        for n in negResults:
            negatives.append(n[3])
        negatives = set(negatives)
        for ign in toIgnore:
            if ign in positives:
                positives.remove(ign)
            if ign in negatives:
                negatives.remove(ign)
        b = getRankedSolutions(positives, negatives, hierarchy)
        if len(b) > 0:
            print(sampleId, ' computed as ', b[0][1])
            queue.put((sampleId, b[0][1]))
        else:
            print(sampleId, ' unable to compute')
            queue.put((sampleId, "?"))

import multiprocessing


allStart = time.time()

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

idchunks = chunks(ids, maxThreads)

cladeMap = {}
for idchunk in idchunks:
    processes = []
    queue = multiprocessing.Queue()
    for sampleId in idchunk:
        pcf = ParallelCladeFind()
        p = multiprocessing.Process(target=pcf.findClade, args=(tb, sampleId, hierarchy, toIgnore, queue))
        processes.append(p)
        p.start()

    rets = []
    for p in processes:
        ret = queue.get()
        rets.append(ret)

    for p in processes:
        p.join()

    for ret in rets:
        cladeMap[ret[0]] = ret[1]



allEnd = time.time()
print('clade finder executed on ' + str(len(ids)) + ' ' + str(round((allEnd-allStart)/60,2)) + " minutes\n")
print(str(round((allEnd-allStart)/len(ids),2)) + " seconds per sample\n")

writeSampleCladesFile(cladeMap)

