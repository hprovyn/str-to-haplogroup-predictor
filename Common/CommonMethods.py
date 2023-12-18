# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 14:41:59 2019

@author: hunte
"""

import numpy as np

def getSTRLabelsFromSets(modesIncluded):
    a_group = ["DYS391","DYS389I","DYS437","DYS439","DYS389II","DYS438","DYS426","DYS393","YCAII","DYS390","DYS385","Y-GATA-H4","DYS388","DYS447","DYS19","DYS392"]
    b_group = ["DYS458","DYS455","DYS454","DYS464","DYS448","DYS449","DYS456","DYS576","CDY","DYS460","DYS459","DYS570","DYS607","DYS442"]
    c_group = ["DYS728","DYS723","DYS711","DYR76","DYR33","DYS727","DYR157","DYS713","DYS531","DYS578","DYF395","DYS590","DYS537","DYS641","DYS472","DYF406S1","DYS511","DYS557","DYS490","DYS446","DYS481","DYS413","DYS534","DYS450","DYS425","DYS594","DYS444","DYS520","DYS436","DYS565","DYS572","DYS617","DYS568","DYS487","DYS640","DYS492"]
    d_group = ["DYR112","DYS518","DYS614","DYS626","DYS644","DYS684","DYS710","DYS485","DYS632","DYS495","DYS540","DYS714","DYS716","DYS717","DYS505","DYS556","DYS549","DYS589","DYS522","DYS494","DYS533","DYS636","DYS575","DYS638","DYS462","DYS452","DYS445","Y-GATA-A10","DYS463","DYS441","Y-GGAAT-1B07","DYS525","DYS712","DYS593","DYS650","DYS532","DYS715","DYS504","DYS513","DYS561","DYS552","DYS726","DYS635","DYS587","DYS643","DYS497","DYS510","DYS434","DYS461","DYS435"]
    ftdna = {}
    ftdna["12"] = ["DYS393","DYS390","DYS19","DYS391","DYS385","DYS426","DYS388","DYS439","DYS389I","DYS392","DYS389II"]
    ftdna["25"] = ftdna["12"] + ["DYS458","DYS459","DYS455","DYS454","DYS447","DYS437","DYS448","DYS449","DYS464"]
    ftdna["37"] = ftdna["25"] + ["DYS460","Y-GATA-H4","YCAII","DYS456","DYS607","DYS576","DYS570","CDY","DYS442","DYS438"]
    ftdna["67"] = ftdna["37"] + ["DYS531","DYS578","DYF395","DYS590","DYS537","DYS641","DYS472","DYF406S1","DYS511","DYS425","DYS413","DYS557","DYS594","DYS436","DYS490","DYS534","DYS450","DYS444","DYS481","DYS520","DYS446","DYS617","DYS568","DYS487","DYS572","DYS640","DYS492","DYS565"]
    ftdna["111"] = ftdna["67"] + ["DYS710","DYS485","DYS632","DYS495","DYS540","DYS714","DYS716","DYS717","DYS505","DYS556","DYS549","DYS589","DYS522","DYS494","DYS533","DYS636","DYS575","DYS638","DYS462","DYS452","DYS445","Y-GATA-A10","DYS463","DYS441","Y-GGAAT-1B07","DYS525","DYS712","DYS593","DYS650","DYS532","DYS715","DYS504","DYS513","DYS561","DYS552","DYS726","DYS635","DYS587","DYS643","DYS497","DYS510","DYS434","DYS461","DYS435"]
    strLabelsUsed = []    
    strs = []
    dubSTRs = []
    quadSTRs = []
    if "ftdna_" in modesIncluded:
        ftdnastrs = modesIncluded.replace("ftdna_","")
        strs = ftdna[ftdnastrs]
        allDubs = ["DYS385","YCAII","DYS459","CDY","DYF395","DYS413"]
        allQuads = ["DYS464"]
        for thestr in strs:
            if thestr in allDubs:
                dubSTRs += [thestr]
            if thestr in allQuads:
                quadSTRs += [thestr]
    else:
        if "a" in modesIncluded:
            strs = a_group
            dubSTRs = ["DYS385","YCAII"]
        if "b" in modesIncluded:
            strs += b_group
            dubSTRs += ["DYS459","CDY"]
            quadSTRs.append("DYS464")
        if "c" in modesIncluded:
            strs += c_group
            dubSTRs += ["DYF395", "DYS413"]
        if "d" in modesIncluded:
            strs += d_group
        
                
    for toremove in dubSTRs + quadSTRs:
        strs.remove(toremove)

    
    for thestr in strs:
        strLabelsUsed.append(thestr)
    
    for thestr in dubSTRs:
        strLabelsUsed.append(thestr + 'a')
        strLabelsUsed.append(thestr + 'b')
        strLabelsUsed.append(thestr + '-length')
        
    for thestr in quadSTRs:
        strLabelsUsed.append(thestr + 'a')
        strLabelsUsed.append(thestr + 'b')
        strLabelsUsed.append(thestr + 'c')
        strLabelsUsed.append(thestr + 'd')
        strLabelsUsed.append(thestr + '-length')
    print(modesIncluded, strs, dubSTRs, quadSTRs)
    return (strs, dubSTRs, quadSTRs)
    
def downstream(clade, hierarchy):
    children = [x for x in hierarchy if hierarchy[x] == clade]        
    downChilds = []
    for child in children:
        downChilds = np.append(downChilds, downstream(child, hierarchy))
 
    return np.append(downChilds, children)


def downstreamAndNotBelowNegative(clade, negatives, hierarchy):
    belowNegatives = []
    for neg in negatives:
        belowNegatives = np.append(belowNegatives, downstream(neg, hierarchy))
        belowNegatives = np.append(belowNegatives,[neg])
    alldowns = downstream(clade, hierarchy)
    downnotneg = []
    print(belowNegatives)
    for down in alldowns:
        if down not in belowNegatives:
            downnotneg.append(down)
    return downnotneg

def is_float(value):
    try:
        float(value)
        return True
    except:
        return False
    
    
def convertSTRrowToMap(ks, sampleIdx, strs, dubSTRs, quadSTRs):
    strmap = {}
    for thestr in set(strs) | set(dubSTRs) | set(quadSTRs):
        if ks[thestr][sampleIdx] != 0:
            strmap[thestr] = ks[thestr][sampleIdx]
    return strmap

def addKits(strs, dubSTRs, quadSTRs, ks, ids, idsToKeep, hgs, thekits, theids, thehgs, rejected, percentMissingSTRThreshold, allowable):
        
    for i in range(len(ids)):
        #if i % 20 == 0:
            #print(i, "/", len(ids))
        
        thehg = hgs[i].strip(" ")
        thisId = ids[i]

        if thehg == "?":
            rejected[(thisId)] = "unknown haplogroup"
        else:
            if str(allowable[i]) != "NaN":
                strmap = convertSTRrowToMap(ks, i, strs, dubSTRs, quadSTRs)
                modelInputFormattedSTRs = getValuesForPredictionFromAlleleArray(strmap, strs, dubSTRs, quadSTRs, percentMissingSTRThreshold)
                if modelInputFormattedSTRs != None:
                
                    thekits.append(np.array(modelInputFormattedSTRs))#[0:thelength])
                    theids.append(str(thisId))
                    thehgs.append(thehg)
    
                else:
                    rejected[str(thisId)] = strmap
        
            
         
import pandas as pd
def parseTrainCSV(strs, dubSTRs, quadSTRs, fil, modesIncluded, thekits, theids, thehgs, rejected, percentMissingSTRThreshold):
    k = pd.read_csv(fil)    
    ids = k["Kit Number"]
    hgs = k["Haplogroup"]
    allowable = k["Allowable Downstream"]
    
    addKits(strs, dubSTRs, quadSTRs, k, ids, [], hgs, thekits, theids, thehgs, rejected, percentMissingSTRThreshold, allowable)

def getClosestCutoff(uncertainidx, cutoff, thekits, thehgs, uncertainKits, uncertainAllowable):
    mindist = 100
    minclass = 0
    for i in range(len(thekits)):
        dist = 0
        for stridx in range(len(thekits[0])):
            dist += abs(thekits[i][stridx] - uncertainKits[uncertainidx][stridx])
        if dist < mindist and dist <= cutoff and thehgs[i] in uncertainAllowable[uncertainidx]:
            minclass = thehgs[i]
            mindist = dist
    if minclass != 0:
        return minclass
    else:
        return None
    
def getRefined(cutoff, thekits, theids, thehgs, uncertainKits, uncertainIds, uncertainAllowable):
    refinedPanelToSubpanel = 0
    refinedUnknownToPanel = 0
    for i in range(len(uncertainIds)):
        refined = getClosestCutoff(i,cutoff, thekits, thehgs, uncertainKits, uncertainAllowable)
        print('closest for', uncertainIds[i], refined)
        if refined is not None:
            if uncertainIds[i] in theids:
                thehgs[theids.index(uncertainIds[i])] = refined
                refinedPanelToSubpanel += 1
            else:
                theids.append(uncertainIds[i])
                thekits.append(uncertainKits[i])
                thehgs.append(refined)
                refinedUnknownToPanel += 1
    print('refined',refinedPanelToSubpanel,'panel to subpanel based on STR cutoff',cutoff)
    print('refined',refinedUnknownToPanel,'unknown to panel based on STR cutoff',cutoff)
    
from sklearn.ensemble import RandomForestClassifier



def getTrainingSamplesFromFile(fil, modesIncluded, percentMissingSTRThreshold):
    thekits = []
    theids = []
    thehgs = []
    
    rejected = {}
    (strs, dubSTRs, quadSTRs) = getSTRLabelsFromSets(modesIncluded)

    parseTrainCSV(strs, dubSTRs, quadSTRs, fil, modesIncluded, thekits, theids, thehgs, rejected, percentMissingSTRThreshold)
    return (thekits, theids, thehgs, rejected)

import random

def createTrainTest(x, y, ids):
    idsToHoldout = []
    hgMap = {}
    xH = []
    yH = []
    xT = []
    yT = []
    for idx in range(len(y)):
        if y[idx] not in hgMap:
            hgMap[y[idx]] = [ids[idx]]
        else:
            hgMap[y[idx]].append(ids[idx])
    
    for hg in hgMap:
        hgElems = len(hgMap[hg])
        if hgElems > 2 and hgElems < 10:
            idsToHoldout.append(hgMap[hg][random.randint(0, hgElems - 1)])
        else:
            for i in hgMap[hg]:
                if random.random() < 0.1:
                    idsToHoldout.append(i)
                    
    correctOrderIdsHoldout = []
    for theid in ids:
        idx = ids.index(theid)
        if theid in idsToHoldout:
            xH.append(x[idx])
            yH.append(y[idx])
            correctOrderIdsHoldout.append(theid)
        else:
            xT.append(x[idx])
            yT.append(y[idx])
    return (xT, yT, xH, yH, correctOrderIdsHoldout)

import time

def fromPredProbaGetPolicyValues(preds, predProba):
    policyValues = []
    for i in range(len(preds)):
        a = np.sort(predProba[i])        
        policyValues.append([preds[i],a[-1], a[-1] / a[-2]])
    return policyValues

def getPanelHier(panelHierFile):
    panelHier = {}
    with open(panelHierFile, "r") as f:
        for line in f.readlines():
            linesplit = line.strip("\n").split(",")
            child = linesplit[0]
            parent = linesplit[1]
            panelHier[child] = parent
    return panelHier

def aIsUpstreamB(a, b, panelHier):
    if a == b:
        return True
    if b in panelHier:
        return aIsUpstreamB(a, panelHier[b], panelHier)
    return False

def getErrorTypesAndPercentCorrect(preds, truth, panelHier):
    overSpecificityError = 0
    underSpecificityError = 0
    flatOutWrong = 0
    correct = 0
    for i in range(len(preds)):
        p = preds[i]
        t = truth[i]
        if p != t:
            if aIsUpstreamB(p, t, panelHier):
                #print("underspecific", p, "should be", t)
                underSpecificityError += 1
            else:
                if aIsUpstreamB(t, p, panelHier):
                    overSpecificityError += 1
                    #print("overspecific", p, "should be", t)
                else:
                    flatOutWrong += 1
                    #print("completely wrong", p, "should be", t)
        else:
            correct += 1
    #print(underSpecificityError, overSpecificityError, flatOutWrong, correct)
    return ([underSpecificityError, overSpecificityError, flatOutWrong, correct], correct / len(preds))

def getRawPredConfidenceMap(preds, predProbas, truth, ids, x):
    policyValues = fromPredProbaGetPolicyValues(preds, predProbas)
    rawPredMap = []
    wrongIndexes = []
    for j in range(len(policyValues)):
        policyValue = policyValues[j]
        ratio = policyValue[2]
        pred = policyValue[0]
        
        if pred == truth[j]:
            rawPredMap.append([ratio, True])
        else:
            rawPredMap.append([ratio, False])
            wrongIndexes.append((ids[j],truth[j],pred, x[j]))
    return rawPredMap, wrongIndexes

def lookupConfidence(rawPredMap, ratio):
    total = 0
    totalCorrect = 0
    for rawPred in rawPredMap:
        theratio, correct = rawPred
        if theratio > ratio * 0.8 and theratio < ratio * 1.2:
            total += 1
            if correct:
                totalCorrect += 1
    if total == 0:
        return -1
    return totalCorrect / total

        
def refinePredictionsPerPolicy(preds, predProbas, panelHier, policyA, policyB):
    policyValues = fromPredProbaGetPolicyValues(preds, predProbas)
    refined = []
    for policyValue in policyValues:
        pred = policyValue[0]
        maxPred = policyValue[1]
        ratio = policyValue[2]
        if ratio < policyB:
            if pred in panelHier:
                refined.append(panelHier[pred])
                #print(pred, 'refined based on policy to', panelHier[pred])
            else:
                refined.append(pred)
        else:
            refined.append(pred)
    return refined

def optimizePolicyParameters(preds, predProbas, truth, panelHier, policyARange, policyBRange, errorWeights):
    bestA = policyARange[0]
    bestB = policyBRange[0]
    refined = refinePredictionsPerPolicy(preds, predProbas, panelHier, bestA, bestB)
    (errorTotals, bestPercentCorrect) = getErrorTypesAndPercentCorrect(refined, truth, panelHier)
    bestUtility = getUtility(errorTotals, errorWeights)
    experimentMap = []
    
    for a in policyARange:
        for b in policyBRange:
            refined = refinePredictionsPerPolicy(preds, predProbas, panelHier, a, b)
            (thisErrTots, thisPercentCorrect) = getErrorTypesAndPercentCorrect(refined, truth, panelHier)
            thisUtility = getUtility(thisErrTots, errorWeights)
            experimentMap.append([a,b,thisErrTots[0], thisErrTots[1], thisErrTots[2], thisErrTots[3], thisPercentCorrect, thisUtility])
            if thisUtility > bestUtility:
                bestA = a
                bestB = b
                bestUtility = thisUtility
                bestPercentCorrect = thisPercentCorrect
    return (bestA, bestB, bestUtility, bestPercentCorrect, experimentMap)

def getUtility(errors, errorWeights):
    return errors[0] * errorWeights[0] + errors[1] * errorWeights[1] + errors[2] * errorWeights[2] + errors[3] * errorWeights[3]

def persistPolicy(policyA, policyB, policyFile):
    with open(policyFile, "w") as w:
        w.write("a," + str(policyA) + "\n")
        w.write("b," + str(policyB) + "\n")
    w.close()
    
def persistRawPredConfidence(rawPredConfidence, rawPredConfidenceFile):
    with open(rawPredConfidenceFile, "w") as w:
        for rawPredConf in rawPredConfidence:
            w.write(str(rawPredConf[0]) + "," + str(int(rawPredConf[1])) + "\n")
    w.close()

def persistRawPredErrors(rawPredErrors, rawPredErrorFile):
    with open(rawPredErrorFile, "w") as w:
        w.write("id,actual,predicted,str_alleles\n")
        for rawPredError in rawPredErrors:
            w.write(str(rawPredError[0]) + "," + rawPredError[1] + "," + rawPredError[2] + "," + " ".join(list(map(str,rawPredError[3]))) + "\n")
    w.close()
def readPolicy(policyFile):
    policyA = None
    policyB = None
    with open(policyFile, "r") as f:
        for line in f.readlines():
            linesplit = line.strip("\n").split(",")
            policyVar = linesplit[0]
            policy = linesplit[1]
            if policyVar == "a":
                policyA = float(policy)
            if policyVar == "b":
                policyB = float(policy)
    return (policyA, policyB)

def readRawPredConfidence(rawPredConfidenceFile):
    rawPredConfidence = []
    with open(rawPredConfidenceFile, "r") as f:
        for line in f.readlines():
            linesplit = line.strip("\n").split(",")
            rawPredConfidence.append([float(linesplit[0]), bool(int(linesplit[1]))])
    return rawPredConfidence


def persistExperimentMap(experimentMapFileStem, modesIncluded, expMap):
    with open(experimentMapFileStem + "_" + modesIncluded + ".csv", "w") as w:
        w.write(",".join(["a","b","under specific errors", "over specific errors", "completely wrong errors", "correct", "percent correct", "utility"]) + "\n")
        for entry in expMap:
            w.write(",".join(str(e) for e in entry) + "\n")
    w.close()

import pickle

import multiprocessing
    
class Experiment:
    def parallelExperiment(self, infile, modesIncluded, panelHier, policyFileStem, modelPickleFileStem, utilityWeights, experimentMapFileStem, percentMissingSTRThreshold, rfEstimators, rfMaxDepth):            
        (thekits, theids, thehgs, rejected) = getTrainingSamplesFromFile(infile, modesIncluded, percentMissingSTRThreshold)
        #if cutoff:
        #    getRefined(cutoff, thekits, theids, thehgs, uncertainKits, uncertainIds, uncertainAllowable)
        (x, ids, y) = excludeQuestionMarks(thekits, theids, thehgs)
        #accuracySum = 0
        (xT, yT, xH, yH, idsToHoldout) = createTrainTest(x, y, ids)
        clf = RandomForestClassifier(n_estimators=rfEstimators, max_depth=rfMaxDepth,random_state=0,class_weight="balanced")

        clf.fit(xT, yT)
        #start = time.time()
        preds = clf.predict(x)
        predProbas = clf.predict_proba(x)
        pklfile = open(getModesPklFile(modelPickleFileStem, modesIncluded), 'wb')
        pickle.dump(clf, pklfile, protocol=2)
        
        #underSpecificityErrorWeight = -1
        #overSpecificityErrorWeight = -3
        #completelyWrongErrorWeight = -5
        #rightWeight = 0
        if utilityWeights is None:
            errorWeights = [0, 0, 0, 0]
        else:
            errorWeights = utilityWeights
        
        #policyValues = fromPredProbaGetPolicyValues(preds, predProbas)
        #print(policyValues)
        
        
        #print(getUtility(errorTotals, errorWeights))
        optimum = optimizePolicyParameters(preds, predProbas, y, panelHier, [1],
                                  [1.01,1.02,1.05,1.08,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4], errorWeights)
        #print(optimum)
        #print(len(preds))
        def sortByPercentCorrect(t):
            return t[6]
        def sortByCompletelyWrong(t):
            return t[4]
        def sortByOverSpecific(t):
            return t[3]
        def sortByUnderSpecific(t):
            return t[2]
        def sortByUtility(t):
            return t[7]
        expMap = optimum[4]
        persistExperimentMap(experimentMapFileStem, modesIncluded, expMap)

        expMap.sort(key=sortByPercentCorrect)
        print(modesIncluded)
        print('best percent correct', expMap[-1])
        if utilityWeights == None:
            print('persisting best correct')
            best = expMap[-1]
            policyA = expMap[-1][0]
            policyB = expMap[-1][1]
        
        expMap.sort(key=sortByCompletelyWrong)
        print('least completely wrong', expMap[0])
        expMap.sort(key=sortByOverSpecific)
        print('least over specific', expMap[0])
        expMap.sort(key=sortByUnderSpecific)
        print('least under specific', expMap[0])
        expMap.sort(key=sortByUtility)
        print('highest utility', expMap[-1])
        if utilityWeights != None:
            print('persisting highest utility')
            policyA = expMap[-1][0]
            policyB = expMap[-1][1]
            best = expMap[-1]

        persistPolicy(policyA, policyB, getPolicyFile(policyFileStem, modesIncluded))
        
        (rawPredConf, errors) = getRawPredConfidenceMap(preds, predProbas, y, ids, x)
        persistRawPredConfidence(rawPredConf, getRawPredConfFile(policyFileStem, modesIncluded))
        
        persistRawPredErrors(errors, getRawPredErrorFile(policyFileStem, modesIncluded))
        print("Raw pred conf of 2.0", str(lookupConfidence(rawPredConf, 2.0)))
        createSpecificModelMetadata(modesIncluded, xT, xH, yT, modelPickleFileStem + "_" + modesIncluded + "_metadata", best)


def experimentErrorPolicy(infile, outfile, panelHierFile, policyFileStem, modelPickleFileStem, utilityWeights, experimentMapFileStem, percentMissingSTRThreshold, rfEstimators, rfMaxDepth):
    panelHier = getPanelHier(panelHierFile)
        
    processes = []
    for modesIncluded in modeCombos:
        e = Experiment()
        p = multiprocessing.Process(target=e.parallelExperiment, args=(infile, modesIncluded, panelHier, policyFileStem, modelPickleFileStem, utilityWeights, experimentMapFileStem, percentMissingSTRThreshold, rfEstimators, rfMaxDepth))
        p.start()
        processes.append(p)
    for p in processes:
        p.join()
    k = pd.read_csv(infile)
    samples = len(k["Haplogroup"])
    hgs = set(k["Haplogroup"])
    createGeneralModelMetadata(samples, hgs, rfEstimators, rfMaxDepth, modelPickleFileStem + "_general_metadata")
    
    
def getAlleleArrayFromFile(fil):
    with open(fil, "r") as f:
        filestrs = f.readline().strip("/n").split(",")
        return filestrs

def validateSTRQuery(strarray):
    for thestr in strarray:
        thesplit = thestr.split("=")
        if len(thesplit) == 1:
            return "STR query format error: '" + thestr + "' needs to be in format $ALLELE=$VALUE"
        for split in thesplit[1].split("-"):
            if thesplit[0] != "DYF399X":
                try:
                    float(split)
                except ValueError:
                    return "STR query format error: '" + thestr + "', allele value '" + split + "' not a float"
    return None
    
def getValuesForPredictionFromAlleleArray(strmap, strs, dubSTRs, quadSTRs, percentMissingSTRThreshold):
    thiskit = []
    missing = 0
    for STR in strs:
        if STR not in strmap or is_float(strmap[STR]) == False:
            #print('missing or not float', STR)
            missing += 1
            thiskit.append(float(0))
        else:    
            thiskit.append(float(strmap[STR]))
    for STR in dubSTRs:
        if STR not in strmap:
            #print('missing', STR)
            missing += 1
            thiskit.append(float(0))
            thiskit.append(float(0))
            thiskit.append(float(2))
        else:
            splits = strmap[STR].split("-")
            sl = len(splits)
            if sl < 2:
                for i in range(2-sl):
                    thiskit.append(float(0))
                for i in range(sl):
                    thiskit.append(float(splits[i]))
                thiskit.append(float(sl))
            else:
                thiskit.append(float(splits[0]))    
                thiskit.append(float(splits[sl-1]))
                thiskit.append(float(sl))
    for STR in quadSTRs:
        if STR not in strmap:
            #print('missing', STR)
            missing += 1
            thiskit.append(float(0))
            thiskit.append(float(0))
            thiskit.append(float(0))
            thiskit.append(float(0))
            thiskit.append(float(4))
        else:
            splits = strmap[STR].split("-")
            sl = len(splits)
            if sl < 4:
                for i in range(4-sl):
                    thiskit.append(float(0))
                for i in range(sl):
                    thiskit.append(float(splits[i]))
                thiskit.append(float(sl))
            else:
                thiskit.append(float(splits[0]))
                thiskit.append(float(splits[1]))
                thiskit.append(float(splits[sl-2]))
                thiskit.append(float(splits[sl-1]))
                thiskit.append(float(sl))
    percentMissing = float(missing) / (len(strs) + len(dubSTRs) + len(quadSTRs))
    if percentMissing > percentMissingSTRThreshold:
        #print("percent missing STRs:", percentMissing, "above threshold", percentMissingSTRThreshold)
        return None
    return thiskit
            

modeCombos = ["abcd","ftdna_111","abc","ftdna_67","ab","ftdna_37","cd","ftdna_25","a","ftdna_12","b","c","d"]

def getModesPklFile(stem, modes):
    return stem + "_" + modes + ".pkl"

def getPolicyFile(stem, modes):
    return stem + "_" + modes + ".csv"
def getRawPredConfFile(stem, modes):
    return stem + "_" + modes + "_prediction_confidence.csv"
def getRawPredErrorFile(stem, modes):
    return stem + "_" + modes + "_errors.csv"
def getSTRmap(queryAlleleArray):
    strmap = {}
    for thestr in queryAlleleArray:
        (strname, strval) = thestr.split("=")
        strmap[strname] = strval
    print(strmap)
    print(len(strmap.keys()))
    return strmap

import os

def exactlyMatchesAnyModeCombo(strmap, modeCombos):
    rejected = True
    modesIdx = 0
    predstrs = None
    while rejected and modesIdx < len(modeCombos):
        modesIncluded = modeCombos[modesIdx]
        (strs, dubSTRs, quadSTRs) = getSTRLabelsFromSets(modesIncluded)
        if len(strs + dubSTRs + quadSTRs) == len(strmap.keys()):
            predstrs = getValuesForPredictionFromAlleleArray(strmap, strs, dubSTRs, quadSTRs, 0)
            if predstrs != None:
                rejected = False
        modesIdx += 1
    return (predstrs, modesIncluded)
    
def predict(strAlleleString, panelHierarchy, policyFileStem, modelPickleFileStem, percentMissingSTRThreshold, haplogroupClassConfigPath, outputVersion = "normal"):
    queryAlleleArray = strAlleleString.split(",")
    validationMessage = validateSTRQuery(queryAlleleArray)
    if validationMessage != None:
        print(validationMessage)
        return validationMessage
    else:
        strmap = getSTRmap(queryAlleleArray)
        (predstrs, modesIncluded) = exactlyMatchesAnyModeCombo(strmap, modeCombos)
        if predstrs != None:
            return loadModelAndPredict(predstrs, panelHierarchy, modesIncluded, policyFileStem, modelPickleFileStem, haplogroupClassConfigPath, outputVersion)
        else:
            rejected = True
            modesIdx = 0
            while rejected and modesIdx < len(modeCombos):
                modesIncluded = modeCombos[modesIdx]
                (strs, dubSTRs, quadSTRs) = getSTRLabelsFromSets(modesIncluded)
                print(modesIncluded)
                predstrs = getValuesForPredictionFromAlleleArray(strmap, strs, dubSTRs, quadSTRs, percentMissingSTRThreshold)
                if predstrs != None:
                    rejected = False
                modesIdx += 1
            print(modesIncluded, strs, dubSTRs, quadSTRs)
            return loadModelAndPredict(predstrs, panelHierarchy, modesIncluded, policyFileStem, modelPickleFileStem, haplogroupClassConfigPath, outputVersion)

import json

def getPredictedHTML(key, haplogroupClassConfigPath):
    thejson = json.load(open(haplogroupClassConfigPath))
#    return thejson[key]["html"]
    return thejson[key]["html"] + '<br><br><img src="migration.jpg" height="50" width="50">&nbsp;' + '<a href="' + thejson[key]["migration"] + '" target="_blank">View Migration</a>&nbsp;Computed by PhyloGeographer from YFull and ancient samples'
    
def loadModelAndPredict(predstrs, panelHierarchy, modesIncluded, policyFileStem, modelPickleFileStem, haplogroupClassConfigPath, outputVersion = "normal"):
    if predstrs == None:
        print("Rejected, not enough STRs in sample to predict")
        return "Rejected, not enough STRs in sample to predict"
    else:
        print(modesIncluded, predstrs)
#        if modesIncluded = "a":
#            cutoff = 3
#        else:
#            cutoff = (len(modesIncluded) - 1) * 12
        #(thekits, theids, thehgs, uncertainKits, uncertainIds, uncertainAllowable, rejected) = getTrainingSamplesFromFile(infile, modesIncluded)
        #if cutoff:
            #getRefined(cutoff, thekits, theids, thehgs, uncertainKits, uncertainIds, uncertainAllowable)
        #(x, ids, y) = excludeQuestionMarks(thekits, theids, thehgs)
        modelPickleFile = getModesPklFile(modelPickleFileStem, modesIncluded)
        start = time.time()
        clf = pickle.load(open(modelPickleFile, 'rb'))
        end = time.time()
        
        
        preds = clf.predict([predstrs])
        predproba = clf.predict_proba([predstrs])
        
        
        print("predicted as", preds[0])
        predprobaclass = []
        print(clf.classes_)
        
        (policyA, policyB) = readPolicy(getPolicyFile(policyFileStem, modesIncluded))
        refined = refinePredictionsPerPolicy(preds, predproba, panelHierarchy, policyA, policyB)
        print('refined as', refined[0])
        
        def sortPredProba(t):
            return t[0]
        
        
        
        for i in range(len(predproba[0])):
            predprobaclass.append((predproba[0][i], clf.classes_[i]))
        predprobaclass.sort(key=sortPredProba)
        predprobaclass.reverse()
        print(predprobaclass)
        if outputVersion == "mini":
            print(createHTML(refined, predprobaclass, readRawPredConfidence(getRawPredConfFile(policyFileStem, modesIncluded))) )
        else:
            print(getPredictedHTML(refined[0], haplogroupClassConfigPath) + "<br><br>" + createHTML(refined, predprobaclass, readRawPredConfidence(getRawPredConfFile(policyFileStem, modesIncluded))) + "<b>Model Information</b>" + getSpecificModelMetadata(modelPickleFileStem, modesIncluded))
        return refined[0]

def createHTML(refined, predprobaclass, rawPredConfidence):
    thehtml = "<table border=\"1\"><tr><td>Haplogroup</td><td>Score</td><td>Score Bar</td></tr>"
    added = False
    for predproba in predprobaclass:
        thehtml += "<tr><td>" + predproba[1] + "</td><td>" + str(round(predproba[0],3)) + "</td><td>"
        if not added:
            accuracy = lookupConfidence(rawPredConfidence, predprobaclass[0][0] / predprobaclass[1][0])
            color = "red"
            if accuracy > 0.5:
                color = "yellow"
            if accuracy > 0.9:
                color = "green"
            if accuracy != -1:
                thehtml += "<font color=\""+color+"\">"
            else:
                thehtml += "<font color=\"red\">"
        for i in range(int(predproba[0] * 100)):
            thehtml += "&#9608;"
        if not added:
            thehtml += "</font>"
            if accuracy != -1:
                thehtml += "&nbsp;" + str(round(100 * accuracy,1))+ "% Confidence"
            else:
                thehtml += "&nbsp;No Confidence Information For This Prediction"
            added = True
            
        thehtml += "</td></tr>"
    thehtml += "</table>"
    return thehtml
    
    
def excludeQuestionMarks(thekits, theids, thehgs):
    removedNonesSTRs = []
    removedNonesIDs = []
    signals = []

    for i in range(len(theids)):
        if thehgs[i] != "?":
            removedNonesSTRs.append(thekits[i])
            removedNonesIDs.append(str(theids[i]))
            signals.append(thehgs[i])
    return (removedNonesSTRs, removedNonesIDs, signals)
    
    
def score(preds, truth, ids):
    classAccuracy = {}
    for sig in set(truth):
        classAccuracy[sig] = 0
    right = 0
    truthMatrix = {}
    wronglyPredictedIdsAs = {}
    for i in range(len(preds)):
        if preds[i] == truth[i]:
            right += 1
            classAccuracy[truth[i]] += 1
        else:
            if truth[i] not in truthMatrix:
                truthMatrix[truth[i]] = {}
            if preds[i] in truthMatrix[truth[i]]:
                truthMatrix[truth[i]][preds[i]] += 1
            else:
                truthMatrix[truth[i]][preds[i]] = 1
            if truth[i] not in wronglyPredictedIdsAs:
                wronglyPredictedIdsAs[truth[i]] = {}
            if preds[i] in wronglyPredictedIdsAs[truth[i]]:
                wronglyPredictedIdsAs[truth[i]][preds[i]].add(ids[i])
            else:
                wronglyPredictedIdsAs[truth[i]][preds[i]] = set([ids[i]])
    return (right / len(preds), classAccuracy, truthMatrix, wronglyPredictedIdsAs)

def writeResults(fil, trainTime, signals, testClasses, removedNonesSTRs, accuracy, classAccuracy, truthMatrix, wronglyPredictedIdsAs):
    iterations = 1
    with open(fil, "w") as w:
        w.write('model train time,' + str(round((trainTime)/60,2)) + " minutes\n")
        w.write('total classes in training,' + str(len(set(signals))) + "\n")
        w.write('classes in training,' + str(set(signals)) + "\n")
        w.write('total classes in test,' + str(len(set(testClasses))) + "\n")
        w.write('classes in test,' + str(set(testClasses)) + "\n")
        w.write('total samples,' + str(len(removedNonesSTRs)) + "\n")
        w.write('accuracy,'+ str(accuracy)+ "\n")
        def sortClassAccuracy(t):
            return t[1]
        classAcc = []
        for classAccuracyKey in classAccuracy:
            if classAccuracyKey in testClasses:
                classAcc.append([classAccuracyKey, classAccuracy[classAccuracyKey]])
        classAcc.sort(key=sortClassAccuracy)
        for i in range(len(classAcc)):
            w.write(classAcc[i][0] + ' accuracy, ' + str(round(classAcc[i][1]/iterations,3))+'\n')
        for i in range(len(classAcc)):
            if classAcc[i][0] in truthMatrix:
                truthclass = []
                for wrongpredkey in truthMatrix[classAcc[i][0]]:
                    truthclass.append([wrongpredkey, truthMatrix[classAcc[i][0]][wrongpredkey], wronglyPredictedIdsAs[classAcc[i][0]][wrongpredkey]])
                truthclass.sort(key=sortClassAccuracy)
                for j in range(len(truthclass)):
                    w.write(classAcc[i][0] + " wrongly predicted as " + truthclass[j][0] + ", " + str(round(truthclass[j][1] / iterations,3)) + "," + " ".join(truthclass[j][2]) + "\n")

def createGeneralModelMetadata(samples, classes, rfEstimators, rfMaxDepth, generalModelMetadataFile):
    if "?" in classes:
        classes.remove("?")
    classes = list(classes)
    classes.sort()
    thehtml = "Training Data:<br><br>"
    thehtml += "<table border=\"1\"><tr><td>Samples</td><td>" + str(samples) + "</td></tr><tr><td>Total Haplogroups</td><td>" + str(len(classes)) + "</td></tr><tr><td>Haplogroups</td><td>" + ", ".join(classes) + "</td></tr></table><br>"
    thehtml += "Random Forest Model Parameters:<br><br>"
    thehtml += "<table border=\"1\"><tr><td>Estimators</td><td>" + str(rfEstimators) + "</td></tr><tr><td>Max Depth</td><td>" + str(rfMaxDepth) + "</td></tr></table>"
    #stub for neural net    
    with open(generalModelMetadataFile, "w") as w:
        w.write(thehtml)
    w.close()

def createSpecificModelMetadata(modesIncluded, train, test, classesTrainedOn, specificModelMetadataFile, highestUtility):
    classesTrainedOn = list(set(classesTrainedOn))
    classesTrainedOn.sort()
    (strs, dubs, quads) = getSTRLabelsFromSets(modesIncluded)
    strLabels = list(set(strs + dubs + quads))
    strLabels.sort()
    errType1 = highestUtility[2]
    errType2 = highestUtility[3]
    errType3 = highestUtility[4]
    correct = highestUtility[5]
    totes = errType1 + errType2 + errType3 + correct
    underSpecificRateString = str(round(errType1 / totes * 100,1)) + "%"
    overSpecificRateString = str(round(errType2 / totes * 100,1)) + "%"
    wrongRateString = str(round(errType3 / totes * 100,1)) + "%"
    correctRateString = str(round(correct / totes * 100,1)) + "%"
    accurateOrUnderspecificRateString = str(round((errType1 + correct)/ totes * 100,1)) + "%"
    thehtml = "<table border=\"1\"><tr><td>Total STRs trained on</td><td>" + str(len(strLabels)) + "</td></tr><tr><td>STRs trained on</td><td>" + ", ".join(strLabels) + "</td></tr><tr><td>Training Samples</td><td>" + str(len(train)) + "</td></tr><tr><td>Test Samples</td><td>" + str(len(test)) + "</td></tr><tr><td>Total Haplogroup Classes Trained</td><td>" + str(len(classesTrainedOn)) + "</td></tr><tr><td>Haplogroup Classes Trained</td><td>" + ", ".join(classesTrainedOn) + "</td></tr><tr><td>Underspecificity Error</td><td>"+ underSpecificRateString +"</td></tr><tr><td>Overspecificity Error</td><td>"+overSpecificRateString+"</td></tr><tr><td>Other Error</td><td>"+wrongRateString+"</td></tr><tr><td>Accuracy</td><td>" +correctRateString+ "</td></tr><tr><td>Accurate or Underspecific</td><td>" +accurateOrUnderspecificRateString+ "</td></tr></table>"
    with open(specificModelMetadataFile, "w") as w:
        w.write(thehtml)
    w.close()
    
def getSpecificModelMetadata(modelPickleStem, modesIncluded):
    with open(modelPickleStem + "_" + modesIncluded + "_metadata", "r") as r:
        return r.readline()
    r.close()
