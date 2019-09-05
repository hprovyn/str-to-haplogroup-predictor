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

    strLabelsUsed = []    
    strs = []
    dubSTRs = []
    quadSTRs = []
    
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
        strmap[thestr] = ks[thestr][sampleIdx]
    return strmap

def addKits(strs, dubSTRs, quadSTRs, ks, ids, idsToKeep, hgs, thekits, theids, thehgs, rejected, percentMissingSTRThreshold):
        
    for i in range(len(ids)):
        #if i % 20 == 0:
            #print(i, "/", len(ids))
        
        thehg = hgs[i].strip(" ")
        thisId = ids[i]

        if thehg == "?":
            rejected[(thisId)] = "unknown haplogroup"
        else:
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
    #allowable = k["Allowable Downstream"]
    
    addKits(strs, dubSTRs, quadSTRs, k, ids, [], hgs, thekits, theids, thehgs, rejected, percentMissingSTRThreshold)

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

clf = RandomForestClassifier(n_estimators=300, max_depth=7,random_state=0,class_weight="balanced")

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
        if hgElems > 2:
            idsToHoldout.append(hgMap[hg][random.randint(0, hgElems - 1)])
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

def experiment(infile, modesIncluded, cutoff, outfile, iterations):
    (thekits, theids, thehgs, uncertainKits, uncertainIds, uncertainAllowable, rejected) = getTrainingSamplesFromFile(infile, modesIncluded)
    if cutoff:
        getRefined(cutoff, thekits, theids, thehgs, uncertainKits, uncertainIds, uncertainAllowable)
    (x, ids, y) = excludeQuestionMarks(thekits, theids, thehgs)
    accuracySum = 0
    for i in range(iterations):
        (xT, yT, xH, yH, idsToHoldout) = createTrainTest(x, y, ids)
        clf.fit(xT, yT)
        start = time.time()
        preds = clf.predict(xH)
        end = time.time()
        (accuracy, classAccuracy, truthMatrix, wronglyPredictedIdsAs) = score(preds, yH, idsToHoldout)
        accuracySum += accuracy
    writeResults(outfile, end - start, y, yH, x, accuracySum / iterations, classAccuracy, truthMatrix, wronglyPredictedIdsAs)


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

def refinePredictionsPerPolicy(preds, predProbas, panelHier, policyA, policyB):
    policyValues = fromPredProbaGetPolicyValues(preds, predProbas)
    refined = []
    for policyValue in policyValues:
        pred = policyValue[0]
        maxPred = policyValue[1]
        ratio = policyValue[2]
        if maxPred < policyA or ratio < policyB:
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

def persistExperimentMap(experimentMapFileStem, modesIncluded, expMap):
    with open(experimentMapFileStem + "_" + modesIncluded + ".csv", "w") as w:
        w.write(",".join(["a","b","under specific errors", "over specific errors", "completely wrong errors", "correct", "percent correct", "utility"]) + "\n")
        for entry in expMap:
            w.write(",".join(str(e) for e in entry) + "\n")
    w.close()

import pickle

import multiprocessing
    
class Experiment:
    def parallelExperiment(self, infile, modesIncluded, panelHier, policyFileStem, modelPickleFileStem, utilityWeights, experimentMapFileStem, percentMissingSTRThreshold):            
        (thekits, theids, thehgs, rejected) = getTrainingSamplesFromFile(infile, modesIncluded, percentMissingSTRThreshold)
        #if cutoff:
        #    getRefined(cutoff, thekits, theids, thehgs, uncertainKits, uncertainIds, uncertainAllowable)
        (x, ids, y) = excludeQuestionMarks(thekits, theids, thehgs)
        #accuracySum = 0
        (xT, yT, xH, yH, idsToHoldout) = createTrainTest(x, y, ids)
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
        
        (errorTotals, percentCorrect) = getErrorTypesAndPercentCorrect(preds, y, panelHier)
        #print(getUtility(errorTotals, errorWeights))
        optimum = optimizePolicyParameters(preds, predProbas, y, panelHier, [.11,.12,.13,.14,.15,.16,.17,.18,.19,.2,.21,.22,.23,.24,.25,.26,.27,.28,.29,.3],
                                 [1.02, 1.05,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.5,3,3.5,4], errorWeights)
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
            policyA = expMap[-1][0]
            policyB = expMap[-1][1]
        
        persistPolicy(policyA, policyB, getPolicyFile(policyFileStem, modesIncluded))

def experimentErrorPolicy(infile, outfile, panelHierFile, policyFileStem, modelPickleFileStem, utilityWeights, experimentMapFileStem, percentMissingSTRThreshold):
    panelHier = getPanelHier(panelHierFile)
        
    processes = []
    for modesIncluded in modeCombos:
        e = Experiment()
        p = multiprocessing.Process(target=e.parallelExperiment, args=(infile, modesIncluded, panelHier, policyFileStem, modelPickleFileStem, utilityWeights, experimentMapFileStem, percentMissingSTRThreshold))
        p.start()        
        processes.append(p)
    for p in processes:
        p.join()
    
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
            print('missing or not float', STR)
            missing += 1
            thiskit.append(float(0))
        else:    
            thiskit.append(float(strmap[STR]))
    for STR in dubSTRs:
        if STR not in strmap:
            print('missing', STR)
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
            print('missing', STR)
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
        print("percent missing STRs:", percentMissing, "above threshold", percentMissingSTRThreshold)
        return None
    return thiskit
            

modeCombos = ["abcd","abc","ab","cd","a","b","c","d"]

def getModesPklFile(stem, modes):
    return stem + "_" + modes + ".pkl"

def getPolicyFile(stem, modes):
    return stem + "_" + modes + ".csv"

def getSTRmap(queryAlleleArray):
    strmap = {}
    for thestr in queryAlleleArray:
        (strname, strval) = thestr.split("=")
        strmap[strname] = strval
    print(strmap)
    print(len(strmap.keys()))
    return strmap

import os
    
def predict(strAlleleString, panelHierarchy, policyFileStem, modelPickleFileStem, percentMissingSTRThreshold):
    queryAlleleArray = strAlleleString.split(",")
    validationMessage = validateSTRQuery(queryAlleleArray)
    if validationMessage != None:
        print(validationMessage)
        return validationMessage
    else:
        strmap = getSTRmap(queryAlleleArray)
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
        return loadModelAndPredict(predstrs, panelHierarchy, modesIncluded, policyFileStem, modelPickleFileStem)
    
def loadModelAndPredict(predstrs, panelHierarchy, modesIncluded, policyFileStem, modelPickleFileStem):
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
# =============================================================================
#         with open(outfile, "w") as w:
#             w.write("STR sets used," + modesIncluded + "\n")
#             w.write("model load seconds," + str(round((end - start),3)) + "\n")
#             for p in predprobaclass:
#                 w.write(p[1] + "," + str(round(p[0],4)) + "\n")
#         w.close()
# =============================================================================
        print(createHTML(refined, predprobaclass))
        return refined[0]

def createHTML(refined, predprobaclass):
    thehtml = "Predicted: " + refined[0] + "<br><br>"
    thehtml += "<table border=\"1\"><tr><td>Haplogroup</td><td>Score</td></tr>"
    for predproba in predprobaclass:
        thehtml += "<tr><td>" + predproba[1] + "</td><td>" + str(round(predproba[0],3)) + "</td></tr>"
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
    
def train(x, y):
    clf.fit(x, y)
    
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
