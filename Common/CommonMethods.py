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
    
def addKits(strs, dubSTRs, quadSTRs, ks, ids, idsToKeep, hgs, allowable, thekits, theids, thehgs, uncertainKits, uncertainIds, uncertainAllowable, rejected):
        
    for i in range(len(ids)):
        if i % 20 == 0:
            print(i, "/", len(ids))
        thehg = hgs[i].strip(" ")
        if str(ks[strs[0]][i])  != "nan":
            thiskit = []
            ignore = False
            thisId = str(ids[i]).strip(" ")
            if len(idsToKeep) > 0 and thisId not in idsToKeep:
                ignore = True
            for STR in strs:
                if str(ks[STR][i]) == "nan" or is_float(ks[STR][i]) == False:
                    print("ignored", ks[STR][i])
                    ignore = True
                else:    
                    thiskit.append(float(ks[STR][i]))
            for STR in dubSTRs:
                splits = str(ks[STR][i]).split("-")
                sl = len(splits)
                if sl < 2:
                    print("ignored", ks[STR][i])
                    ignore = True
                else:
                    thiskit.append(float(splits[0]))    
                    thiskit.append(float(splits[sl-1]))
                    thiskit.append(float(sl))
            for STR in quadSTRs:
                splits = str(ks[STR][i]).split("-")
                sl = len(splits)
                if sl < 4:
                    ignore = True
                else:
                    thiskit.append(float(splits[0]))
                    thiskit.append(float(splits[1]))
                    thiskit.append(float(splits[sl-2]))
                    thiskit.append(float(splits[sl-1]))
                    thiskit.append(float(sl))

            if not ignore:
                if thehg != "?":
                    thekits.append(np.array(thiskit))#[0:thelength])
                    theids.append(str(thisId))
                    thehgs.append(thehg)
                if str(allowable[i]) != "nan":
                    uncertainKits.append(np.array(thiskit))
                    uncertainIds.append(str(thisId))
                    uncertainAllowable.append(allowable[i].split(":"))
                #print('added', ids[i],thiskit)
            else:
                rejected[ids[i]] = thiskit
                #print('ignored', ids[i],thiskit)
        else:
            rejected[ids[i]] = "TOTAL REJECT " + str(ks[strs[0]][i])
         
import pandas as pd
def parseTrainCSV(strs, dubSTRs, quadSTRs, fil, modesIncluded, thekits, theids, thehgs, uncertainKits, uncertainIds, uncertainAllowables, rejected):
    k = pd.read_csv(fil)    
    ids = k["Kit Number"]
    hgs = k["Haplogroup"]
    allowable = k["Allowable Downstream"]
    
    addKits(strs, dubSTRs, quadSTRs, k, ids, [], hgs, allowable, thekits, theids, thehgs, uncertainKits, uncertainIds, uncertainAllowables, rejected)

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

clf = RandomForestClassifier(n_estimators=100, max_depth=2,random_state=0,class_weight="balanced")

def getTrainingSamplesFromFile(fil, modesIncluded):
    thekits = []
    theids = []
    thehgs = []
    
    uncertainKits = []
    uncertainIds = []
    uncertainAllowable = []
    rejected = {}
    (strs, dubSTRs, quadSTRs) = getSTRLabelsFromSets(modesIncluded)

    parseTrainCSV(strs, dubSTRs, quadSTRs, fil, modesIncluded, thekits, theids, thehgs, uncertainKits, uncertainIds, uncertainAllowable, rejected)
    return (thekits, theids, thehgs, uncertainKits, uncertainIds, uncertainAllowable, rejected)

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
    for theid in ids:
        idx = ids.index(theid)
        if theid in idsToHoldout:
            xH.append(x[idx])
            yH.append(y[idx])
        else:
            xT.append(x[idx])
            yT.append(y[idx])
    return (xT, yT, xH, yH, idsToHoldout)

import time

def experiment(infile, modesIncluded, cutoff, outfile):
    (thekits, theids, thehgs, uncertainKits, uncertainIds, uncertainAllowable, rejected) = getTrainingSamplesFromFile(infile, modesIncluded)
    if cutoff:
        getRefined(cutoff, thekits, theids, thehgs, uncertainKits, uncertainIds, uncertainAllowable)
    (x, ids, y) = excludeQuestionMarks(thekits, theids, thehgs)
    (xT, yT, xH, yH, idsToHoldout) = createTrainTest(x, y, ids)
    clf.fit(xT, yT)
    start = time.time()
    preds = clf.predict(xH)
    end = time.time()
    (accuracy, classAccuracy, truthMatrix, wronglyPredictedIdsAs) = score(preds, yH, idsToHoldout)
    writeResults(outfile, end - start, y, yH, x, accuracy, classAccuracy, truthMatrix, wronglyPredictedIdsAs)

    
def getValuesForPrediction(fil, strs, dubSTRs, quadSTRs):
    with open(fil, "r") as f:
        filestrs = f.readline().split(",")
        strmap = {}
        for thestr in filestrs:
            (strname, strval) = thestr.split("=")
            strmap[strname] = strval
        ignore = False
        thiskit = []
        for STR in strs:
            if STR not in strmap or is_float(strmap[STR]) == False:
                ignore = True
            else:    
                thiskit.append(float(strmap[STR]))
        for STR in dubSTRs:
            if STR not in strmap:
                ignore = True
            else:
                splits = strmap[STR].split("-")
                sl = len(splits)
                if sl < 2:
                    ignore = True
                else:
                    thiskit.append(float(splits[0]))    
                    thiskit.append(float(splits[sl-1]))
                    thiskit.append(float(sl))
        for STR in quadSTRs:
            if STR not in strmap:
                ignore = True
            else:
                splits = strmap[STR].split("-")
                sl = len(splits)
                if sl < 4:
                    ignore = True
                else:
                    thiskit.append(float(splits[0]))
                    thiskit.append(float(splits[1]))
                    thiskit.append(float(splits[sl-2]))
                    thiskit.append(float(splits[sl-1]))
                    thiskit.append(float(sl))
        if ignore:
            return None
        else:
            return thiskit
            
import os

def predict(infile, dataDir, sampleId):
    predfile = os.path.join(dataDir,sampleId,"str")
    outfile = os.path.join(dataDir,sampleId,"prediction")
    modeCombos = ["abcd","abc","ab","cd","a","b","c","d"]
    rejected = True
    modesIdx = 0
    while rejected and modesIdx < len(modeCombos):
        modesIncluded = modeCombos[modesIdx]
        (strs, dubSTRs, quadSTRs) = getSTRLabelsFromSets(modesIncluded)
        predstrs = getValuesForPrediction(predfile, strs, dubSTRs, quadSTRs)
        if predstrs != None:
            rejected = False
        modesIdx += 1
    print(strs, dubSTRs, quadSTRs)
    if rejected:
        print("Rejected, not enough STRs in sample to predict")
    else:
        print(predstrs)
#        if modesIncluded = "a":
#            cutoff = 3
#        else:
#            cutoff = (len(modesIncluded) - 1) * 12
        (thekits, theids, thehgs, uncertainKits, uncertainIds, uncertainAllowable, rejected) = getTrainingSamplesFromFile(infile, modesIncluded)
        #if cutoff:
            #getRefined(cutoff, thekits, theids, thehgs, uncertainKits, uncertainIds, uncertainAllowable)
        (x, ids, y) = excludeQuestionMarks(thekits, theids, thehgs)
        start = time.time()
        clf.fit(x, y)
        end = time.time()
        
        preds = clf.predict([predstrs])
        predproba = clf.predict_proba([predstrs])
        print("predicted as", preds[0])
        predprobaclass = []
        print(clf.classes_)
        
        def sortPredProba(t):
            return t[0]
        
        
        
        for i in range(len(predproba[0])):
            predprobaclass.append((predproba[0][i], clf.classes_[i]))
        predprobaclass.sort(key=sortPredProba)
        predprobaclass.reverse()
        print(predprobaclass)
        with open(outfile, "w") as w:
            w.write("STR sets used," + modesIncluded + "\n")
            w.write("model train minutes," + str(round((end - start)/60,3)) + "\n")
            for p in predprobaclass:
                w.write(p[1] + "," + str(round(p[0],4)) + "\n")
        w.close()
    
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
