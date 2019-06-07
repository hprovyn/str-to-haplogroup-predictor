# -*- coding: utf-8 -*-
"""
Created on Wed May 22 14:35:06 2019

@author: hunte
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 15:47:45 2018

@author: hunte
"""

import numpy as np

import pandas as pd

import time

start = time.time()
def parseKitsCSV(projectName, fil):    
    k = pd.read_csv(fil)
    
    ids = k["id"]
    lat = k["latitude"]
    lon = k["longitude"]
    clade = k["clade"]
    neg = k["negatives"]
    ybp = k["ybp"]
    thekits = []
    for i in range(len(ids)):
        thiskit = {"id": ids[i], "clade": clade[i], "lat": lat[i], "lon": lon[i] }
        if projectName != "J-M241":
            thiskit["id"] = "YF_" + thiskit["id"]
        n = neg[i]
        if n != "" and str(n) != 'nan':
            thiskit["neg"] = n.split(":")
        else:
            thiskit["neg"] = []
        y = ybp[i]
        if y!= "" and str(y) != 'nan':
            thiskit["ybp"] = y
        else:
            thiskit["ybp"] = 65
        thekits.append(thiskit)
    return thekits


def parseFTDNACSV(fil):
    k = pd.read_csv(fil)
    
    ids = k["Kit Number"]
    hgs = k["Haplogroup"]
    allowable = k["Allowable Downstream"]
    addKits(k, ids, [], hgs, allowable)

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

a_group = ["DYS391","DYS389I","DYS437","DYS439","DYS389II","DYS438","DYS426","DYS393","YCAII","DYS390","DYS385","Y-GATA-H4","DYS388","DYS447","DYS19","DYS392"]
b_group = ["DYS458","DYS455","DYS454","DYS464","DYS448","DYS449","DYS456","DYS576","CDY","DYS460","DYS459","DYS570","DYS607","DYS442"]
c_group = ["DYS728","DYS723","DYS711","DYR76","DYR33","DYS727","DYR157","DYS713","DYS531","DYS578","DYF395","DYS590","DYS537","DYS641","DYS472","DYF406S1","DYS511","DYS557","DYS490","DYS446","DYS481","DYS413","DYS534","DYS450","DYS425","DYS594","DYS444","DYS520","DYS436","DYS565","DYS572","DYS617","DYS568","DYS487","DYS640","DYS492"]
d_group = ["DYR112","DYS518","DYS614","DYS626","DYS644","DYS684","DYS710","DYS485","DYS632","DYS495","DYS540","DYS714","DYS716","DYS717","DYS505","DYS556","DYS549","DYS589","DYS522","DYS494","DYS533","DYS636","DYS575","DYS638","DYS462","DYS452","DYS445","Y-GATA-A10","DYS463","DYS441","Y-GGAAT-1B07","DYS525","DYS712","DYS593","DYS650","DYS532","DYS715","DYS504","DYS513","DYS561","DYS552","DYS726","DYS635","DYS587","DYS643","DYS497","DYS510","DYS434","DYS461","DYS435"]

import sys
if len(sys.argv) > 1:
    strFile1 = sys.argv[1]
    outFile = sys.argv[2]
    modesIncluded = sys.argv[3]
    cutoff = int(sys.argv[4])
    iterations = int(sys.argv[5])

strLabelsUsed = []

rejected = {}

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
    
def is_float(value):
    try:
        float(value)
        return True
    except:
        return False

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

thekits = []
theids = []
thehgs = []

uncertainKits = []
uncertainIds = []
uncertainAllowable = []
def addKits(ks, ids, idsToKeep, hgs, allowable):    
    
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


parseFTDNACSV(strFile1)

def getRefined(cutoff):
    refinedPanelToSubpanel = 0
    refinedUnknownToPanel = 0
    for i in range(len(uncertainIds)):
        refined = getClosestCutoff(i,cutoff)
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

def getClosestCutoff(uncertainidx, cutoff):
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

def writeRFPreds(updated):
    w = open("rfPreds","w+")
    for update in updated:
        w.write(update + "," + ",".join(updated[update]))
        w.write("\n")
    w.close()

from sklearn.ensemble import RandomForestClassifier
clf = RandomForestClassifier(n_estimators=100, max_depth=2,random_state=0,class_weight="balanced")

if cutoff:
    getRefined(cutoff)

removedNonesSTRs = []
removedNonesIDs = []
signals = []

for i in range(len(theids)):
    if thehgs[i] != "?":
        removedNonesSTRs.append(thekits[i])
        removedNonesIDs.append(str(theids[i]))
        signals.append(thehgs[i])

import random

mutationRates = [0.00076,0.00311,0.00151,0.00265,0.00226,0.00009,0.00022,0.00477,0.00186,0.00052,0.00056,0.00814,0.00132,0.00016,0.00016,0.00264,0.00099,0.00135,0.00838,0.00566,0.00402,0.00208,0.00123,0.00735,0.00411,0.01022,0.0079,0.03531,0.00324,0.00055,0.00037,0.00008,0.00031,0.00054,0.00057,0.00018,0.00001,0.00154,0.00128,0.00018,0.00202,0.00321,0.00029,0.00018,0.00019,0.00832,0.0002,0.00321,0.00544,0.00245,0.00365,0.00042,0.00053,0.00097,0.00212,0.00034,0.00042,0.00087,0.01722,0.00204,0.00001,0.00135,0.00077,0.00502,0.00042,0.00265,0.00355,0.00142,0.00432,0.00192,0.00103,0.00008,0.00362,0.00099,0.00018,0.00063,0.00046,0.00223,0.00082,0.00278,0.00095,0.00238,0.00081,0.00121,0.02194,0.00022,0.00662,0.00259,0.00364,0.00414,0.00275,0.00198,0.00253,0.0001,0.00373,0.00097,0.0018,0.00093,0.00298,0.00018,0.0021,0.00005]
mrateLookupSTRs = ["DYS393","DYS390","DYS19","DYS391","DYS385","DYS426","DYS388","DYS439","DYS389I","DYS392","DYS389II","DYS458","DYS459","DYS455","DYS454","DYS447","DYS437","DYS448","DYS449","DYS464","DYS460","Y-GATA-H4","YCAII","DYS456","DYS607","DYS576","DYS570","CDY","DYS442","DYS438","DYS531","DYS578","DYF395","DYS590","DYS537","DYS641","DYS472","DYF406S1","DYS511","DYS425","DYS413","DYS557","DYS594","DYS436","DYS490","DYS534","DYS450","DYS444","DYS481","DYS520","DYS446","DYS617","DYS568","DYS487","DYS572","DYS640","DYS492","DYS565","DYS710","DYS485","DYS632","DYS495","DYS540","DYS714","DYS716","DYS717","DYS505","DYS556","DYS549","DYS589","DYS522","DYS494","DYS533","DYS636","DYS575","DYS638","DYS462","DYS452","DYS445","Y-GATA-A10","DYS463","DYS441","Y-GGAAT-1B07","DYS525","DYS712","DYS593","DYS650","DYS532","DYS715","DYS504","DYS513","DYS561","DYS552","DYS726","DYS635","DYS587","DYS643","DYS497","DYS510","DYS434","DYS461","DYS435"]

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

totalAccuracy = 0

truthMatrix = {}
wronglyPredictedIdsAs = {}

classAccuracy = {}
for sig in set(signals):
    classAccuracy[sig] = 0
    
def score(preds, truth, ids):
    right = 0
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
    return right / len(preds)
testClasses = None


for i in range(iterations):
    
    (xT, yT, xH, yH, idsH) = createTrainTest(removedNonesSTRs, signals, removedNonesIDs)
    clf.fit(xT, yT)
    preds = clf.predict(xH)
    totalAccuracy += score(preds, yH, idsH)
    testClasses = set(yH)

print(totalAccuracy / iterations, "overall accuracy")

def writeResults(fil):
    end = time.time()

    with open(fil, "w") as w:
        w.write('elapsed time,' + str(round((end-start)/60,2)) + " minutes\n")
        w.write('trainClassesLen,' + str(len(set(signals))) + "\n")
        w.write('trainClasses,' + str(set(signals)) + "\n")
        w.write('testClassesLen,' + str(len(set(testClasses))) + "\n")
        w.write('testClasses,' + str(set(testClasses)) + "\n")
        w.write('samples,' + str(len(removedNonesSTRs)) + "\n")
        w.write('overall accuracy,'+ str(totalAccuracy / iterations)+ "\n")
        def sortClassAccuracy(t):
            return t[1]
        classAcc = []
        for classAccuracyKey in classAccuracy:
            if classAccuracyKey in yH:
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

writeResults(outFile)