# -*- coding: utf-8 -*-
"""
Created on Fri May 24 00:23:13 2019

@author: hunte
"""

import json
import tabix

hierarchy = {}
snps = {}

def parseTreeJSON(fil):
    root = json.load(open(fil))
    recurseTreeJson(root, hierarchy, snps)
    return (root["id"], hierarchy, snps)

def parseSNPsString(snpsString):
    thesnps = []
    for snps in snpsString.split(", "):
        for snp in snps.split("/"):
            thesnps.append(snp)
    return thesnps
            
def recurseTreeJson(node, hierarchy, snps):
    if "children" in node:
        for child in node["children"]:
            hierarchy[child["id"]] = node["id"]
            snps[child["id"]] = parseSNPsString(child["snps"])
            recurseTreeJson(child, hierarchy, snps)

def createCSVforRF(outfile, headers, kits):
    with open(outfile, 'w') as f:
        f.write(",".join(["Kit Number", "Haplogroup", "Allowable Downstream"] + headers) + "\n")
        for kit in kits:
            strArray = []
            for h in headers:
                if h in kits[kit]["STRs"]:
                    strArray.append(kits[kit]["STRs"][h])
                else:
                    strArray.append("0")
            f.write(",".join([kit, kits[kit]["Haplogroup"], ":".join(kits[kit]["allowableDownstream"])] + strArray) + "\n")
        f.close()

class ParallelGetKit():
    def getKit(self, tabixFilePath, theid, specificHg, panelMap, panelHier, queue):
        kitToAdd = {}
        kitToAdd["Haplogroup"] = specificHg
        kitToAdd["STRs"] = {}
        tb = tabix.open(tabixFilePath)

        strResults = tb.querys("str:" + theid + "-" + theid)
        for strResult in strResults:
            thestr = strResult[3]
            allele = strResult[4].replace("\\","")         
            kitToAdd["STRs"][thestr] = allele
        if len(kitToAdd["STRs"]) >= 15:
            kitToAdd["allowableDownstream"] = getAllowableDownstream([], kitToAdd["Haplogroup"], panelMap, tb, theid, panelHier)

            queue.put((theid, kitToAdd))
        else:
            queue.put((0,0))

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

import multiprocessing

import time
            
def parseKits(haplogroupFile, tabixFilePath, hierarchy, panelMap, panelHier, maxThreads):
    kits = {}
    start = time.time()
    theids = []
    thehgs = []
    with open(haplogroupFile, 'r') as f:
        for line in f.readlines():
            (theid, clade) = line.split(",")
            thehg = clade.replace("\n","")
            if thehg != "?":
                theids.append(theid)
                thehgs.append(thehg)
    f.close()
    mostSpecificPanels = {}
    for hg in set(thehgs):
        if hg != "?":
            panels = getPanels(hg, panelMap)
            if len(panels) > 0:
                #print(panels)
                panel = getMostSpecificPanel(panels, panelMap)
                print('most specific', panel, 'for', hg)                    
                mostSpecificPanels[hg] = panel
            else:
                mostSpecificPanels[hg] = "?"
        else:
            mostSpecificPanels[hg] = "?"
    
    idChunks = list(chunks(theids, maxThreads))
    cladeChunks = list(chunks(thehgs, maxThreads))
    for i in range(len(idChunks)):
        idChunk = idChunks[i]
        cladeChunk = cladeChunks[i]
        print (idChunk, cladeChunk)
        processes = []
        queue = multiprocessing.Queue()
        for j in range(len(idChunk)):
            pgk = ParallelGetKit()
            p = multiprocessing.Process(target=pgk.getKit, args=(tabixFilePath, idChunk[j], mostSpecificPanels[cladeChunk[j]], panelMap, panelHier, queue))
            p.start()
            processes.append(p)
        rets = []
        for p in processes:
            ret = queue.get()
            rets.append(ret)
        for p in processes:
            p.join()
        for ret in rets:
            if ret[0] != 0:
            	kits[ret[0]] = ret[1]
    
    end = time.time()
    print("elapsed time processing specific and allowable downstream SNPs for " + str(len(theids)) + " samples " + str(round((end - start) / 60,2)) + " min")
    print(str(round((end - start) / len(theids),2)) + " sec per sample")

    return kits

def parseKitsMulti(haplogroupFile, tb, hierarchy, panelMap):
    kits = {}
    theids = []
    theclades = []
        
    with open(haplogroupFile, 'r') as f:
        for line in f.readlines():
            (theid, clade) = line.split(",")
            theids.append(theid)
            theclades.append(clade)
        f.close()
    idChunks = list(chunks(theids, 5))
    cladeChunks = list(chunks(theclades, 5))
    for i in range(len(idChunks)):
        idChunk = idChunks[i]
        cladeChunk = cladeChunks[i]
        print (idChunk, cladeChunk)
        processes = []
        queue = multiprocessing.Queue()
        for j in range(len(idChunk)):
            pgk = ParallelGetKit()
                    
            p = multiprocessing.Process(target=pgk.getKit, args=(tb, idChunk[j], cladeChunk[j], panelMap, queue))
            p.start()
            processes.append(p)
        rets = []
        for p in processes:
            ret = queue.get()
            rets.append(ret)
        for p in processes:
            p.join()
        for ret in rets:
            if ret[0] != 0:
            	kits[ret[0]] = ret[1]
    print(str(len(kits)) + " kits converted to RF input format via multiprocessing")
    return kits

def parseHeaders(kits):
    headers = set([])
    for kit in kits:
        for thestr in kits[kit]["STRs"]:
            headers.add(thestr)
    return list(headers)

def getUpstream(branch):
    thebranch = branch
    seq = [thebranch]
    while thebranch in hierarchy:
        thebranch = hierarchy[thebranch]
        seq.append(thebranch)
    if '' in seq:
        seq.remove('')
    return seq

def getUpstreamStop(branch, stop):
    thebranch = branch
    seq = [thebranch]
    while thebranch in hierarchy and thebranch != stop:
        thebranch = hierarchy[thebranch]
        seq.append(thebranch)
    if '' in seq:
        seq.remove('')
    return seq

def writePanelTree(panelMap, panelHierFile):
    justBranches = []
    branchToPanelMap = {}
    panelHier = {}
    for p in panelMap:
        thebranch = panelMap[p][0][0]
        justBranches.append(thebranch)
        branchToPanelMap[thebranch] = p
    for branch in branchToPanelMap:
        added = False
        ups = getUpstream(branch)
        ups.remove(branch)
        for up in ups:
            if not added:
                if up in branchToPanelMap:
                    panelHier[branchToPanelMap[branch]] = branchToPanelMap[up]
                    added = True
    with open(panelHierFile, "w") as w:
        for p in panelHier:
            w.write(p + "," + panelHier[p] + "\n")
    w.close()            
                
def getAllowableDownstream(negs, snpPredictedClade, panelMap, tb, theid, panelHier):
    allowableDownstream = list(panelMap.keys())
    if snpPredictedClade != "?":
        allowableDownstream = []
        for panel in panelMap:
            if panel != snpPredictedClade:
                if CommonMethods.aIsUpstreamB(snpPredictedClade, panel, panelHier):
                    allowableDownstream.append(panel)
    if len(allowableDownstream) > 0:
        toremove = []
        negs = []
        negResults = tb.querys("neg:" + theid + "-" + theid)
        for negResult in negResults:
            negs.append(negResult[3])
        for allowable in allowableDownstream:
            negated = False
            for upstreambranch in panelMap[allowable][0]:
                upstop = 'null'
                if snpPredictedClade != "?":
                    upstop = panelMap[snpPredictedClade][0][0]
                for upstream in getUpstreamStop(upstreambranch, upstop):
                    if not negated and any(n in snps[upstream] for n in negs):
                        #print(allowable, 'disallowed due to a negative', negs, 'in upstream',upstream,snps[upstream])
                        negated = True
            if negated:
                toremove.append(allowable)
        for torem in toremove:
            allowableDownstream.remove(torem)
    return allowableDownstream
        
    
def getPanels(branch, panelMap):
    panels = []
    branchReplacedAsterisk = branch.replace("*","")
    for panel in panelMap:
        positive = False
        for panelBranch in panelMap[panel][0]:
            if branch == panelBranch or panelBranch in getUpstream(branchReplacedAsterisk):
                positive = True
        if positive:
            panels.append(panel)
    return panels

def getMostSpecificPanel(panels, panelMap):
    mostSpecific = panelMap[panels[0]][0][0]
    mostSpecificLabel = panels[0]
    for i in range(len(panels) - 1):        
        if mostSpecific in getUpstream(panelMap[panels[i+1]][0][0]):
            mostSpecific = panelMap[panels[i+1]][0][0]
            mostSpecificLabel = panels[i+1]
    return mostSpecificLabel

import sys


if len(sys.argv) > 1:
    treeFile = sys.argv[1]
    tabixFilePath = sys.argv[2]
    haplogroupFile = sys.argv[3]
    csvoutrffile = sys.argv[4]
    panelHierFile = sys.argv[5]
    maxThreads = int(sys.argv[6])
    haplogroupClassFilePath = sys.argv[7]
    
    
thetreestuff = parseTreeJSON(treeFile)
hierarchy = thetreestuff[1]

haplogroupClassesJson = json.load(open(haplogroupClassFilePath))
panelMap = {}
for key in haplogroupClassesJson:
    panelMap[key] = [haplogroupClassesJson[key]["branches"],0]
    
writePanelTree(panelMap, panelHierFile)


from Common import CommonMethods

panelHier = CommonMethods.getPanelHier(panelHierFile)
kits = parseKits(haplogroupFile, tabixFilePath, hierarchy, panelMap, panelHier, maxThreads)

headers = parseHeaders(kits)


createCSVforRF(csvoutrffile, headers, kits)

    
