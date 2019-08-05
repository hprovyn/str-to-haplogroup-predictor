# -*- coding: utf-8 -*-
"""
Created on Fri May 24 00:23:13 2019

@author: hunte
"""

import json


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
    def getKit(self, tb, theid, clade, panelMap, queue):
        kitToAdd = {}

        if clade != "?\n":
            hg = clade.replace("\n","") #TODO keep * for A super panel
            panels = getPanels(hg, panelMap)
            if len(panels) > 0:
                print(panels)
                panel = getMostSpecificPanel(panels, panelMap)
                print(theid, 'most specific', panel, 'for', clade)                    
                kitToAdd["Haplogroup"] = panel
            else:
                kitToAdd["Haplogroup"] = "?"
        else:
            kitToAdd["Haplogroup"] = "?"
        allowableDownstream = getAllowableDownstream([], kitToAdd["Haplogroup"], panelMap, tb, theid)
        kitToAdd["STRs"] = {}
        strResults = tb.querys("str:" + theid + "-" + theid)
        strs = []
        for strResult in strResults:
            strs.append(strResult[3] + "=" + strResult[4])
        for strallele in strs:
            (thestr, allele) = strallele.split("=")
            kitToAdd["STRs"][thestr] = allele.replace("\\","")
            kitToAdd["allowableDownstream"] = allowableDownstream
        if len(kitToAdd["STRs"]) >= 15:
            queue.put((theid, kitToAdd))
        else:
            queue.put((0,0))

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

import multiprocessing

import time
            
def parseKits(haplogroupFile, tb, hierarchy, panelMap, panelHier):
    subs = {}
    kits = {}
    start = time.time()
    with open(haplogroupFile, 'r') as f:
        theids = []
        thehgs = []
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
        for i in range(len(thehgs)):
            print('processing' + theids[i])
            hg = thehgs[i]
            theid = theids[i]
            kitToAdd = {}

            kitToAdd["Haplogroup"] = mostSpecificPanels[hg]
            allowableDownstream = getAllowableDownstream([], kitToAdd["Haplogroup"], panelMap, tb, theid, panelHier)
            kitToAdd["STRs"] = {}
            strResults = tb.querys("str:" + theid + "-" + theid)
            strs = []
            for strResult in strResults:
                strs.append(strResult[3] + "=" + strResult[4])
            for strallele in strs:
                (thestr, allele) = strallele.split("=")
                strtouse = thestr
                if thestr in subs:
                    strtouse = subs[thestr]
                kitToAdd["STRs"][strtouse] = allele.replace("\\","")
                kitToAdd["allowableDownstream"] = allowableDownstream
            if len(kitToAdd["STRs"]) >= 15:
                kits[theid] = kitToAdd
    end = time.time()
    print("elapsed time processing specific and allowable downstream SNPs" + str(round((end - start) / 60,2)) + " min")
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

panelMap = {"A": [["A00", "A1a", "A1b1", "A1*", "A1b*"],0],
             "B": [["B"],0],
             "C": [["C"],0],
             "D": [["D"],0],
             "E": [["E"],0],
             "E1a-M132": [["E-M132"],0],
             "E1b-L19": [["E-L19"],0],
             "E1b-M191": [["E-M191"],0],
             "E1b-U175": [["E-U175"],0],
             "E1b-V12": [["E-V12"],0],
             "E1b-V13": [["E-V13"],0],
             "E1b-Z834": [["E-Z830"],0],
             "G": [["G"],0],
             "G2a-L497":[["G-L497"],0],
             "G2a-L1259":[["G-L1259"],0],
             "G2a-Z6552":[["G-Z6552"],0],
             "H-L901":[["H"],0],
             "H-M82":[["H-M82"],0],
             "I1":[["I1"],0],
             "I1-L22":[["I-L22"],0],
             "I1-Z63":[["I-Z63"],0],
             "I1-Z140":[["I-Z140"],0],
             "I1-Z2539":[["I-CTS7362"],0],
             "I2":[["I2"],0],
             "I2a-CTS595":[["I-CTS595"],0],
             "I2a-CTS10100":[["I-CTS10057"],0],
             "I2a-M284":[["I-M284"],0],
             "I2a-M423":[["I-M423"],0],
             "J1":[["J1"],0],
             "J2":[["J2"],0],
             "J2a-L24":[["J-L24"],0],
             "J2a-M67":[["J-M67"],0],
             "J2a-PF5197":[["J-PF5197"],0],
             "J2b":[["J-M102"],0],
             "L":[["L"],0],
             "MS-P397":[["K2b1"],0],
             "N":[["N"],0],
             "N1a-VL29":[["N-VL29"],0],
             "N1a-Z1936":[["N-Z1936"],0],
             "O":[["O"],0],
             "O-F100":[["O-M134"],0],
             "O-F145":[["O-L465"],0],
             "O-M268":[["O-M268"],0],
             "Q":[["Q"],0],
             "Q1b-L53":[["Q-L53"],0],
             "Q1b-L275":[["Q-L275"],0],
             "R1a":[["R1a"],0],
             "R1a-L664":[["R-L664"],0],
             "R1a-M458":[["R-M458"],0],
             "R1a-Z93":[["R-Z93"],0],
             "R1a-Z280":[["R-Z280"],0],
             "R1a-Z284":[["R-Z284"],0],
             "R1b-DF19":[["R-DF19"],0],
             "R1b-DF21":[["R-DF21"],0],
             "R1b-DF27":[["R-DF27"],0],
             "R1b-DF41":[["R-CTS2501"],0],
             "R1b-DF49":[["R-DF49"],0],             
             "R1b-DF63":[["R-DF63"],0],
             #"":[["R-FGC5494"],0],
             "R1b-FGC11134":[["R-FGC11134"],0],
             "R1b-L21":[["R-L21"],0],
             "R1b-L48":[["R-L48"],0],
             "R1b-L513":[["R-DF1"],0],
             "R1b-L1335":[["R-L1335"],0],
             "R1b":[["R1b"],0],
             "R1b-S1051":[["R-S1051"],0],
             "R1b-U106":[["R-U106"],0],
             "R1b-U152":[["R-U152"],0],
             "R1b-Z251":[["R-Z251"],0],
             "R1b-Z253":[["R-Z253"],0],
             "R1b-Z255":[["R-Z255"],0],
             "R1b-Z2103":[["R-Z2103"],0],
             "R2":[["R2"],0],
             "T":[["T"],0],
             "T-L131":[["T-L131"],0],
             "T-P77":[["T-P77"],0]}

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
    
thetreestuff = parseTreeJSON(treeFile)
hierarchy = thetreestuff[1]

writePanelTree(panelMap, panelHierFile)

import tabix
tb = tabix.open(tabixFilePath)

from Common import CommonMethods

panelHier = CommonMethods.getPanelHier(panelHierFile)
kits = parseKits(haplogroupFile, tb, hierarchy, panelMap, panelHier)

headers = parseHeaders(kits)


createCSVforRF(csvoutrffile, headers, kits)

    
