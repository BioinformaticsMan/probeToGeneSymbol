#!/usr/bin/python

import sys, os
import numpy as np

matrixFile = sys.argv[1]
familyFile = sys.argv[2]
outputFile = sys.argv[3]
traitFile  = sys.argv[4]
print (sys.argv)

def readMatrixFile(matrixFile,matrixBegin,sampleID,trait):
    idExprDict={}
    with open (matrixFile,"r") as f_in:
        lines = f_in.readlines()
        for line in lines[matrixBegin:-1]:
            l = line.strip("\n").split("\t")
            ID = l[0].strip("\"")
            exprList = l[1:]
            if ID not in idExprDict:
                idExprDict[ID] = []
                idExprDict[ID].append(exprList)
            else:
                print("ID matched with multiExprValue")
                print(l)
    #headlineList = lines[65].strip("\n").split("\t")
    #headline = "\t".join(headlineList) +"\n"
    headline  = lines[sampleID]
    traitInfo = lines[trait]

    return [idExprDict,headline,traitInfo]

def readFamilyFile(familyFile,tableBegin,tableEnd,geneSymbolCol):
    idGeneDict={}
    with open (familyFile,"r") as f_in:
        lines = f_in.readlines()
        for line in lines[tableBegin:tableEnd]:
        #for line in lines[232:235]:
            l = line.strip("\n").split("\t")
            ID = l[0]
            if len(l) >= 11:
                geneName = "" 
                if "//" in l[geneSymbolCol]:
                    c = l[geneSymbolCol].split(" // ")
                    geneName = c[1]
                    #if "NM_" in c[0]:
                    #    geneName = c[1]
                    #if "ENST" in c[0]:
                    #    geneName = c[1]
            else:
                geneName = ""
            #print (geneName)
            if ID not in idGeneDict:
                if geneName != "":
                    idGeneDict[ID] = geneName
            else:
                continue
                #print ("ID matched with multigene ")
                #print (l)
    return idGeneDict

def removeDulProbes(idExprDict,idGeneDict):
    geneExprDict = {}
    for symbolID in idExprDict:
        #print (symbolID)
        if symbolID in idGeneDict:
            #print(symbolID)
            geneName = idGeneDict[symbolID]
            #print (geneName)
            if geneName not in geneExprDict:
                geneExprDict[geneName] = []
                geneExprDict[geneName].extend(idExprDict[symbolID])
            else:
                geneExprDict[geneName].extend(idExprDict[symbolID])

    for geneID in geneExprDict:
    	i = len(geneExprDict[geneID])
    	exp = np.array(geneExprDict[geneID])
    	exp = exp.astype(float)
    	a = np.sum(exp,axis=0)/i
    	geneExprDict[geneID] = a.tolist()

    return geneExprDict

def outputMatrixFile(outputFile,headline,geneExprDict):
    with open (outputFile, "w") as f_out:
        f_out.write(headline)
        #f_out.write(traitInfo)
        for key in geneExprDict:
    	    content = key + "\t"+"\t".join(str(i) for i in geneExprDict[key])+"\n"
    	    #print (content)
    	    f_out.write(content)

def outputTraitFile(traitFile,headline,traitInfo):
    with open (traitFile, "w") as f_out:
    	f_out.write(headline)
    	f_out.write(traitInfo)

## 01 read the matrix file
print ("Loading matrix file")
l = readMatrixFile(matrixFile,88,87,57)
idExprDict = l[0]
headline = l[1]
traitInfo = l[2]
print ("idExprDict done!")

## 02 read the family file
print ("Loading family file")
idGeneDict = readFamilyFile(familyFile,554,317549,9)
print ("idGeneDict done!")

## 03 remove the dulplicate gene symbol
print ("Removing the dulplicate gene symbol......")
geneExprDict = removeDulProbes(idExprDict,idGeneDict)
print ("Gene symbols are unique")
print (len(idExprDict))
print (len(idGeneDict))
print (len(geneExprDict))

## 04 output matrixFile and traitFile
outputMatrixFile(outputFile,headline,geneExprDict)
outputTraitFile(traitFile,headline,traitInfo)
