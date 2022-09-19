activateErrorFunctions(activate=True,activatefirsttime=True)

""" 
Below each of the functions are tested seperately to see whether "
a) The functions produce the same result when the error rate is 0
b) The values produced when the error rate is used make sense. 
"""
#-----------------------------
'MERGE VECTORS'
# following 2vectors should cover all
vect1WithFlags = [(4, 8, 0.0007, True), (4, 9), (5, 10), (4,11,0.000, True), (4, 12, 0.0004, True), (6, 13, [0.8, 0.2, 0, 0]), (6, 14, [0.8, 0.2, 0, 0]), (4, 1481)]
vect2WithFlags = [(4,7), (2, 8), (4, 9, 0.0008, True), (4, 10, 0.0004, True), (5, 11), (6, 12, [0.8, 0.2, 0, 0]), (0, 13, 0.00071, True), (6, 14, [0.8, 0.2, 0, 0]), (4, 1481)]
vect1WithoutFlags = [(4, 8, 0.0007), (4, 9), (5, 10), (4,11), (4, 12, 0.0004), (6, 13, [0.8, 0.2, 0, 0]), (6, 14, [0.8, 0.2, 0, 0]), (4, 1481)]
vect2WithoutFlags = [(4,7), (2, 8), (4, 9, 0.0008), (4, 10, 0.0004), (5,11), (6, 12, [0.8, 0.2, 0, 0]), (0, 13, 0.00071), (6, 14, [0.8, 0.2, 0, 0]), (4, 1481)]

bLen1 = 0.0001
bLen2 = 0.0002

errorRate=0
resWith = mergeVectorsError(vect1WithFlags, bLen1, vect2WithFlags, bLen2, mutMatrix, returnLK=True, useRateVariation=rateVariation,
                 mutMatrices=mutMatrices, node1isleaf=False, node2isleaf=False)
resWithout = mergeVectorsOriginal(vect1WithoutFlags, bLen1, vect2WithoutFlags, bLen2, mutMatrix, returnLK=True, useRateVariation=rateVariation,
                 mutMatrices=mutMatrices, node1isleaf=False, node2isleaf=False)
assert(areVectorsDifferentDebugging(resWith[0], resWithout[0])==False) #should be false
print(resWith)
print(resWithout)
errorRate=0.05
resWith = mergeVectorsError(vect1WithFlags, bLen1, vect2WithFlags, bLen2, mutMatrix, returnLK=True, useRateVariation=rateVariation,
                 mutMatrices=mutMatrices, node1isleaf=False, node2isleaf=False)
assert(areVectorsDifferentDebugging(resWith[0], resWithout[0])==True) #should be true
print(resWith)
print(resWithout)
#-----------------------------
'MERGE VECTORS UP DOWN'
vect1WithFlags = [(4, 8, 0.0007, 0.0007, True), (4, 9), (5, 10), (4,11,0.0007,0.0007, True), (2, 12, 0.0004, 0.0007, True), (6, 13, [0.8, 0.2, 0, 0]), (6, 14, [0.8, 0.2, 0, 0]), (4, 1481)]
vect2WithFlags = [(4,7), (2, 8), (4, 9, 0.0008, True), (4, 10, 0.0004, True), (5, 11), (6, 12, [0.8, 0.2, 0, 0]), (0, 13, 0.00071, True), (6, 14, [0.8, 0.2, 0, 0]), (4, 1481)]
vect1WithoutFlags = [(4, 8, 0.0007,0.0007), (4, 9), (5, 10), (4,11,0.0007,0.0007), (2, 12, 0.0004,0.0007), (6, 13, [0.8, 0.2, 0, 0]), (6, 14, [0.8, 0.2, 0, 0]), (4, 1481)]
vect2WithoutFlags = [(4,7), (2, 8), (4, 9, 0.0008), (4, 10, 0.0004), (5,11), (6, 12, [0.8, 0.2, 0, 0]), (0, 13, 0.00071), (6, 14, [0.8, 0.2, 0, 0]), (4, 1481)]
errorRate=0
resWith = mergeVectorsUpDownError(vect1WithFlags, bLen1, vect2WithFlags, bLen2, mutMatrix,  useRateVariation=False,
                 mutMatrices=None, node1isleaf=False, node2isleaf=False)
assert(areVectorsDifferentDebugging(resWith, resWithout)==False) #should be false
print(resWith)
print(resWithout)
errorRate=0.05
resWith = mergeVectorsUpDownError(vect1WithFlags, bLen1, vect2WithFlags, bLen2, mutMatrix, useRateVariation=False,
                 mutMatrices=None, node1isleaf=False, node2isleaf=False)
resWithout = mergeVectorsUpDownOriginal(vect1WithoutFlags, bLen1, vect2WithoutFlags, bLen2, mutMatrix, useRateVariation=False,
                 mutMatrices=None, node1isleaf=False, node2isleaf=False)
assert(areVectorsDifferentDebugging(resWith[0], resWithout[0])==True) #should be true
print(resWith)
print(resWithout)
#----------------
'APPEND PROB NODE'
errorRate=0
resWith = appendProbNodeErrorRate(vect1WithFlags, vect2WithFlags, bLen1,  mutMatrix,  useRateVariation=False, mutMatrices=None, node2isleaf=False)
resWithout = appendProbNodeOriginal(vect1WithoutFlags, vect2WithoutFlags, bLen1, mutMatrix,  useRateVariation=False, mutMatrices=None, node2isleaf=False)
assert(round(resWith,7)==round(resWithout,7)) #should be true
print(resWith)
print(resWithout)
errorRate=0.05
resWith = appendProbNodeErrorRate(vect1WithFlags, vect2WithFlags, bLen1, mutMatrix,  useRateVariation=False, mutMatrices=None, node2isleaf=False)
resWithout = appendProbNodeOriginal(vect1WithoutFlags, vect2WithoutFlags, bLen1, mutMatrix,  useRateVariation=False, mutMatrices=None, node2isleaf=False)
assert(round(resWith,7)!=round(resWithout,7)) #should be false
print(resWith) # is -8, whereas without error rate it is -18. In some cases in appendprobnode the LKcost would get more negative (eg. A-->A), but for other cases it can get either more negative or more positive.
print(resWithout)
#----------------
'FIND ROOT PROB'
errorRate=0
print(findProbRootError(vect1WithFlags))
print(findProbRootOriginal(vect1WithoutFlags))
print(findProbRootError(vect1WithFlags) == findProbRootOriginal(vect1WithoutFlags))
errorRate=0.05
print(findProbRootError(vect1WithFlags)) # gets a bit lower, which is what should happen
print(findProbRootOriginal(vect1WithoutFlags))
#----------------
'ESTIMATE BlEN'
errorRate=0
resWith =errorRateEstimateBranchLengthWithDerivative(vect1WithFlags, vect2WithFlags, mutMatrix)
resWithout = estimateBranchLengthWithDerivativeOriginal(vect1WithoutFlags, vect2WithoutFlags, mutMatrix)
assert(round(resWith,7)==round(resWithout,7)) #should be true
print(resWith)
print(resWithout)
errorRate=0.0005
resWith =errorRateEstimateBranchLengthWithDerivative(vect1WithFlags, vect2WithFlags, mutMatrix) #branch length becomes 0. if you let the error rate approach 0, it becomes similar to the previous value (with normal MAPLE functions)
assert(round(resWith,7)!=round(resWithout,7)) #should be false
print(resWith)
print(resWithout)

#Another case to test the branch length derivative.
vect2WithFlags = [(5, 5), (4, 8), (2, 9), (4, 945), (6, 946, [2.9754581065168545e-13, 0.9999990762239054, 2.899054813462519e-16, 9.237757967939945e-07]), (4, 1122), (6, 1123, [1.2048219778560687e-13, 3.4972858641174724e-06, 3.603821044511927e-15, 0.9999965027140117]), (4, 1413), (5, 1414), (4, 1443), (4, 1444, 0.0055191203096563776, False), (5, 1445), (4, 1454), (5, 1455), (4, 1467), (5, 1481)]
vect2WithoutFlags = [(5, 5), (4, 8), (2, 9), (4, 945), (6, 946, [2.9754581065168545e-13, 0.9999990762239054, 2.899054813462519e-16, 9.237757967939945e-07]), (4, 1122), (6, 1123, [1.2048219778560687e-13, 3.4972858641174724e-06, 3.603821044511927e-15, 0.9999965027140117]), (4, 1413), (5, 1414), (4, 1443), (4, 1444, 0.0055191203096563776), (5, 1445), (4, 1454), (5, 1455), (4, 1467), (5, 1481)]
vect1WithFlags =[(5, 25), (4, 1083, 5.055016065236067e-07, 0.0, True), (3, 1084, 5.055016065236067e-07, 0.0, True), (4, 1413, 5.055016065236067e-07, 0.0, True), (5, 1414), (4, 1436, 5.055016065236067e-07, 0.0, True), (5, 1481)]
vect1WithoutFlags =[(5, 25), (4, 1083, 5.055016065236067e-07, 0.0), (3, 1084, 5.055016065236067e-07, 0.0), (4, 1413, 5.055016065236067e-07, 0.0), (5, 1414), (4, 1436, 5.055016065236067e-07, 0.0), (5, 1481)]
errorRate=0
resWith =errorRateEstimateBranchLengthWithDerivative(vect1WithFlags, vect2WithFlags, mutMatrix)
resWithout = estimateBranchLengthWithDerivativeOriginal(vect1WithoutFlags, vect2WithoutFlags, mutMatrix)
assert(round(resWith,7)==round(resWithout,7)) #should be true
print(resWith)
print(resWithout)
errorRate=0.0005
resWith =errorRateEstimateBranchLengthWithDerivative(vect1WithFlags, vect2WithFlags, mutMatrix) #branch length becomes 0. if you let the error rate approach 0, it becomes similar to the previous value (with normal MAPLE functions)
assert(round(resWith,7)!=round(resWithout,7)) #should be false
print(resWith)
print(resWithout)

#Another case to test the branch length derivative.
vect2WithoutFlags = [(5, 25), (4, 1083), (3, 1084), (4, 1413), (5, 1414), (4, 1436), (5, 1481)]
vect2WithFlags = vect2WithoutFlags
vect1WithFlags = [(5, 5), (4, 8, 0.0007112327864252331, 0.0, False), (2, 9, 0.0007112327864252331, 0.0, False), (4, 945, 0.0007112327864252331, 0.0, False), (6, 946, [0.00018386197763338815, 0.9989795877350927, 1.516124784132128e-05, 0.0008213890394325605]), (4, 1122, 0.0007112327864252331, 0.0, False), (6, 1123, [7.08674650514798e-05, 0.0005874520216743095, 3.2378409262364464e-05, 0.9993093021040118]), (4, 1413, 0.0007112327864252331, 0.0, False), (5, 1414), (4, 1443, 0.0007112327864252331, 0.0, False), (4, 1444, 0.006230353096081611, 0.0, False), (5, 1445), (4, 1454, 0.0007112327864252331, 0.0, False), (5, 1455), (4, 1467, 0.0007112327864252331, 0.0, False), (5, 1481)]
vect1WithoutFlags = [(5, 5), (4, 8, 0.0007112327864252331, 0.0), (2, 9, 0.0007112327864252331, 0.0), (4, 945, 0.0007112327864252331, 0.0), (6, 946, [0.00018386197763338815, 0.9989795877350927, 1.516124784132128e-05, 0.0008213890394325605]), (4, 1122, 0.0007112327864252331, 0.0), (6, 1123, [7.08674650514798e-05, 0.0005874520216743095, 3.2378409262364464e-05, 0.9993093021040118]), (4, 1413, 0.0007112327864252331, 0.0), (5, 1414), (4, 1443, 0.0007112327864252331, 0.0), (4, 1444, 0.006230353096081611, 0.0), (5, 1445), (4, 1454, 0.0007112327864252331, 0.0), (5, 1455), (4, 1467, 0.0007112327864252331, 0.0), (5, 1481)]
errorRate=0
resWithout = estimateBranchLengthWithDerivativeOriginal(vect1WithoutFlags, vect2WithoutFlags, mutMatrix)
resWith =errorRateEstimateBranchLengthWithDerivative(vect1WithFlags, vect2WithFlags, mutMatrix, node2isleaf=True)
assert(round(resWith,7)==round(resWithout,7)) #should be true
print(resWith)
print(resWithout)
errorRate=0.0005
resWith =errorRateEstimateBranchLengthWithDerivative(vect1WithFlags, vect2WithFlags, mutMatrix, node2isleaf=True) #branch length becomes 0. if you let the error rate approach 0, it becomes similar to the previous value (with normal MAPLE functions)
assert(round(resWith,7)!=round(resWithout,7)) #should be false
print(resWith)
print(resWithout)
#---------------------------------------------------------------


"""
The following tests that when creating a tree with flags and errorrate functions, but an error rate of 0, and optimizing this,
the result will be the same as a tree without the incorporation of flags or usuage of error rate functions
"""
import copy
errorRate = 0
reCalculateAllGenomeLists(t1,mutMatrix)
traverseTreeToOptimizeBranchLengths(t1,mutMatrix,testing=False,useRateVariation=rateVariation,mutMatrices=mutMatrices)
t0_1 = copy.deepcopy(t1)
reCalculateWithErrors(t1,mutMatrix, errorRate,useRateVariation=rateVariation,mutMatrices=mutMatrices,firstTimeError=True)
traverseTwoTopologies(t1,t0_1) # mergeVectors and reCalculateAllGenomeLists seem to work.
activateErrorFunctions(activate=False)
improvement, best_lengths_t0 = traverseTreeToOptimizeBranchLengths(t0_1,mutMatrix,testing=True,useRateVariation=rateVariation,mutMatrices=mutMatrices)
activateErrorFunctions(activate=True)
test_if_different=True
improvement, best_lengths_t1 = traverseTreeToOptimizeBranchLengths(t1,mutMatrix,testing=True,useRateVariation=rateVariation,mutMatrices=mutMatrices) # Initial border parameters don't fit expectations
traverseTwoTopologies(t1,t0_1)
test_if_different = False

# ------------------------------
"""
The following code snippet creates a tree with errors, and compares certain metrics with the corresponding tree without the incorporation of error rates, e.g.:
a) the number different types of entries
b) the sum of all branch lengths in the tree
c) the number of flags stored and set to true. 
Don't run this code when running also the above piece. 
"""

import copy
t0 = copy.deepcopy(t1)
reCalculateAllGenomeLists(t1,mutMatrix)
traverseTreeToOptimizeBranchLengths(t1,mutMatrix,testing=False,useRateVariation=rateVariation,mutMatrices=mutMatrices)
countEntriesAll(t1)
oldTotBLen = counttotBLenAll(t1)
t0_1 = copy.deepcopy(t1)
oldTotBLen = counttotBLenAll(t0_1)
reCalculateAllGenomeLists(t0_1, mutMatrix,useRateVariation=rateVariation,mutMatrices=mutMatrices,)
activateErrorFunctions(activate=True, activatefirsttime=True) #"start recalculating with ErrorRate"
errorRate = 0.0005
reCalculateWithErrors(t1, mutMatrix, errorRate, useRateVariation=rateVariation, mutMatrices=mutMatrices, firstTimeError=True)
print("ErrorRate included")
countEntriesAll(t1) #checking if node compositions remain the same. .. Difficult to conlcude anything based on this, since node.dist is often False.
traverseTreeToOptimizeBranchLengths(t1,mutMatrix,testing=False,useRateVariation=rateVariation,mutMatrices=mutMatrices)
newTotBLen = counttotBLenAll(t1)
print('diff tot blen '+ str(newTotBLen - oldTotBLen))
print( 'expected dif = errorRate*nsamples/subtitutionrate ' + str(errorRate*178) ) #the number of sample nodes are 178 in the case I have been testing for. This was the change I would expect in branch lengths.
countFlagsAll(t1)

#----------------------------------------------
"""Shorten MAPLE genome length """
def shortenVect(vect, pos=None, direction=1):
    if not pos:
        pos = vect[-1][1]/2 #halve size original
    pos = int(pos-1)
    newVect =[]
    for entry in vect[::direction]:
        if entry[1] < pos:
            newVect.append(entry)
        else:
            entry= list(entry)
            entry[1] = pos
            newVect.append(tuple(entry))
            return newVect

def shortenGenomeLengthNode(node, pos=None, direction=1):
    for vects in [ 'probVect',  'probVectTotUp',  'probVectUpLeft',  'probVectUpRight']:
        if hasattr(node, vects):
            vect =getattr(node, vects)
            if vect:
                vect = shortenVect(vect, pos=pos, direction=direction)
                setattr(node, vects, vect)
traverseTopology(t3,shortenGenomeLengthNode)

def shortenGenomeLengthFile(mapleFile, lRef=None):
