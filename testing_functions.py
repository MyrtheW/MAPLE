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
# def shortenVect(vect, pos=None, direction=1):
#     if not pos:
#         pos = vect[-1][1]/2 #halve size original
#     pos = int(pos-1)
#     newVect =[]
#     for entry in vect[::direction]:
#         if entry[1] < pos:
#             newVect.append(entry)
#         else:
#             entry= list(entry)
#             entry[1] = pos
#             newVect.append(tuple(entry))
#             return newVect
#
# def shortenGenomeLengthNode(node, pos=None, direction=1):
#     for vects in [ 'probVect',  'probVectTotUp',  'probVectUpLeft',  'probVectUpRight']:
#         if hasattr(node, vects):
#             vect =getattr(node, vects)
#             if vect:
#                 vect = shortenVect(vect, pos=pos, direction=direction)
#                 setattr(node, vects, vect)
# traverseTopology(t3,shortenGenomeLengthNode)
#
#
#


#--------------------------------------------

def findBestParentTopologyDebugging(appendedToNode, bestBranchLengths, node,child,bestLKdiff,removedBLen,mutMatrix,compensanteForBranchLengthChange =False, strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,useRateVariation=False,mutMatrices=None):
    """Force a placement of the 'node' to the appendedToNode. Force it with bestBranchLengths Tuple"""
    #currentLK = bestLKdiff
    bestNode=node
    bestNodes=[]
    nodesToVisit=[]
    removedPartials=node.children[child].probVect
    removedPartialsIsLeaf = (node.children[child].children==[])
    if node.up!=None:
        if node.up.children[0]==node:
            childUp=1
            vectUpUp=node.up.probVectUpRight
        else:
            childUp=2
            vectUpUp=node.up.probVectUpLeft
        #the list nodesToVisit keeps trak of the nodes of the tree we wil need to traverse, (first element of each entry),
        # the direction from where we are visitng them (0=from parent, and 1,2=from one of the children),
        # an updated genome list from the direction where we come from (taking into account the removal of the given subtree),
        # a branch length value separating the node from this updated genome list (useful for the fact that the removal of the subtree changes the branch length at the removal node),
        # a flag that says if the updated genome list passed needs still updating, or if it has become identical to the pre-existing genome list in the tree (which usually happens after a while),
        # the likelihood cost of appending at the last node encountered in this direction,
        # a number of consecutively failed traversal steps since the last improvement found (if this number goes beyond a threshold, traversal in the considered direction might be halted).
        nodesToVisit.append((node.up,childUp,node.children[1-child].probVect,node.children[1-child].dist+node.dist,True,bestLKdiff,0, node.children[1-child].children==[]))
        nodesToVisit.append((node.children[1-child],0,vectUpUp,node.children[1-child].dist+node.dist,True,bestLKdiff,0, False))
    else:
        # case node is root
        if node.children[1-child].children: # case there is only one sample outside of the subtree doesn't need to be considered
            child1=node.children[1-child].children[0]
            child2=node.children[1-child].children[1]
            vectUp1=rootVector(child2.probVect,child2.dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,
                               isLeaf=child2.children==[])
            nodesToVisit.append((child1,0,vectUp1,child1.dist,True,bestLKdiff,0, False)) #isLeaf is false for vectup1
            vectUp2=rootVector(child1.probVect,child1.dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,
                               isLeaf=child1.children==[])
            nodesToVisit.append((child2,0,vectUp2,child2.dist,True,bestLKdiff,0, False)) #isLeaf is false

    while nodesToVisit:
        t1,direction,passedPartials,distance,needsUpdating,lastLK,failedPasses, passedPartialsIsLeaf=nodesToVisit.pop()
        if direction==0:
            #consider the case we are moving from a parent to a child
            if t1.dist and (not (t1.up==node or t1.up==None)):
                if needsUpdating:
                    midTot=mergeVectorsUpDown(passedPartials,distance/2,t1.probVect,distance/2,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,
                                              node1isleaf=False, #Because it is a vectup, mergeVectorsUpDown
                                              node2isleaf=t1.children==[])

                    if not areVectorsDifferent(midTot,t1.probVectTotUp):
                        needsUpdating=False
                else:
                    midTot=t1.probVectTotUp
                if midTot==None:
                    continue


            #keep crawling down into children nodes unless the stop criteria for the traversal are satisfied.

            if t1.children:
                #if (failedPasses<=allowedFailsTopology or nodeProb>(bestLKdiff-thresholdLogLKtopology) ) and len(t1.children)==2:
                child=t1.children[0]
                otherChild=t1.children[1]
                if needsUpdating:
                    vectUpRight=mergeVectorsUpDown(passedPartials,distance,otherChild.probVect,otherChild.dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,
                                                   node2isleaf=otherChild.children==[] )
                else:
                    vectUpRight=t1.probVectUpRight
                if vectUpRight!=None:
                    nodesToVisit.append((child,0,vectUpRight,child.dist,needsUpdating,bestLKdiff,failedPasses, False))
                child=t1.children[1]
                otherChild=t1.children[0]
                if needsUpdating:
                    vectUpLeft=mergeVectorsUpDown(passedPartials,distance,otherChild.probVect,otherChild.dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,
                                                  node2isleaf=otherChild.children == [])
                else:
                    vectUpLeft=t1.probVectUpLeft
                if vectUpLeft!=None:
                    nodesToVisit.append((child,0,vectUpLeft,child.dist,needsUpdating,bestLKdiff,failedPasses, False)) #no leaf

        else: #case when crawling up from child to parent
            otherChild=t1.children[2-direction]
            midBottom=None
            if t1.dist and t1.up!=None: #try appending mid-branch
                if needsUpdating:
                    midBottom=mergeVectors(otherChild.probVect,otherChild.dist,passedPartials,distance,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,
                                           node1isleaf=otherChild.children==[], node2isleaf= passedPartialsIsLeaf)
                    if midBottom==None:
                        continue
                    if t1==t1.up.children[0]:
                        vectUp=t1.up.probVectUpRight
                    else:
                        vectUp=t1.up.probVectUpLeft
                    midTot=mergeVectorsUpDown(vectUp,t1.dist/2,midBottom,t1.dist/2,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,
                                              node2isleaf=False) #midBottom is not from a leaf node! it has just been created above.
                    if not areVectorsDifferent(midTot,t1.probVectTotUp):
                        needsUpdating=False
                else:
                    midTot=t1.probVectTotUp
                if midTot==None:
                    continue



            if True:
                # keep crawling up into parent and sibling node
                if t1.up!=None: #case the node is not the root
                    #first pass the crawling down the other child
                    if t1==t1.up.children[0]:
                        upChild=0
                        if needsUpdating:
                            vectUpUp=t1.up.probVectUpRight
                    else:
                        upChild=1
                        if needsUpdating:
                            vectUpUp=t1.up.probVectUpLeft
                    if needsUpdating:
                        vectUp=mergeVectorsUpDown(vectUpUp,t1.dist,passedPartials,distance,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices, node2isleaf=passedPartialsIsLeaf)
                    else:
                        if direction==1:
                            vectUp=t1.probVectUpLeft
                        else:
                            vectUp=t1.probVectUpRight

                    if vectUp==None:
                        continue
                    else:
                        nodesToVisit.append((otherChild,0,vectUp,otherChild.dist,needsUpdating,bestLKdiff,failedPasses, False))
                    #now pass the crawling up to the parent node
                    if needsUpdating:
                        if midBottom==None:
                            midBottom=mergeVectors(otherChild.probVect,otherChild.dist,passedPartials,distance,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,
                                                   node1isleaf=otherChild.children==[], node2isleaf=passedPartialsIsLeaf )#, node1isleaf= node1isleaf,node2isleaf=node2isleaf )
                            if midBottom==None:
                                continue
                    else:
                        midBottom=t1.probVect
                    nodesToVisit.append((t1.up,upChild+1,midBottom,t1.dist,needsUpdating,bestLKdiff,failedPasses, False))
                #now consider case of root node
                else:
                    if needsUpdating:
                        vectUp=rootVector(passedPartials,distance,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices, isLeaf=passedPartialsIsLeaf)
                    else:
                        if direction==1:
                            vectUp=t1.probVectUpLeft
                        else:
                            vectUp=t1.probVectUpRight
                    nodesToVisit.append((otherChild,0,vectUp,otherChild.dist,needsUpdating,bestLKdiff,failedPasses, False))

    #Initial exploration is finished.
    #Now, for each branch within threshold likelihood distance from the best found, optimize branch lengths.
    #Use optimized scores to select final best branch

    # IMPROVMENT ESTIMATE OF NEW NODE
    t1=appendedToNode

    if t1==t1.up.children[0]:
        upVect=t1.up.probVectUpRight
    else:
        upVect=t1.up.probVectUpLeft
    downVect=t1.probVect
    distance=t1.dist
    #midTot=t1.probVectTotUp
    downVectIsLeaf = (t1.children==[])


    bestTopLength,bestBottomLength,bestAppendingLength = bestBranchLengths
    #now optimize appending location
    # midLowerVector=mergeVectors(downVect,distance/2,removedPartials,bestAppendingLength,mutMatrix,
    #                             node1isleaf=downVectIsLeaf,node2isleaf=removedPartialsIsLeaf) # When changing node1isleaf or node2isleaf, it remains 0.
    # midTopVector=mergeVectorsUpDown(upVect,bestTopLength,removedPartials,bestAppendingLength,mutMatrix,
    #                                 node2isleaf=removedPartialsIsLeaf)
    newMidVector=mergeVectorsUpDown(upVect,bestTopLength,downVect,bestBottomLength,mutMatrix,node2isleaf=downVectIsLeaf)
    appendingCost=appendProbNode(newMidVector,removedPartials,bestAppendingLength,mutMatrix, node2isleaf=removedPartialsIsLeaf) #shouldn't rate variation be used here and the original one?

    if compensanteForBranchLengthChange:  # if wanted, do a more thorough examination of the appending cost, taking into account a possible change in branch length of the branch on which to be appended.
        initialCost = appendProbNode(upVect, downVect, distance, mutMatrix, node2isleaf=downVectIsLeaf)
        newPartialCost = appendProbNode(upVect, downVect, bestBottomLength + bestTopLength, mutMatrix,
                                        node2isleaf=downVectIsLeaf)
        appendingCost = appendingCost + newPartialCost - initialCost
    return bestNode, appendingCost, bestBranchLengths


def traverseTopology3(node, leafsOnly=False):
    """
    Traverse a topology and perform a function
    :param node: the root node.
    :param function: the function to be performed at every node. It should take only 'node' as argument
    :param leafsOnly: If True, the function will only be performed at the leaf node
    """
    nodesToVisit=[node]
    nodesToVisit2 =[]
    while nodesToVisit:
        newNode=nodesToVisit.pop()
        for c in newNode.children:
            nodesToVisit.append(c)
        if not newNode.children or not leafsOnly:
            nodesToVisit2.append(newNode)
    return nodesToVisit2

def getRoot(tree):
    root = tree
    while root.up != None:
        root = root.up
    return root

def shortenVect(vect):
    # if not begin_pos:
    #     begin_pos = 916#349 #vect[-1][1]/10 #halve size original
    # if not end_pos:
    #     end_pos = 919 #351#lRef
    # # begin_pos = int(begin_pos-0.45)
    # # end_pos= int(end_pos-0.45)
    newVect =[(4, begin_pos-1)] if begin_pos -1 >0 else []
    for entry in vect:
        if (entry[1] < begin_pos):
            continue
        elif (entry[1] > end_pos):
            while len(newVect):
                if newVect[-1][0] ==4:
                    newVect.pop()
                else:
                    break
            newVect.append((4,lRef))
            return newVect
        else:
            newVect.append(entry)
    return newVect



def shortenGenomeLengthNode(node, pos=None, direction=1):
    for vects in [ 'probVect',  'probVectTot', 'probVectTotUp',  'probVectUpLeft',  'probVectUpRight']:
        if hasattr(node, vects):
            vect =getattr(node, vects)
            if vect:
                vect = shortenVect(vect)
                setattr(node, vects, vect)
                # if not node.children:
                #     print(node.name)
                #     print(vect)

traverseTreeToOptimizeBranchLengths(getRoot(t1), mutMatrix, mutMatrices=mutMatrices) #remove. at the moment to prevent the branch length leading to S121 to become 1.

"""TESTING CODE: CREATING A WORSE TREE AND HOPING FOR SPR DOWNSTREAM"""
import copy
import numpy as np
t0 = copy.deepcopy(t1)



def shortenSPR():
    t1=copy.deepcopy(t0)
    traverseTopology(t1, shortenGenomeLengthNode)
    reCalculateAllGenomeLists(getRoot(t1), mutMatrix, checkExistingAreCorrect=False, useRateVariation=False,mutMatrices=mutMatrices)  # remove               #print("Post-SPR tree: "+createBinaryNewick(root))
    traverseTreeToOptimizeBranchLengths(getRoot(t1), mutMatrix, mutMatrices=mutMatrices) #remove. at the moment to prevent the branch length leading to S121 to become 1.
    reCalculateAllGenomeLists(getRoot(t1), mutMatrix, checkExistingAreCorrect=False, useRateVariation=False,mutMatrices=mutMatrices)  # remove               #print("Post-SPR tree: "+createBinaryNewick(root))

    branchLengths = (1e-8, 1e-8, 1e-8)
    beforeForcedSPRLK = calculateTreeLikelihood(getRoot(t1), mutMatrix, useRateVariation=False, mutMatrices=mutMatrices)
    for seed in [1]:# range(1,100): #[80]: #
        t2= copy.deepcopy(t1)
        # recalculate LK --> check if correct
        # select some node as best node
        nodeList = traverseTopology3(t2, leafsOnly=True)
        np.random.seed(seed); subTree= np.random.choice(nodeList) # select some subtree
        np.random.seed(200-seed); appendedToNode = np.random.choice(nodeList)
        subTreeNodeList = traverseTopology3(subTree)  # remove. at the moment to prevent the branch length leading to S121 to become 1.
        if appendedToNode in subTreeNodeList: # subTree should not be appended to a node in the subTree
            continue
        newRoot = cutAndPasteNode(subTree, appendedToNode, branchLengths, 0, mutMatrix, useRateVariation=False, mutMatrices=mutMatrices)
        root = getRoot(t2)
        # reCalculateAllGenomeLists(root, mutMatrix, checkExistingAreCorrect=False, useRateVariation=False, mutMatrices=mutMatrices)
        # traverseTreeToOptimizeBranchLengths(root, mutMatrix, mutMatrices=mutMatrices)  # remove. at the moment to prevent the branch length leading to S121 to become 1.
        reCalculateAllGenomeLists(root, mutMatrix, checkExistingAreCorrect=True, useRateVariation=False, mutMatrices=mutMatrices)  # remove               #print("Post-SPR tree: "+createBinaryNewick(root))
        oldTreeLK = calculateTreeLikelihood(root, mutMatrix, useRateVariation=False, mutMatrices=mutMatrices)

        parentNode = subTree.up
        if parentNode.children[0] == subTree:
            child = 0
            vectUp = parentNode.probVectUpRight
        else:
            child = 1
            vectUp = parentNode.probVectUpLeft
        bestCurrenBLen = subTree.dist

        #ORIGINAL appendprobnode value   #optionally optimize Blen before: #when this is not used newLK is positive for the first case, does that make sense?
        #bestCurrenBLen=estimateBranchLengthWithDerivative(vectUp,subTree.probVect,mutMatrix,useRateVariation=rateVariation,mutMatrices=mutMatrices, node2isleaf=(subTree.children==[]))
        originalLK=appendProbNode(vectUp,subTree.probVect,bestCurrenBLen,mutMatrix,useRateVariation=rateVariation,mutMatrices=mutMatrices, node2isleaf=(subTree.children==[]))
        # NEW appendprobnode value
        bestNodeSoFar, newLK, bestBranchLengths = findBestParentTopology(parentNode, child, originalLK, bestCurrenBLen, mutMatrix, strictTopologyStopRules=strictTopologyStopRules, allowedFailsTopology=allowedFailsTopology, thresholdLogLKtopology=thresholdLogLKtopology, useRateVariation=rateVariation, mutMatrices=mutMatrices)# search for SPR move on that node at 'subTree'
        newRoot = cutAndPasteNode(subTree, bestNodeSoFar, bestBranchLengths, newLK, mutMatrix,useRateVariation=False, mutMatrices=mutMatrices)

        root= getRoot(t2)
        reCalculateAllGenomeLists(root,mutMatrix, checkExistingAreCorrect=True,useRateVariation=False,mutMatrices=mutMatrices) #remove               #print("Post-SPR tree: "+createBinaryNewick(root))
        newTreeLK=calculateTreeLikelihood(root,mutMatrix,useRateVariation=False,mutMatrices=mutMatrices)
        print(seed)
        print(subTree)
        print(appendedToNode)
        if appendedToNode != bestNodeSoFar:
            print("subtree is appeneded to other node than that it was originally taken from:")
            print(bestNodeSoFar)
        print("actual improvement: " + str(newTreeLK-oldTreeLK))
        print("supposed improvement: " + str(newLK-originalLK))
        if abs((newTreeLK-oldTreeLK) - (newLK-originalLK))>1:
            print("diference is bigger than 1.0, either the actual improvment is smaller or larger than the supposed improvemnet")
            print((newTreeLK-oldTreeLK) - (newLK-originalLK))
            if abs((newTreeLK - oldTreeLK) - (newLK - originalLK)) > 5:
                print(
                    "diference is bigger than 5.0, either the actual improvment is smaller or larger than the supposed improvemnet")
            if (newTreeLK - oldTreeLK) + 1 < (newLK - originalLK):
                print("diference is bigger than 1.0, either the actual improvment is smaller or larger than the supposed improvemnet")
        if newTreeLK + 1 < beforeForcedSPRLK :
            print("LK has gotten worse than the original tree: new, old:")
            print(newTreeLK)
            print(beforeForcedSPRLK)
        #check if the LK of the tree is equal or better to that before doing the SPR.
        return (newTreeLK-oldTreeLK) - (newLK-originalLK)


# shortening the genomelists notes: when also updating the branch lengths.
start = 1
end = lRef
begin_pos = start
end_pos = end
old_discrepency = shortenSPR()
while end - start - 1 > 0:
    middle = int((-start + end) / 2)

    # left
    begin_pos = start
    end_pos = middle
    left_discrepency = shortenSPR()

    # right
    begin_pos = middle
    end_pos = end
    right_discrepency = shortenSPR()

    if abs(left_discrepency) < abs(right_discrepency):
        new_discrepency = right_discrepency
        start = middle
    else:
        new_discrepency = left_discrepency
        end = middle
    if abs(new_discrepency) < abs(0.5 * old_discrepency):
        print('discrepency got smaller')
    old_discrepency = new_discrepency

#with errors:
# In cutAndPasteNode() removing subtree from the tree, subtree root partials:
# [(5, 25), (4, 57), (1, 58), (4, 80), (3, 81), (4, 179), (0, 180), (4, 186), (1, 187), (4, 221), (0, 222), (4, 306), (1, 307), (4, 348), (3, 349), (4, 456), (1, 457), (4, 690), (3, 691), (4, 855), (2, 856), (4, 945), (3, 946), (4, 960), (1, 961), (4, 1266), (0, 1267), (4, 1413), (5, 1414), (4, 1436), (5, 1481)]
# likelihoods to which it is attached:
# [(5, 5), (4, 8), (2, 9), (4, 57), (1, 58), (4, 80), (3, 81), (4, 179), (0, 180), (4, 186), (1, 187), (4, 221), (0, 222), (4, 306), (1, 307), (4, 348), (6, 349, [1.2336139915254472e-05, 0.49954358563728135, 6.557106703654755e-07, 0.500443422512133]), (4, 456), (1, 457), (4, 520), (6, 521, [0.5003096937558481, 3.58268724681325e-06, 0.4996825891169774, 4.134439927645592e-06]), (4, 615), (6, 616, [7.25390891247358e-08, 0.9999981275949229, 6.96463358705483e-10, 1.7991695246199856e-06]), (4, 690), (3, 691), (4, 855), (2, 856), (4, 861), (6, 862, [7.253949339448537e-08, 0.9999981107395971, 6.959099467408271e-10, 1.8160249994914466e-06]), (4, 945), (6, 946, [1.3817488687477626e-07, 0.0033288099731677593, 1.098360316322003e-08, 0.9966710408683422]), (4, 960), (1, 961), (4, 1041), (6, 1042, [0.00029289425574509587, 3.1347873941137613e-09, 0.9997070935367032, 9.072764323561272e-09]), (4, 1155), (6, 1156, [0.5003096937558481, 3.58268724681325e-06, 0.4996825891169774, 4.134439927645592e-06]), (4, 1266), (6, 1267, [0.5003096937558481, 3.58268724681325e-06, 0.4996825891169774, 4.134439927645592e-06]), (4, 1413), (5, 1414), (4, 1443), (4, 1444, 0.002062316386627106, False), (5, 1445), (4, 1454), (5, 1455), (4, 1467), (5, 1481)]
# In cutAndPasteNode() removing subtree from the tree, subtree root partials:
# [(5, 25), (4, 57), (1, 58), (4, 80), (3, 81), (4, 179), (0, 180), (4, 186), (1, 187), (4, 221), (0, 222), (4, 306), (1, 307), (4, 348), (3, 349), (4, 456), (1, 457), (4, 690), (3, 691), (4, 855), (2, 856), (4, 945), (3, 946), (4, 960), (1, 961), (4, 1266), (0, 1267), (4, 1413), (5, 1414), (4, 1436), (5, 1481)]
# likelihoods to which it is attached:
# [(5, 5), (4, 8, 0.001399137536186032, False), (2, 9, 0.001399137536186032, False), (4, 17, 0.001399137536186032, False), (4, 25, 1e-08, False), (4, 57), (1, 58), (4, 80), (3, 81), (4, 179), (0, 180), (4, 186), (1, 187), (4, 221), (0, 222), (4, 306), (1, 307), (4, 348), (6, 349, [6.062244982425829e-07, 0.9983817858667566, 5.050322006704491e-08, 0.001617557405525141]), (4, 456), (1, 457), (4, 690), (3, 691), (4, 855), (2, 856), (4, 915), (6, 916, [3.391654845289658e-07, 0.9975650044128566, 1.565511697507912e-07, 0.002434499870489097]), (4, 918), (6, 919, [6.062244982425829e-07, 0.9983817858667566, 5.050322006704491e-08, 0.001617557405525141]), (4, 960), (1, 961), (4, 1210), (6, 1211, [0.013364120438165655, 0.9866229268446985, 1.1453454279424498e-05, 1.4992628563336996e-06]), (4, 1413), (5, 1414), (4, 1436), (4, 1443, 1e-08, False), (4, 1444, 0.001399137536186032, False), (5, 1445), (4, 1453, 1e-08, False), (4, 1454, 0.001399137536186032, False), (5, 1455), (4, 1467, 0.001399137536186032, False), (5, 1481)]
# 80
# S104
# S106
# subtree is appeneded to other node than that it was originally taken from:
# actual improvement: 25.406066410150743
# supposed improvement: 25.427065203938383
# LK has gotten worse than the original tree: new, old:
# -4629.5528664597105
# -4629.55282513119

# without errors:
# In cutAndPasteNode() removing subtree from the tree, subtree root partials:
# [(5, 25), (4, 57), (1, 58), (4, 80), (3, 81), (4, 179), (0, 180), (4, 186), (1, 187), (4, 221), (0, 222), (4, 306), (1, 307), (4, 348), (3, 349), (4, 456), (1, 457), (4, 690), (3, 691), (4, 855), (2, 856), (4, 945), (3, 946), (4, 960), (1, 961), (4, 1266), (0, 1267), (4, 1413), (5, 1414), (4, 1436), (5, 1481)]
# likelihoods to which it is attached:
# [(5, 5), (4, 8), (2, 9), (4, 57), (1, 58), (4, 80), (3, 81), (4, 179), (0, 180), (4, 186), (1, 187), (4, 221), (0, 222), (4, 306), (1, 307), (4, 348), (6, 349, [1.2336139915254472e-05, 0.49954358563728135, 6.557106703654755e-07, 0.500443422512133]), (4, 456), (1, 457), (4, 520), (6, 521, [0.5003096937558481, 3.58268724681325e-06, 0.4996825891169774, 4.134439927645592e-06]), (4, 690), (3, 691), (4, 855), (2, 856), (4, 945), (6, 946, [5.5512477129951514e-08, 0.0015858533091860345, 6.555562236380453e-09, 0.9984140846227745]), (4, 960), (1, 961), (4, 1041), (6, 1042, [0.0008408745283214221, 6.82100407586884e-09, 0.9991591070056989, 1.164497566481238e-08]), (4, 1155), (6, 1156, [0.5003096937558481, 3.58268724681325e-06, 0.4996825891169774, 4.134439927645592e-06]), (4, 1266), (6, 1267, [0.5003096937558481, 3.58268724681325e-06, 0.4996825891169774, 4.134439927645592e-06]), (4, 1413), (5, 1414), (4, 1443), (4, 1444, 0.002062316386627106), (5, 1445), (4, 1454), (5, 1455), (4, 1467), (5, 1481)]
# In cutAndPasteNode() removing subtree from the tree, subtree root partials:
# [(5, 25), (4, 57), (1, 58), (4, 80), (3, 81), (4, 179), (0, 180), (4, 186), (1, 187), (4, 221), (0, 222), (4, 306), (1, 307), (4, 348), (3, 349), (4, 456), (1, 457), (4, 690), (3, 691), (4, 855), (2, 856), (4, 945), (3, 946), (4, 960), (1, 961), (4, 1266), (0, 1267), (4, 1413), (5, 1414), (4, 1436), (5, 1481)]
# likelihoods to which it is attached:
# [(5, 5), (4, 8, 0.001399137536186032), (2, 9, 0.001399137536186032), (4, 17, 0.001399137536186032), (4, 25, 1e-08), (4, 57), (1, 58), (4, 80), (3, 81), (4, 179), (0, 180), (4, 186), (1, 187), (4, 221), (0, 222), (4, 306), (1, 307), (4, 456), (1, 457), (4, 690), (3, 691), (4, 855), (2, 856), (4, 915), (1, 916), (4, 960), (1, 961), (4, 1210), (1, 1211), (4, 1413), (5, 1414), (4, 1436), (4, 1443, 1e-08), (4, 1444, 0.001399137536186032), (5, 1445), (4, 1453, 1e-08), (4, 1454, 0.001399137536186032), (5, 1455), (4, 1467, 0.001399137536186032), (5, 1481)]
# 80
# S104 #--> normal, no O entries. origanlly a dist of 0 .
# S106 #--> normal, no O entries.
# subtree is appeneded to other node than that it was originally taken from: [S92,]. Has 1 O entry, and one Blen entry.
# actual improvement: 88.8849341418927
# supposed improvement: 90.41477542630506 #between length 347 and 551 there is still a difference, with small lengths, not specifically we have:
# likelihoods to which it was attached before falsely cutting (doesnt matter so much) 6, 349, [1.2336139915254472e-05, 0.49954358563728135, 6.557106703654755e-07, 0.500443422512133]
# likelihood probvect at the subtree  (3, 349)
# to which it was attached (4, 1481)
# to which it is is reattached. (3, 349) (the same as our subtree).

# -3992.2089306319276
# -3992.208898233922

"""TESTING CODE: FORCING A CUT AND PASTE NODE MOVEMENT, AND CALCULATING THE EXPECTED AND ACTUAL IMPROVEMENT OR DECREASE"""
# problem of this is that I had adapted Nicola's original functions to achieve this. To rule out making any possible mistakes by changes I make,

 # remove               #print("Post-SPR tree: "+createBinaryNewick(root))
for seed in range(1,100):

    t2= copy.deepcopy(t1)
    # recalculate LK --> check if correct
    # select some node as best node
    nodeList = traverseTopology3(t2, leafsOnly=True)
    np.random.seed(seed); subTree= np.random.choice(nodeList) # select some subtree
    np.random.seed(200-seed); appendedToNode = np.random.choice(nodeList)
    subTreeNodeList = traverseTopology3(subTree)  # remove. at the moment to prevent the branch length leading to S121 to become 1.
    if appendedToNode in subTreeNodeList: # subTree should not be appended to a node in the subTree
        continue

    bestBranchLengths = (appendedToNode.dist/2, appendedToNode.dist/2, subTree.dist) # select some branch lengths as best branch lengths
    if not subTree.dist or not appendedToNode.dist: # should at least not both have 0 distance, disturbs branch length part.
        bestBranchLengths = (0.00000001, 0.00000001, 0.00000001) #todo --> this induces the LK difference to rise above 0.# continue

    root = getRoot(t1)
    reCalculateAllGenomeLists(root,mutMatrix, checkExistingAreCorrect=True,useRateVariation=False,mutMatrices=mutMatrices) #remove               #print("Post-SPR tree: "+createBinaryNewick(root))
    oldTreeLK=calculateTreeLikelihood(root,mutMatrix,useRateVariation=False,mutMatrices=mutMatrices)

    # calculate improvement

    #ORIGINAL appendprobnode value
    vectUpSubTree = subTree.up.probVectUpLeft if subTree == subTree.up.children[1] else subTree.up.probVectUpRight #todo check if this is the right vector to use
    originalLK=appendProbNode(vectUpSubTree,subTree.probVect,subTree.dist,mutMatrix,useRateVariation=False,mutMatrices=mutMatrices, node2isleaf=(subTree.children==[]))

    #NEW appendprobnode value
    parentNode = subTree.up
    if parentNode.children[0] == subTree:
        child = 0
        vectUp = parentNode.probVectUpRight
    else:
        child = 1
        vectUp = parentNode.probVectUpLeft
    bestNode, newLK, bestBranchLengths = findBestParentTopologyDebugging(appendedToNode, bestBranchLengths, parentNode, child,originalLK,subTree.dist,mutMatrix,
                                                                         compensanteForBranchLengthChange=True,strictTopologyStopRules=True,allowedFailsTopology=1,thresholdLogLKtopology=40,useRateVariation=False,mutMatrices=mutMatrices)
    #what does allowsfailledtopology mean?
    LK=0 #so that it will be placed where I want it.
    newRoot = cutAndPasteNode(subTree, appendedToNode, bestBranchLengths, LK, mutMatrix,useRateVariation=False, mutMatrices=mutMatrices)

    root= getRoot(t2)
    reCalculateAllGenomeLists(root,mutMatrix, checkExistingAreCorrect=True,useRateVariation=False,mutMatrices=mutMatrices) #remove               #print("Post-SPR tree: "+createBinaryNewick(root))
    newTreeLK=calculateTreeLikelihood(root,mutMatrix,useRateVariation=False,mutMatrices=mutMatrices)
    print(seed)
    print(subTree)
    print(appendedToNode)
    print(newTreeLK-oldTreeLK)
    print(newLK-originalLK)
    if abs((newTreeLK-oldTreeLK) - (newLK-originalLK))>1:
        print("diference is bigger than 1.0")



















