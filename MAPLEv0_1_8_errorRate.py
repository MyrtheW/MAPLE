#from MAPLEv0_1_8_experimental import *
def addErrorTerminalNode(node, errorRateOneThird):
    for entry in node.probVect: # loop over the lower likelihood entries in the genome list
        #for position specific errors; check if errorRate[i] / 3 > thresholdProb: in such case the error rate doesnt have to be considered for the position
        # for position specific errors; check if errorRate[i] > thresholdHighError: store as O vector.
        # if entry[0] <4: #Case CATGR
        #     entry.append(True)
        # You dont'have to to add an extra entry representing the flag to the genome list, because it can be checked whether a node is a child node
        if entry[0] == 6: # we are dealing with ambiguity characters other than N. We should not encounter this if onlyNambiguities==True:
            # derive the new probVector using the errorRate.
            sumUnnormalizedVector =  sum(bool(entry[-1]))
            if sumUnnormalizedVector ==2:
                for i in range4:             # M, instead of [0.5, 0.5, 0, 0] we will now get [0.5- ⅓ε, 0.5- ⅓ε,  ⅓ε,  ⅓ε]
                    if entry[-1][i]==0:
                        entry[-1][i] = errorRateOneThird
                    else: #if entry[-1][i]==0.5:
                        entry[-1][i] -= errorRateOneThird
            else: #sumUnnormalizedVector == 3:
                for i in range4:             # for V instead of [⅓, ⅓, ⅓, 0] we get, [ ⅓ - ε/9,  ⅓ -ε/9,  ⅓ -ε/9,  ⅓ε]
                    if entry[-1][i] == 0:
                        entry[-1][i] = errorRateOneThird
                    else:  # if entry[-1][i]==0.5:
                        entry[-1][i] -= errorRateOneThird / 3
            # this can be done faster with numpy
        #else: #entry[0]==5 of type N, we don't take this into account now.

def getFlags(entry1, entry2, node1isleaf, node2isleaf):
    if len(entry1) == 4 and entry1[0] < 5:
        flag1 = entry1[3]
    else:
        if node1isleaf == []:
            flag1 = True
        elif entry1[0] >= 5:  # cases of O or N
            flag1 = False
        else:  # other cases with 0 branch length:
            flag1 = False

    if len(entry2) >= 4 and entry2[0] < 5:  # or len(entry) >=3:
        flag2 = entry2[-1]
    else:
        if node2isleaf == []:
            flag2 = True
        elif entry2[0] >= 5:  # cases of O or N
            flag2 = False
        else:  # other cases with 0 branch length:
            flag2 = False
    assert (type(flag1) == bool)
    assert (type(flag2) == bool)
    return flag1, flag2

    # 2) add extra entry representing the flag to the genome list.
    # As a consequence, the entries with a flag will be of length 4 or 5 (if its branch length element bridges the root).
    # merge two lower child partial likelihood vectors to create a new one
    # (and also calculate the logLk of the merging if necessary, which is currently only useful for the root but could also be useful to calculate the overall total likelihood of the tree).
def mergeVectorsError(probVect1, bLen1, probVect2, bLen2, mutMatrix, returnLK=False, useRateVariation=False,
                 mutMatrices=None, node1isleaf=False, node2isleaf=False):
    indexEntry1, indexEntry2, pos = 0, 0, 0
    probVect = []
    cumulPartLk = 0.0
    entry1 = probVect1[indexEntry1]
    entry2 = probVect2[indexEntry2]
    end = min(entry1[1], entry2[1])
    while True:
        if entry1[0] == 5: # Case τ1 = N
            if entry2[0] == 5: # Case τ1 = τ2 = N
                pos = min(entry1[1], entry2[1])
                probVect.append((5, pos))
            elif entry2[0] < 5: #Case τ1 = N,  τ2 ϵ {A, C, T, G, R}
                pos = min(entry1[1], entry2[1])
                if len(entry2) == 2: # entry 2 is a child node or has bLength element of 0 with its predecessor.
                    if bLen2:
                        probVect.append((entry2[0], pos, bLen2, node2isleaf)) #flag = flag2 = node2isleaf in this case. If node 2 was not a leaf, it means it is a node with a 0 branch length element.
                    else: #0-branch length: no need to set the flag
                        probVect.append((entry2[0], pos))
                else:
                    assert (type(entry2[-1]) == bool)
                    assert (len(entry2) == 4)
                    if bLen2:
                        probVect.append((entry2[0], pos, entry2[2] + bLen2, entry2[3])) #set flag = flag2 (entry2[-1])
                    else:
                        probVect.append((entry2[0], pos, entry2[2], entry2[3])) #set flag = flag2 (entry2[3] = entry2[-1])
            else:  # case entry2 is "O" and entry1 is "N"; no flags needed.
                pos += 1
                if len(entry2) == 3:
                    if bLen2:
                        probVect.append((6, pos, bLen2, entry2[-1]))
                    else:
                        probVect.append((6, pos, entry2[-1]))
                else:
                    if bLen2:
                        probVect.append((6, pos, entry2[2] + bLen2, entry2[-1]))
                    else:
                        probVect.append((6, pos, entry2[2], entry2[-1]))
        elif entry2[0] == 5:  # entry2 is N
            if entry1[0] < 5: #Case τ2 = N,  τ1 ϵ {A, C, T, G, R}
                pos = min(entry1[1], entry2[1])
                if len(entry1) == 2:
                    if bLen1:
                        probVect.append((entry1[0], pos, bLen1, node1isleaf)) #flag = flag1 = node1isleaf in this case.
                    else: #0-branch length: no need to set the flag
                        probVect.append((entry1[0], pos))
                else:
                    assert (type(entry1[-1]) == bool)
                    assert (len(entry1) == 4)
                    if bLen1:
                        probVect.append((entry1[0], pos, entry1[2] + bLen1, entry1[3])) #set flag = flag1 (entry1[3])
                    else:
                        probVect.append((entry1[0], pos, entry1[2], entry1[3])) #set flag = flag1 (entry1[3])
            else:  # entry1 is "O" and entry2 is "N": no need for flags.
                pos += 1
                if len(entry1) == 3:
                    if bLen1:
                        probVect.append((6, pos, bLen1, entry1[-1]))
                    else:
                        probVect.append((6, pos, entry1[-1]))
                else:
                    if bLen1:
                        probVect.append((6, pos, entry1[2] + bLen1, entry1[-1]))
                    else:
                        probVect.append((6, pos, entry1[2], entry1[-1]))

        else:  # entry1 and entry2 are not "N"
            if entry1[0] < 5:
                if len(entry1) == 2:
                    totLen1 = bLen1
                else:
                    totLen1 = entry1[2]
                    if bLen1:
                        totLen1 += bLen1
            else:
                if len(entry1) == 3:
                    totLen1 = bLen1
                else:
                    totLen1 = entry1[2] #Question how can the entry be of type "N" or "O" and have a branch length element?
                    if bLen1:
                        totLen1 += bLen1

            if entry2[0] < 5:
                if len(entry2) == 2:
                    totLen2 = bLen2
                else:
                    totLen2 = entry2[2]
                    if bLen2:
                        totLen2 += bLen2
            else:
                if len(entry2) == 3:
                    totLen2 = bLen2
                else:
                    totLen2 = entry2[2]
                    if bLen2:
                        totLen2 += bLen2

            flag1, flag2 = getFlags(entry1, entry2, node1isleaf, node2isleaf)

            if entry2[0] == entry1[0] and entry2[0] < 5:  # entry1 and entry2 are two identical nucleotides #Case τ1 = τ2 ϵ {A, C, T, G, R}
                end = min(entry1[1], entry2[1])
                probVect.append((entry2[0], end)) #no flag is needed, it can be deduced that this is a 0  branch length entry
                if returnLK:
                    if entry2[0] == 4: # case entry2 is R
                        cumulPartLk += (totLen1 + totLen2) * (cumulativeRate[end] - cumulativeRate[pos])
                    else:
                        if useRateVariation:
                            cumulPartLk += mutMatrices[pos][entry1[0]][entry1[0]] * (totLen1 + totLen2)
                        else:
                            cumulPartLk += nonMutRates[entry1[0]] * (totLen1 + totLen2)
                    if flag1 + flag2:
                        # cumErrorRate  = cumulativeErrorRate[end]-cumulativeErrorRate[pos]  #cumulativeErrorRate should be log scaled already.
                        cumErrorRate = log(1 - errorRate) * (end - pos)  # OR simplified: cumErrorRate =-errorRate*(end-pos)
                        cumulPartLk -= cumErrorRate * (flag1 + flag2)
                pos = end
            elif (not totLen1) and (not totLen2) and entry1[0] < 5 and entry2[0] < 5:  # 0 distance between different nucleotides: merge is not possible
                if returnLK:
                    return None, float("-inf")
                else:
                    return None
            elif entry1[0] < 5:  # entry1 is a nucleotide
                if useRateVariation:
                    mutMatrix = mutMatrices[pos]
                if entry1[0] == 4:
                    i1 = refIndeces[pos]
                else:
                    i1 = entry1[0]
                # represent entry1 as an O entry

                def obtainPartialVec(i1, flag1, errorRate, method1=True):
                    if flag1:
                        partialVec1 = [flag1 * errorRate / 3] * 4  # without error rate [0.0, 0.0, 0.0, 0.0]
                        partialVec1[i1] = 1.0 - errorRate * flag1 #without error rate: 0
                        if totLen1:
                            mutatedPartialVec1 = []
                            if method1: # extensive calculation method 1
                                for i in range4:
                                    for j in range4:
                                        tot = 0.0
                                        for i in range4:
                                            tot += mutMatrix[j][i] * partialVec1[i]
                                        tot *= totLen1
                                        tot += partialVec1[-1][j]
                                        mutatedPartialVec1[j] = tot
                                partialVec1 = mutatedPartialVec1
                            else: # approximate calculation method 2
                                for i in range4:
                                    if i == i1:
                                        partialVec1[i1] *= (1.0 + mutMatrix[i][i] * totLen1)
                                    else:
                                        partialVec1[i1] *= (mutMatrix[i][i1] * totLen1)
                    else: # no flag 1
                        if totLen1:
                            partialVec1 = []
                            for i in range4:
                                if i == i1:
                                    partialVec1.append(1.0 + mutMatrix[i][i] * totLen1)
                                else:
                                    partialVec1.append(mutMatrix[i][i1] * totLen1)
                        else:
                            partialVec1 = [0.0, 0.0, 0.0, 0.0]
                            partialVec1[i1] = 1.0
                    return partialVec1

                newVec = obtainPartialVec(i1, flag1, errorRate)

                if entry2[0] == 6:  # entry1 is a nucleotide and entry2 is "O"
                    if totLen2:
                        for j in range4:
                            tot = 0.0
                            for i in range4:
                                tot += mutMatrix[j][i] * entry2[-1][i] # Question: Where is the chronicler delta here? shouldn't you add 1, in case i == j, as described for cases T1=O and T2=O
                            tot *= totLen2
                            tot += entry2[-1][j]
                            newVec[j] *= tot
                    else:
                        for j in range4:
                            newVec[j] *= entry2[-1][j]
                    sumV = sum(newVec)
                    if not sumV:
                        if returnLK:
                            return None, float("-inf")
                        else:
                            return None
                    for i in range4:
                        newVec[i] = newVec[i] / sumV
                    state = simplfy(newVec, refIndeces[pos])
                    pos += 1
                    if state == 6:
                        probVect.append((6, pos, newVec)) # no flag, because O type
                    else:
                        probVect.append((state, pos)) #no flag, because 0-branch length element
                    if returnLK:
                        cumulPartLk += log(sumV)
                else:  # entry1 and entry2 are nucleotides
                    if entry2[0] == 4:
                        i2 = refIndeces[pos]
                    else:
                        i2 = entry2[0]
                    #TODO: continue here;
                    partialVec2 = obtainPartialVec(i2, flag2, errorRate)
                    if totLen2:
                        for i in range4:
                            newVec[i] *= partialVec2
                            # if i == i2:
                            #     newVec[i] *= 1.0 + mutMatrix[i][i] * totLen2
                            # else:
                            #     newVec[i] *= mutMatrix[i][i2] * totLen2
                        sumV = sum(newVec)
                        for i in range4: #normalize
                            newVec[i] = newVec[i] / sumV
                        state = simplfy(newVec, refIndeces[pos])
                        pos += 1
                        if state == 6:
                            probVect.append((6, pos, newVec))
                        else:
                            probVect.append((state, pos)) #no flag, because 0-branch length element
                        # probVect.append((6,pos,newVec))
                        if returnLK:
                            cumulPartLk += log(sumV)
                    else:
                        pos += 1
                        probVect.append((entry2[0], pos))
                        if returnLK:
                            cumulPartLk += log(newVec[i2])

            else:  # entry1[0]==6:
                if useRateVariation:
                    mutMatrix = mutMatrices[pos]
                if totLen1:
                    newVec = []
                    for i in range4:
                        tot = 0.0
                        for j in range4:
                            tot += mutMatrix[i][j] * entry1[-1][j]
                        tot *= totLen1
                        tot += entry1[-1][i]
                        newVec.append(tot)
                else:
                    newVec = list(entry1[-1])

                if entry2[0] == 6: # Case τ1= τ2= 'O'
                    if totLen2:
                        for i in range4:
                            tot = 0.0
                            for j in range4:
                                tot += mutMatrix[i][j] * entry2[-1][j]
                            tot *= totLen2
                            tot += entry2[-1][i] #is this the cronicler delta part ?
                            newVec[i] *= tot
                    else:
                        for i in range4:
                            newVec[i] *= entry2[-1][i]
                    sumV = sum(newVec)
                    if not sumV:
                        if returnLK:
                            return None, float("-inf")
                        else:
                            return None
                    for i in range4:
                        newVec[i] = newVec[i] / sumV
                    state = simplfy(newVec, refIndeces[pos])
                    pos += 1
                    if state == 6:
                        probVect.append((6, pos, newVec))
                    else:
                        probVect.append((state, pos))
                    if returnLK:
                        cumulPartLk += log(sumV)
                else:  # entry2 is a nucleotide and entry1 is "O"
                    if entry2[0] == 4:
                        i2 = refIndeces[pos]
                    else:
                        i2 = entry2[0]
                    if totLen2:
                        partialVec2 = obtainPartialVec(i2, flag2, errorRate)
                        for i in range4:
                            newVec[i] *= partialVec2
                    # if totLen2:
                    #     for i in range4:
                    #         if i == i2:
                    #             newVec[i] *= (1.0 + mutMatrix[i][i] * totLen2)
                    #         else:
                    #             newVec[i] *= mutMatrix[i][i2] * totLen2
                        sumV = sum(newVec)
                        for i in range4:
                            newVec[i] = newVec[i] / sumV
                        state = simplfy(newVec, refIndeces[pos])
                        pos += 1
                        if state == 6:
                            probVect.append((6, pos, newVec))
                        else:
                            probVect.append((state, pos))
                        if returnLK:
                            cumulPartLk += log(sumV)
                    else:
                        if not newVec[i2]:
                            if returnLK:
                                return None, float("-inf")
                            else:
                                return None
                        pos += 1
                        probVect.append((entry2[0], pos))
                        if returnLK:
                            cumulPartLk += log(newVec[i2])

        if pos == lRef:
            break
        if pos == entry1[1]:
            indexEntry1 += 1
            entry1 = probVect1[indexEntry1]
        if pos == entry2[1]:
            indexEntry2 += 1
            entry2 = probVect2[indexEntry2]

    if verbose:
        print("Merged vector at root ")
        print(probVect)
    # check if the final  probVect can be simplified by merging consecutive entries
    shorten(probVect)
    if verbose:
        print("Shortened root merging ")
        print(probVect)
    if returnLK:
        return probVect, cumulPartLk
    else:
        return probVect

# def probVectTerminalNode(diffs):
#     if diffs is None:
#         probVect=[(5,lRef)]
#         return probVect
#     pos=1
#     probVect=[]
#     for m in diffs:
#             currPos=m[1]
#             if currPos>pos:
#                 #region where the node with branch length bLen is identical to the ref.
#                 if errorRate / 3 < thresholdProb:
#                     probVect.append((4,currPos-1))
#                 else:
#                     probVect.append((6, currPos-1, errLikelihoodMatrix[allelesListLow[refIndeces[i]]])) #Myrthe
#             pos=currPos
#             if m[0]=="n" or m[0]=="-":
#                 if len(m)>2:
#                     length=m[2]
#                 else:
#                     length=1
#                 #region with no info, store last position and length.
#                 probVect.append((5,currPos+length-1))
#                 pos=currPos+length
#             elif m[0] in allelesListLow:
#                 #position at which node allele is sure but is different from the reference.
#                 if errorRate/3 < thresholdProb: #Myrthe #Only for uniform error rate.
#                     probVect.append((allelesLow[m[0]],currPos))
#                 else:
#                     probVect.append((6,currPos,errLikelihoodMatrix[m[0]]))
#                 pos=currPos+1
#             else:
#                 # non-"n" ambiguity character; for now interpret this as ambiguity instead of as a polymorphism.
#                 if onlyNambiguities:
#                     # if user asks to, to make things easier, interpret any ambiguity as an "n".
#                     probVect.append((5,currPos))
#                 else:
#                     #otherwise, store as an "other" scenario, where each nucleotide has its own partial likelihood.
#                     probVect.append((6,currPos,errLikelihoodMatrix[m[0]])) #Myrthe
#                 pos=currPos+1
#     if pos<=lRef:
#         probVect.append((4,lRef))
#     return probVect



def reCalculateWithErrors(root,mutMatrix, errorRate, checkExistingAreCorrect=False,countNodes=False,
                          countPseudoCounts=False,pseudoMutCounts=None,
                          data=None,useRateVariation=False,mutMatrices=None):

    # this function is a copy of recalculate all genome lists. Similarly to that, we loop in post-order, starting with calculating the lower-likelihoods of the leaves.
    # after calculating all lowers, we will update other genome lists.
    # this is also the point where we will first add the error rates.
    node = root # direction 0 means from parent, direction 1 means from a child
    lastNode = None
    direction = 0
    errorRateOneThird = errorRate / 3
    while node != None:
        if direction == 0:
            if node.children:
                node = node.children[0]
            else: #leaf node:
                addErrorTerminalNode(node, errorRateOneThird)
                if data != None: #if the initial tree is formed, data should be an empty list or None, since all nodes have been added.
                    raise Exception("Error: Data should be none ")
                lastNode = node
                node = node.up
                direction = 1
        else:
            if lastNode == node.children[0]: # arrived here from the left child node.
                node = node.children[1] # visit the right child node first.
                direction = 0
            else: # arrived here from the left child node, so both child nodes have been visited.
                newLower = mergeVectorsError(node.children[0].probVect, node.children[0].dist, node.children[1].probVect,
                                        node.children[1].dist, mutMatrix, returnLK=False,
                                        useRateVariation=useRateVariation, mutMatrices=mutMatrices)
                # if checkExistingAreCorrect:
                #     if areVectorsDifferentDebugging(newLower, node.probVect):
                #         print(
                #             "Inside reCalculateAllGenomeLists(), new lower at node is different from the old one, and it shouldn't be.")
                #         print(newLower)
                #         print(node.probVect)
                #         exit()
                if newLower == None:
                    if not node.children[0].dist:
                        nodeList = []
                        updateBLen(nodeList, node, mutMatrix, useRateVariation=useRateVariation,mutMatrices=mutMatrices) #TODO: write a version for error Rates.
                        updatePartials(nodeList, mutMatrix, useRateVariation=useRateVariation, mutMatrices=mutMatrices)
                    elif not node.children[1].dist:
                        nodeList = []
                        updateBLen(nodeList, node.children[1], mutMatrix, useRateVariation=useRateVariation, mutMatrices=mutMatrices)
                        updatePartials(nodeList, mutMatrix, useRateVariation=useRateVariation, mutMatrices=mutMatrices)
                    else:
                        raise Exception( "Strange, distances>0 but inconsistent lower genome list creation in reCalculateAllGenomeLists()")
                else:
                    node.probVect = newLower
                lastNode = node
                node = node.up
                direction = 1

    #now update the other genome lists for the root
    node=root
    if node.children:
        newVect=rootVector(node.children[1].probVect,node.children[1].dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
        if checkExistingAreCorrect:
            if areVectorsDifferentDebugging(newVect,node.probVectUpRight):
                print("new probVectUpRight at root is different from the old one, and it shouldn't be.")
                print(newVect)
                print(node.probVectUpRight)
                exit()
        node.probVectUpRight=newVect
        newVect=rootVector(node.children[0].probVect,node.children[0].dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
        if checkExistingAreCorrect:
            if areVectorsDifferentDebugging(newVect,node.probVectUpLeft):
                print("new probVectUpLeft at root is different from the old one, and it shouldn't be; distance of child node "+str(node.children[0].dist))
                print(node.dist)
                print(node.children)
                print(len(node.children))
                print(node.children[0].dist)
                print(node.children[1].dist)
                print(node.children[0].probVect)
                print(createBinaryNewick(node))
                print(newVect)
                print(node.probVectUpLeft)
                exit()
        node.probVectUpLeft=newVect

        #now traverse the tree downward and update the non-lower genome lists for all other nodes of the tree.
        lastNode=None
        node=node.children[0]
        direction=0
        while node!=None:
            if direction==0:
                if node==node.up.children[0]:
                    vectUp=node.up.probVectUpRight
                else:
                    vectUp=node.up.probVectUpLeft
                if node.dist:
                    #newVect=getTot(node,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                    if countPseudoCounts:
                        updatePesudoCounts(vectUp,node.probVect,pseudoMutCounts)
                    newVect=mergeVectorsUpDown(vectUp,node.dist/2,node.probVect,node.dist/2,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                    if checkExistingAreCorrect:
                        if areVectorsDifferentDebugging(newVect,node.probVectTotUp):
                            print("new probVectTotUp at node is different from the old one, and it shouldn't be.")
                            print(newVect)
                            print(node.probVectTotUp)
                            exit()
                    node.probVectTotUp=newVect
                    # if node.dist>=2*minBLenForMidNode:
                    # 	createFurtherMidNodes(node,vectUp,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                if node.children:
                    newUpRight=mergeVectorsUpDown(vectUp,node.dist,node.children[1].probVect,node.children[1].dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                    if newUpRight==None:
                        if (not node.children[1].dist):
                            nodeList=[]
                            updateBLen(nodeList,node.children[1],mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                            updatePartials(nodeList,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                        elif (not node.dist):
                            nodeList=[]
                            updateBLen(nodeList,node,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                            updatePartials(nodeList,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                        else:
                            print("Strange, distances>0 but inconsistent upRight genome list creation in reCalculateAllGenomeLists()")
                            exit()
                    else:
                        if checkExistingAreCorrect:
                            if areVectorsDifferentDebugging(newUpRight,node.probVectUpRight):
                                print("new probVectUpRight at node is different from the old one, and it shouldn't be.")
                                print(node.children)
                                print(node.dist)
                                print(node.children[1].dist)
                                print()
                                print(newUpRight)
                                print()
                                print(node.probVectUpRight)
                                print()
                                print(vectUp)
                                print()
                                print(node.children[1].probVect)
                                exit()
                        node.probVectUpRight=newUpRight
                    newUpLeft=mergeVectorsUpDown(vectUp,node.dist,node.children[0].probVect,node.children[0].dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                    if newUpLeft==None:
                        if(not node.children[0].dist):
                            nodeList=[]
                            updateBLen(nodeList,node.children[0],mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                            updatePartials(nodeList,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                        elif (not node.dist):
                            nodeList=[]
                            updateBLen(nodeList,node,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                            updatePartials(nodeList,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                        else:
                            print("Strange, distances>0 but inconsistent upLeft genome list creation in reCalculateAllGenomeLists()")
                            exit()
                    else:
                        if checkExistingAreCorrect:
                            if areVectorsDifferentDebugging(newUpLeft,node.probVectUpLeft):
                                print("new probVectUpLeft at node is different from the old one, and it shouldn't be.")
                                print(node.children)
                                print(node.dist)
                                print(node.children[0].dist)
                                print(newUpLeft)
                                print()
                                print(node.probVectUpLeft)
                                print()
                                print(vectUp)
                                print()
                                print(node.children[0].probVect)
                                exit()
                        node.probVectUpLeft=newUpLeft
                    node=node.children[0]
                else:
                    lastNode=node
                    node=node.up
                    direction=1
            else:
                if lastNode==node.children[0]:
                    node=node.children[1]
                    direction=0
                else:
                    lastNode=node
                    node=node.up
                    direction=1


