import sys
from sys import exit
from math import log
import argparse
from time import time
import os.path

#©EMBL-European Bioinformatics Institute, 2021
# MAPLE code to estimate a tree by maximum likelihood from a MAPLE format input.

#things that it makes sense to try soon:
#TODO use analytical branch length optimization to improve placement and topological moves?

#things in progress
#TODO C++ implementation
#TODO Parallelization

#longer-term plans
#TODO Model sequencing errors
#TODO codon models, selection
#TODO indels, alignment
#TODO Timed trees
#TODO Bayesian implementation in BEAST
#TODO Recombination
#TODO phylogeography
#TODO phylodynamics/selection

#things one could MAYBE also try in the future:
#TODO use slots (https://www.geeksforgeeks.org/slots-in-python/) to reduce memory cost?
#TODO develop more nuanced exploration of topology space: reattempt to re-place only subtrees whose root have been affected in the last round of topology search (which is now already done),
# and their neighbours (which is not implemented yet)?
#TODO for increased accuracy, one could calculate probVectTot genome lists in a way that would not depend by the choice of which two genome lists to merge first.
#For now implemented consistent calculation of tot likelihoods so to reduce inconsinstencies met.
#TODO would it be convenient to remove the genome position element from entries of type ACGTO ? Sounds like a lot of changes, and might make things slower and less readable?
#TODO if the model is reversible, no need to consider complicated case of root position, so the genome vector lists can be simpler and I can avoid additional calculations?
#TODO does it make sense to have createFurtherMidNodes() or maybe rather save on the memory and create these genome lists anew each time? Maybe rather find optimal values for the threshold.
#TODO try to replace appendProb() with appendProbNode()? it might be slightly slower, but it would be one fewer function in the code.
#TODO do multiple rounds of SPR search with increasing thresholds? Now performing a quick short-range SPR tree traverse before the proper one.
#TODO create an initial stage for the placement with very low thresholds; then, a second placement stage would use the default thresholds and work
# like the SPR search but starting the placement search from the node found by the first stage; finally, the last stage would be as usual.
# I tried this but it comes at little computational advantage and at large accuracy cost.


parser = argparse.ArgumentParser(description='Estimate a tree from a diff format and using iterative approximate maximum likelihood sample placement.')
parser.add_argument('--input',default="/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/phylogenetic_inference/2021-03-31_unmasked_differences_reduced.txt_consensus-based.txt", help='input MAPLE file name; should contain first the reference genome and then the difference of all samples with respet to the reference.')
parser.add_argument('--reference',default="", help='optional input reference file name. By default it assumes instead that the reference is part of the MAPLE format input.')
parser.add_argument('--output',default="/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/phylogenetic_inference/MAPLE", help='output path and identifier to be used for newick output file.')
parser.add_argument('--inputTree',default="", help='input newick tree file name; this is optional, and is used for online inference (or for Robinson-Foulds distance calculation if option --inputRFtrees is also used).')
parser.add_argument("--largeUpdate", help="When using option --inputTree to do online inference, by default the tree is only updated locally where the new sequences are inserted. Use this option instead to perform a thorough update of the phylogeny.", action="store_true")
parser.add_argument('--inputRFtrees',default="", help='file name with input newick trees to be compared with the input tree in option --inputTree to calculate Robinson-Foulds distances; this option will turn off the normal MAPLE estimation - only RF distances will be calculated. Each newick tree contained in this file will be compared to the single tree specified with option --inputTree .')
parser.add_argument("--onlyNambiguities", help="Treat all ambiguities as N (total missing information).", action="store_true")
parser.add_argument("--thresholdProb",help="relative probability threshold used to ignore possible states with very low probabilities.",  type=float, default=0.00000001)
parser.add_argument("--thresholdLogLK",help="logLK difference threshold to consider a logLk close to optimal.",  type=float, default=200.0)
parser.add_argument("--thresholdLogLKtopology",help="logLK difference threshold to consider a logLk close to optimal when looking for topology improvements.",  type=float, default=160.0)
parser.add_argument("--allowedFails",help="Number of times one can go down the tree without inclreasing placement likelihood before the tree traversal is stopped (only applies to non-0 branch lengths).",  type=int, default=5)
parser.add_argument("--allowedFailsTopology",help="Number of times one can crawl along the tree decreasing placement likelihood before the tree traversal is stopped during topology search (only applies to non-0 branch lengths).",  type=int, default=4)
#parser.add_argument("--bLenAdjustment",help="If >1, try placing also with a longer bLen than the standard bLen.",  type=int, default=1)
parser.add_argument("--verbose", help="Print to screen a lot of stuff.", action="store_true")
parser.add_argument("--debugging", help="Test that likelihoods are calculated and updated as expected - time consuming and only meant for small trees for debugging purposes.", action="store_true")
parser.add_argument("--model", help="Which substitution model should be used. Allowed models so far are JC, GTR (default) or UNREST.", default="GTR")
parser.add_argument("--overwrite", help="Overwrite previous results if already present.", action="store_true")
parser.add_argument("--nonBinaryTree", help="Write output tree with multifurcations - by default the tree is written as binary so to avoid problems reading the tree in other software.", action="store_true")
parser.add_argument("--numTopologyImprovements",help="Number of times we traverse the tree looking for topological improvements.",  type=int, default=1)
parser.add_argument("--thresholdTopologyPlacement",help="Don't try to re-place nodes that have current appending logLK cost above this threshold.",  type=float, default=-0.01)
#parser.add_argument("--minBLenForMidNode",help="Don't try to place samples mid-branch if the considered branch is below this length in number of mutations genome-wide.",  type=float, default=4.1)
parser.add_argument("--updateSubstMatrixEveryThisSamples",help="How many new samples to place before each update of the substitution rate matrix.",  type=int, default=25)
parser.add_argument("--nonStrictInitialStopRules", help="If specified, then during the initial placement stage, a slower, non-strict rule for stopping the placement search is applied: the search is stopped if enough many consencutive LK worsening are observed, AND if LK is below the considered threshold.", action="store_true")
parser.add_argument("--strictTopologyStopRules", help="If specified, then during the topological improvement stage, a slower, non-strict rule for stopping the SPR search is applied: the search is stopped if enough many consencutive LK worsening are observed, AND if LK is below the considered threshold.", action="store_true")
parser.add_argument("--thresholdDiffForUpdate",help="Consider the probability of a new partial changed if the difference between old and new if above this threshold.",  type=float, default=0.0000001)
parser.add_argument("--thresholdFoldChangeUpdate",help="Consider the probability of a new partial changed, if the fold difference between old and new if above this threshold.",  type=float, default=1.001)
#parser.add_argument("--minBLen",help="Don't attempt branch lengths below this number of mutations (genome-wide): just consider 0 length instead.",  type=float, default=0.2)
#parser.add_argument("--maxBLen",help="Don't attempt branch lengths above this number of mutations (genome-wide).",  type=float, default=40.0)
parser.add_argument("--thresholdLogLKconsecutivePlacement",help="logLK difference threshold to consider something as a significant decrease in log-LK when considering consecutive likelihood decreases.",  type=float, default=0.01)
parser.add_argument("--thresholdLogLKwholeTopologyImprovement",help="logLK difference threshold to consider something as a significant decrease in log-LK when considering whole rounds of SPR moves.",  type=float, default=1.0)
parser.add_argument("--calculateLKfinalTree", help="Calculate the log-LK of the final infered tree.", action="store_true")
parser.add_argument("--fast", help="Set parameters so to run tree inference faster; this will be less accurate in cases of high complexity, for example with recombination, sequencing errors, etc. It will overrule user choices for options --thresholdLogLK , --thresholdLogLKtopology , --allowedFails , --allowedFailsTopology .", action="store_true")
#parser.add_argument("--fastPreliminaryPlacement", help="When placing each sample run first a fast search and then start the proper search from there. Not recommended as this comes at accuracy cost and brings little to no computational advantage.", action="store_true")
parser.add_argument("--noFastTopologyInitialSearch", help="Don't run a fast short-range topology search before the extensive one.", action="store_true")
parser.add_argument("--noOptimizeBranchLengths", help="Don't run a final round of detailed branch length optimization.", action="store_true")
parser.add_argument("--rateVariation", help="Estimate and use rate variation: the model assumes one rate per site, and the rates are assumed independently (no rate categories). This might cause overfitting if the dataset is not large enough, but in any case one would probably only use MAPLE for large enough datasets.", action="store_true")
parser.add_argument("--minBLenSensitivity",help="Fraction of a mutation to be considered as a precision for branch length estimation (default 0.001, which means branch lengths estimated up to a 1000th of a mutation precision).",  type=float, default=0.001)
parser.add_argument("--factorOptimizePlacementLKvsSearchLK",help="Parameter that determines for how many nodes visited to perform branch length optimization before performing SPR etc. Default 0.04 means we use a log-LK threshold 25 times smaller than in the default tree traversal.",  type=float, default=0.04)

# Error rate additional parameters:
parser.add_argument("--errorRate", help="", type=float, default=0, help="Provide an error rate to be used during the inference")
parser.add_argument("--errorRateSiteSpecific", help="provide a file directory to a file that contains the siteSpecific error rates", type=str, default=None) #errorRateSiteSpecific
parser.add_argument("--genomeLength", help="", type=float, default=0, help="When you wish to calculate the RF(L) distances after your tree(s) have been inferred, provide the genomelength used to scale the branch lengths in the newick file of the true tree.")
parser.add_argument("--benchmarkingFile", help="", type=str, default=None, help="If you want the results to be stored in a comparable format, provide a tsv file and the path to it, e.g. 'results/benchmarkingFile.tsv'. For this analysis you should also provide the true tree. ")
parser.add_argument("--trueTree", help="", type=str, default=None, help="Provide the path and filename to newick file of the true tree that should be used for calculating the RF and RFL distance")

args = parser.parse_args()

onlyNambiguities=args.onlyNambiguities
thresholdProb=args.thresholdProb
verbose=args.verbose
debugging=args.debugging
inputFile=args.input
outputFile=args.output
refFile=args.reference
allowedFails=args.allowedFails
allowedFailsTopology=args.allowedFailsTopology
model=args.model
#bLenAdjustment=args.bLenAdjustment
#if bLenAdjustment>1.0+thresholdProb:
#	tryOtherBLen=True
#else:
#	tryOtherBLen=False
thresholdLogLK=args.thresholdLogLK
thresholdLogLKtopology=args.thresholdLogLKtopology
overwrite=args.overwrite
binaryTree=(not args.nonBinaryTree)
numTopologyImprovements=args.numTopologyImprovements
thresholdTopologyPlacement=args.thresholdTopologyPlacement
#minBLenForMidNode=args.minBLenForMidNode
updateSubstMatrixEveryThisSamples=args.updateSubstMatrixEveryThisSamples
strictInitialStopRules=(not args.nonStrictInitialStopRules)
strictTopologyStopRules=args.strictTopologyStopRules
thresholdDiffForUpdate=args.thresholdDiffForUpdate
thresholdFoldChangeUpdate=args.thresholdFoldChangeUpdate
#minBLen=args.minBLen
#maxBLen=args.maxBLen
thresholdLogLKconsecutivePlacement=args.thresholdLogLKconsecutivePlacement
thresholdLogLKwholeTopologyImprovement=args.thresholdLogLKwholeTopologyImprovement
calculateLKfinalTree=args.calculateLKfinalTree
#fastPreliminaryPlacement=args.fastPreliminaryPlacement
fastTopologyInitialSearch=(not args.noFastTopologyInitialSearch)
optimizeBranchLengths=(not args.noOptimizeBranchLengths)
example=False
runFluTest=False
runFast=args.fast
if runFast:
	thresholdLogLK=160.0
	allowedFails=4
	allowedFailsTopology=2
	thresholdLogLKtopology=80.0
	thresholdTopologyPlacement=-1.0

strictInitialStopRulesInitial=True
allowedFailsInitial=3
thresholdLogLKinitial=120.0
strictTopologyStopRulesInitial=True
allowedFailsTopologyInitial=1
thresholdLogLKtopologyInitial=40.0
thresholdTopologyPlacementInitial=-1.0

inputTree=args.inputTree
inputRFtrees=args.inputRFtrees
largeUpdate=args.largeUpdate
rateVariation=args.rateVariation
genomeLength = args.genomeLength
benchmarkingFile=args.benchmarkingFile
trueTree=args.trueTree
errorRate = args.errorRate
errorRateSiteSpecific = args.errorRateSiteSpecific
#When site specific error rates are not used in the end, simply comment out all 12 linescontainging  'errorRateSiteSpecific'

#ratio of likelihood cost tree search threshold vs likelihood threshold for placement optimization search
factorOptimizePlacementLKvsSearchLK=args.factorOptimizePlacementLKvsSearchLK
minBLenSensitivity=args.minBLenSensitivity


#Class defininf nodes of the tree
class Tree(object):
    def __init__(self, name='', children=None, dist=1.0):
        if name!='':
            self.name = name
        self.dist = dist
        self.children = []
        self.up=None
        self.dirty=True
        if children is not None:
            for child in children:
                self.add_child(child)
    def __repr__(self):
        try:
            return str(self.name)
        except AttributeError:
            try:
                return self.name
            except AttributeError:
                return ""
    def add_child(self, node):
        assert isinstance(node, Tree)
        self.children.append(node)


#function to read input newick string
def readNewick(nwFile,multipleTrees=False,dirtiness=True, divideBranchLengthsBy=1):
    phyloFile=open(nwFile)
    trees=[]
    line=phyloFile.readline()
    while line!="":
        while line=="\n":
            line=phyloFile.readline()
        if line=="":
            break
        nwString=line.replace("\n","")

        index=0
        node=Tree()
        node.dirty=dirtiness
        name=""
        distStr=""
        finished=False
        while index<len(nwString):
            if nwString[index]=="(":
                newNode=Tree()
                newNode.minorSequences=[]
                newNode.dirty=dirtiness
                node.add_child(newNode)
                newNode.up=node
                node=newNode
                index+=1
            elif nwString[index]==";":
                trees.append(node)
                finished=True
                break
                #return node
            elif nwString[index]=="[":
                while nwString[index]!="]":
                    index+=1
                index+=1
            elif nwString[index]==":":
                index+=1
                while nwString[index]!="," and nwString[index]!=")" and nwString[index]!=";":
                    distStr+=nwString[index]
                    index+=1
            elif nwString[index]==",":
                if name!="":
                    node.name=name
                    name=""
                if distStr!="":
                    node.dist=float(distStr)/divideBranchLengthsBy
                    distStr=""
                newNode=Tree()
                newNode.minorSequences=[]
                newNode.dirty=dirtiness
                node=node.up
                node.add_child(newNode)
                newNode.up=node
                node=newNode
                index+=1
            elif nwString[index]==")":
                if name!="":
                    node.name=name
                    name=""
                if distStr!="":
                    node.dist=float(distStr)/divideBranchLengthsBy
                    distStr=""
                node=node.up
                index+=1
            else:
                name+=nwString[index]
                index+=1
        if not finished:
            print("Error, final character ; not found in newick string in file "+nwFile+".")
            raise Exception("exit")

        if not multipleTrees:
            break
        line=phyloFile.readline()

    phyloFile.close()
    return trees


#function that changes multifurcating tree structure into a binary tree by adding 0-length branches/nodes
def makeTreeBinary(root):
    nodesToVisit=[root]
    while nodesToVisit:
        node=nodesToVisit.pop()
        if node.children:
            while len(node.children)>2:
                child2=node.children.pop()
                child1=node.children.pop()
                newParent=Tree(dist=False)
                newParent.add_child(child1)
                newParent.add_child(child2)
                child1.up=newParent
                child2.up=newParent
                newParent.up=node
                node.children.append(newParent)
            nodesToVisit.append(node.children[0])
            nodesToVisit.append(node.children[1])




def prepareTreeComparison(t1,rooted=False,minimumBLen=0.000006, addRootRFL = False):
    """
    #Robinson-Foulds distance (1981) using a simplification of the algorithm from Day 1985.
    #this function prepare the data to compare trees to a reference one t1.
    #I split in two functions (the second one is RobinsonFouldsWithDay1985() below) so that I don't have to repeat these steps for the reference tree if I compare to the same reference tree multiple times.

    I extended this to calculate the RFL distance (see the Robinson Foulds paper) where they describe this.
    Basically, for every branch that is shared among the two trees, add the absolute difference of the two branches to the RFL.
    (This is first stored seperately as the KL value first).
    For every branch that is only in either of the two trees, add the length to the RFL distance.
    For the tree that is 'prepared for comparison', usually the true tree, I will store two hash tables (python dictionaries) together containing all branch lengths.
    One contains all leaf-branches. These are indexed by the leaf's name, usually a number (node.name).
    The other dictionary contains all inner branches, indexed by the <Left, Right> values of the subtree that is attached by the branch of interest.
    (One could improve the implementation by writing one's own hash maps and functions. Given that the indices have certain limits to what values they can become, it would be easy to obtain a perfect hash function)
    Furthermore, I calculate the sum of all the inner branch lengths. Whenever in the second function of the Day's algorithm (the comparison with the estimated tree),
    a branch is found in both trees, I add the absolute difference between the two branches, and substract the branch length of the true tree.
    I don't do this for leaf branches, as we can be sure that a leaf branch exists in both trees (if higher than the minimum length).
    (I do this because it is less straight forward to add branch lengths of the true tree when they are not found during the comparison with the estimated tree)

    Similar to the normal RF distance function, I do not take into account branches with length smaller than the minimum branch length.
    However, cases may occur that a branch exists in both trees, but in one of them it is lower than the minimum branch length, making the RFL less accurate.

    :param t1: input tree
    :param rooted:
    :param minimumBLen:
    :param addRootRFL: Should the distance from the root node up be taken into account? I think this has usually no significance, and is always set to 1.0.
    :return: Various metrics, among which the RFL.
    """
    #dictionary of values given to sequence names
    leafNameDict={}
    #list of sequence names sorted according to value
    leafNameDictReverse=[]
    #table containing clusters in the tree
    nodeTable=[]
    #dictionaries with branch lengths
    branchLengthDict, leafDistDict = {}, {}
    sumBranchLengths = 0 # sum of all branch lengths, except from the leaves
    #if comparing as unrooted trees, calculate tot num of leaves (using a postorder traversal), which will become useful later
    if not rooted:
        nLeaves=0
        node=t1
        movingFrom=0
        while node!=t1.up:
            if movingFrom==0: #0 means reaching node from parent, 1 means coming back from a child
                if len(node.children)==0:
                    nLeaves+=1
                    nextNode=node.up
                    movingFrom=1
                    nodeTable.append([0,0])
                else:
                    nextNode=node.children[0]
                    movingFrom=0
                    node.exploredChildren=0
            else:
                nChildren=len(node.children)
                node.exploredChildren+=1
                if node.exploredChildren==nChildren:
                    nextNode=node.up
                    movingFrom=1
                else:
                    nextNode=node.children[node.exploredChildren]
                    movingFrom=0
            node=nextNode

    #implementing a non-recursive postorder traversal to assign values to internal nodes to fill nodeTable
    leafCount=0
    node=t1
    movingFrom=0
    lastL=float("inf")
    lastR=float("-inf")
    lastDesc=0
    numBranches=0
    while node!=t1.up:
        if movingFrom==0: #0 means reaching node from parent, 1 means coming back from a child
            if len(node.children)==0:
                node.name=(node.name).replace("?","_").replace("&","_")
                leafNameDict[node.name]=leafCount
                leafNameDictReverse.append(node.name)
                if rooted:
                    nodeTable.append([0,0])
                lastL=leafCount
                lastR=leafCount
                lastDesc=1
                leafCount+=1
                nextNode=node.up
                movingFrom=1
                leafDistDict[node.name] = node.dist
            else:
                node.exploredChildren=0
                node.maxSoFar=float("-inf")
                node.minSoFar=float("inf")
                node.nDescendants=0
                nextNode=node.children[0]
                movingFrom=0
        else:
            nChildren=len(node.children)
            node.exploredChildren+=1
            if lastL<node.minSoFar:
                node.minSoFar=lastL
            if lastR>node.maxSoFar:
                node.maxSoFar=lastR
            node.nDescendants+=lastDesc
            if node.exploredChildren==nChildren:
                sumBranchLengths += node.dist
                nextNode=node.up
                movingFrom=1
                lastL=node.minSoFar
                lastR=node.maxSoFar
                lastDesc=node.nDescendants
                if node==t1:
                    nodeTable[lastR][0]=lastL
                    nodeTable[lastR][1]=lastR
                    branchLengthDict[(lastL, lastR)] = node.dist  # , node.nDescendants)]
                    if not rooted or not addRootRFL: sumBranchLengths-=node.dist
                else:
                    if node.dist>minimumBLen:
                        numBranches+=1
                        if rooted or lastL>0:
                            if node==node.up.children[-1]:
                                nodeTable[lastL][0]=lastL
                                nodeTable[lastL][1]=lastR
                                branchLengthDict[(lastL, lastR)] = node.dist #, node.nDescendants)]
                            else:
                                nodeTable[lastR][0]=lastL
                                nodeTable[lastR][1]=lastR
                                branchLengthDict[(lastL, lastR)] = node.dist
                        else: # re-root at leaf 0, so flip the values for the current branch if it contains leaf 0.
                                flippedL=lastR+1
                                flippedR=nLeaves-1
                                nodeTable[flippedL][0]=flippedL
                                nodeTable[flippedL][1]=flippedR
                                branchLengthDict[(flippedL, flippedR)] = node.dist

            else:
                nextNode=node.children[node.exploredChildren]
                movingFrom=0
        node=nextNode
    return leafNameDict, nodeTable, leafCount, numBranches, leafDistDict, branchLengthDict, sumBranchLengths

#Robinson-Foulds distance (1981) using a simplification of the algorithm from Day 1985.
#this function compares the current tree t2 to a previous one for which prepareTreeComparison() was run.
# Example usage: leafNameDict, nodeTable, leafCount, numBranches = prepareTreeComparison(phyloTrue,rooted=False)
# numDiffs, normalisedRF, leafCount, foundBranches, missedBranches, notFoundBranches = RobinsonFouldsWithDay1985(phyloEstimated,leafNameDict,nodeTable,leafCount,numBranches,rooted=False)
def RobinsonFouldsWithDay1985(t2,leafNameDict,nodeTable,leafCount,numBranches,leafDistDict, branchLengthDict, sumBranchLengths, rooted=False,minimumBLen=0.000006, addRootRFL = False):
    #implementing a non-recursive postorder traversal to check branch existance in the reference tree
    node=t2
    #branches in reference tree that are also in t2
    foundBranches=0
    #branches in t2 that are not found in the reference
    missedBranches=0
    movingFrom=0
    lastL=float("inf")
    lastR=float("-inf")
    lastDesc=0
    visitedLeaves=0
    RFL = sumBranchLengths
    KF = 0
    while node!=t2.up:
        if movingFrom==0: #0 means reaching node from parent, 1 means coming back from a child
            if len(node.children)==0:
                node.name=(node.name).replace("?","_").replace("&","_")
                if node.name in leafNameDict:
                    leafNum=leafNameDict[node.name]

                else:
                    print(node.name+" not in reference tree - aborting RF distance")
                    return None, None, None, None, None, None
                lastL=leafNum
                lastR=leafNum
                lastDesc=1
                nextNode=node.up
                movingFrom=1
                visitedLeaves+=1
                trueDist = leafDistDict[node.name]
                KF += abs(trueDist- node.dist)   #As described, I have not added branch lengths of leaf nodes to sumBranchLengths, so no need for: RFL=- trueDist
            else:
                node.exploredChildren=0
                node.maxSoFar=float("-inf")
                node.minSoFar=float("inf")
                node.nDescendants=0
                nextNode=node.children[0]
                movingFrom=0
        else:
            nChildren=len(node.children)
            node.exploredChildren+=1
            if lastL<node.minSoFar:
                node.minSoFar=lastL
            if lastR>node.maxSoFar:
                node.maxSoFar=lastR
            node.nDescendants+=lastDesc
            if node.exploredChildren==nChildren:
                nextNode=node.up
                movingFrom=1
                lastL=node.minSoFar
                lastR=node.maxSoFar
                lastDesc=node.nDescendants
                if node!=t2:
                    if node.dist>minimumBLen:
                        if (lastR+1-lastL)==lastDesc:
                            if rooted or lastL>0:
                                if nodeTable[lastL][0]==lastL and nodeTable[lastL][1]==lastR: # ISCLUST(X,L,R --> is clust <L,R> or <R, L> part of the original tree?
                                    foundBranches+=1
                                    trueDist = branchLengthDict[(lastL, lastR)]
                                    KF += abs(trueDist - node.dist)
                                    RFL-= trueDist
                                elif nodeTable[lastR][0]==lastL and nodeTable[lastR][1]==lastR:
                                    foundBranches+=1
                                    branchLengthDict[(lastR, lastL)] = node.dist
                                    trueDist = branchLengthDict[(lastL, lastR)]
                                    KF += abs(trueDist - node.dist)
                                    RFL-= trueDist
                                else:
                                    missedBranches+=1
                                    RFL += node.dist
                            else: # re-root at leaf 0, so flip the values for the current branch if it contains leaf 0.
                                flippedL=lastR+1
                                flippedR=leafCount-1
                                if nodeTable[flippedL][0]==flippedL and nodeTable[flippedL][1]==flippedR:
                                    foundBranches+=1
                                    trueDist = branchLengthDict[(flippedL, flippedR)]
                                    KF += abs(trueDist - node.dist)
                                    RFL -= trueDist
                                elif nodeTable[flippedR][0]==flippedL and nodeTable[flippedR][1]==flippedR:
                                    foundBranches+=1
                                    trueDist = branchLengthDict[(flippedL, flippedR)]
                                    KF += abs(trueDist - node.dist)
                                    RFL -= trueDist
                                else:
                                    missedBranches+=1
                                    RFL += node.dist
                        else:
                            missedBranches+=1
                            RFL += node.dist
                elif rooted and addRootRFL: #if rooted take into account difference between roots.
                    if nodeTable[lastR][0] == lastL and nodeTable[lastR][1] == lastR:
                        trueDist = branchLengthDict[(lastL, lastR)]
                        KF += abs(trueDist - node.dist)
                        RFL -= trueDist
                    else: RFL += node.dist
            else:
                nextNode=node.children[node.exploredChildren]
                movingFrom=0
        node=nextNode
    if visitedLeaves<leafCount:
        print("There are leaves in the reference that have not been found in this new tree")
        return None, None, None, None, None, None
    #first value is number of differences, second value is max number of differences just in case one wants the normalized values;
    #the other values are there just in case on wants more detail.
    numDiffs=((numBranches-foundBranches)+missedBranches)
    RFL += KF
    return numDiffs, float(numDiffs)/(2*(leafCount-3)), leafCount, foundBranches, missedBranches, (numBranches-foundBranches), RFL

def traverseTopology(node, function, leafsOnly=False):
    """
    Traverse a topology and perform a function
    :param node: the root node.
    :param function: the function to be performed at every node. It should take only 'node' as argument
    :param leafsOnly: If True, the function will only be performed at the leaf node
    """
    nodesToVisit=[node]
    while nodesToVisit:
        newNode=nodesToVisit.pop()
        for c in newNode.children:
            nodesToVisit.append(c)
        if not newNode.children or not leafsOnly:
            function(newNode)

def getRoot(tree):
    root = tree
    while root.up != None:
        root = root.up
    return root

def counttotBLenAll(tree):
    """
    Traverse the tree to count the sum of all the branch lengths.
    """
    root = getRoot(tree)
    global totBLen
    totBLen = 0
    def counttotBLen(node):
        global totBLen
        totBLen += node.dist
    traverseTopology(root, counttotBLen)
    # print("Tot branch length: "+str(totBLen))
    return totBLen - root.dist



#run Robinson-Foulds distance calculations
if inputRFtrees!="":
    if not os.path.isfile(inputTree):
        print("Input tree in newick format "+inputTree+" not found, quitting MAPLE RF distance calculation. Use option --inputTree to specify a valid input newick tree file.")
        raise Exception("exit")
    if not os.path.isfile(inputRFtrees):
        print("Input trees in newick format "+inputRFtrees+" not found, quitting MAPLE RF distance calculation. Use option --inputRFtrees to specify valid file with input newick trees.")
        raise Exception("exit")
    lineSplit=outputFile.split("/")
    lineSplit[-1]=""
    outFolder="/".join(lineSplit)
    if not os.path.isdir(outFolder):
        print("Path to output file "+outFolder+" does not exist, quitting MAPLE RF calculation. Use option --output to specify a valid output file path and output file name.")
        raise Exception("exit")
    if os.path.isfile(outputFile+"_RFdistances.txt")  and (not overwrite):
        print("File "+outputFile+"_RFdistances.txt already exists, quitting MAPLE RF calculation. Use option --overwrite if you want to overwirte previous inference.")
        raise Exception("exit")
    if genomeLength: # provide args.genomeLength
        tree1 = readNewick(inputTree, divideBranchLengthsBy=genomeLength)[0]
    else:
        tree1=readNewick(inputTree)[0]
    print("Read input newick tree")
    leafNameDict, nodeTable, leafCount, numBranches = prepareTreeComparison(tree1,rooted=False)
    otherTrees=readNewick(inputRFtrees,multipleTrees=True)
    print("Read other input newick trees to be compared to the first one")
    file=open(outputFile+"_RFdistances.txt","w")
    file.write("RF\t"+"normalisedRF\t"+"leaves\t"+"foundBranches\t"+"missedBranches\t"+"notFoundBranches\n")
    for tree in otherTrees:
        numDiffs, normalisedRF, leafCount, foundBranches, missedBranches, notFoundBranches = RobinsonFouldsWithDay1985(tree,leafNameDict, nodeTable, leafCount, numBranches,rooted=False)
        file.write(str(numDiffs)+"\t"+str(normalisedRF)+"\t"+str(leafCount)+"\t"+str(foundBranches)+"\t"+str(missedBranches)+"\t"+str(notFoundBranches)+"\n")
    print("Comparison ended")
    raise Exception("exit")






minimumCarryOver=sys.float_info.min*(1e50)

if os.path.isfile(outputFile+"_tree.tree")  and (not overwrite):
    print("File "+outputFile+"_tree.tree already exists, quitting MAPLE tree inference. Use option --overwrite if you want to overwirte previous inference.")
    raise Exception("exit")
if not os.path.isfile(inputFile):
    print("Input file in Maple format "+inputFile+" not found, quitting MAPLE tree inference. Use option --input to specify a valid input file.")
    raise Exception("exit")
if refFile!="" and (not os.path.isfile(refFile)):
    print("Input reference fasta file "+refFile+" not found, quitting MAPLE tree inference. Use option --reference to specify a valid input reference file.")
    raise Exception("exit")
lineSplit=outputFile.split("/")
lineSplit[-1]=""
outFolder="/".join(lineSplit)
if not os.path.isdir(outFolder):
    print("Path to output file "+outFolder+" does not exist, quitting MAPLE tree inference. Use option --output to specify a valid output file path and output file name.")
    raise Exception("exit")
if inputTree!="":
    if not os.path.isfile(inputTree):
        print("Input tree in newick format "+inputTree+" not found, quitting MAPLE. Use option --inputTree to specify a valid input newick tree file.")
        raise Exception("exit")
    tree1=readNewick(inputTree,dirtiness=largeUpdate)[0]
    print("Read input newick tree")
    makeTreeBinary(tree1)


alleles={"A":0,"C":1,"G":2,"T":3}
allelesList=["A","C","G","T"]
allelesLow={"a":0,"c":1,"g":2,"t":3}
allelesUpOrLow={"a":0,"c":1,"g":2,"t":3,"A":0,"C":1,"G":2,"T":3}
allelesListLow=["a","c","g","t"]
ambiguities={"y":[0.0,0.5,0.0,0.5],"r":[0.5,0.0,0.5,0.0],"w":[0.5,0.0,0.0,0.5],"s":[0.0,0.5,0.5,0.0],"k":[0.0,0.0,0.5,0.5],"m":[0.5,0.5,0.0,0.0],"d":[1.0/3,0.0,1.0/3,1.0/3],"v":[1.0/3,1.0/3,1.0/3,0.0],"h":[1.0/3,1.0/3,0.0,1.0/3],"b":[0.0,1.0/3,1.0/3,1.0/3]}



#collect reference and calculate vector of cumulative numbers of nucleotides and calculate root frequencies.
def collectReference(fileName):
    file=open(fileName)
    line=file.readline()
    ref=""
    while line!="":
        line=file.readline()
        ref+=line.replace("\n","")
    ref=ref.lower()
    file.close()
    return ref



#Read input file
def readConciseAlignment(fileName,extractReference=True,ref="",extractNames=False):
    fileI=open(fileName)
    line=fileI.readline()
    if extractReference:
        line=fileI.readline()
        ref=""
        while line!="" and line[0]!=">":
            ref+=line.replace("\n","")
            line=fileI.readline()
        ref=ref.lower()
    nSeqs=0
    if extractNames:
        data={}
    else:
        data=[]
    while line!="" and line!="\n":
        seqList=[]
        if extractNames:
            name=line.replace(">","").replace("\n","")
        line=fileI.readline()
        pos=0
        while line!="" and line!="\n" and line[0]!=">":
            linelist=line.split()
            if len(linelist)>2:
                entry=(linelist[0].lower(),int(linelist[1]),int(linelist[2]))
            else:
                entry=(linelist[0].lower(),int(linelist[1]))
            if ref[entry[1]-1]==entry[0]:
                print("Mutation observed into reference nucleotide at position "+str(entry[1])+" , nucleotide "+entry[0]+". Wrong reference and/or diff file?")
                raise Exception("exit")
            if entry[1]<=pos:
                print("WARNING, at sample number "+str(nSeqs+1)+" found entry")
                print(line.replace("\n",""))
                print("which is inconsistent since the position is already represented by another entry:")
                print(seqList[-1])
                #print(" ignoring the entry in the file and proceeding with the analysis.")
                raise Exception("exit")
            else:
                seqList.append(entry)
                if len(entry)==2:
                    pos=entry[1]
                else:
                    pos=entry[1]+entry[2]-1
            line=fileI.readline()
        if extractNames:
            data[name]=seqList
        else:
            data.append(seqList)
        nSeqs+=1
    fileI.close()
    print(str(nSeqs)+" sequences in diff file.")
    if extractReference:
        return ref, data
    else:
        return data

runOnlyExample=False
if not runOnlyExample:
    #read sequence data from file
    if refFile=="":
        ref, data=readConciseAlignment(inputFile,extractNames=(inputTree!=""))
    else:
        ref=collectReference(refFile)
        data=readConciseAlignment(inputFile, extractReference=False, ref=ref,extractNames=(inputTree!=""))

lRef=len(ref)
print("Length of reference genome: "+str(lRef))
#vector to count how many bases of each type are cumulatively in the reference genome up to a certain position
cumulativeBases=[[0,0,0,0]]
for i in range(lRef):
    cumulativeBases.append(list(cumulativeBases[i]))
    cumulativeBases[i+1][allelesUpOrLow[ref[i]]]+=1
rootFreqs=[0.0,0.0,0.0,0.0]
rootFreqsLog=[0.0,0.0,0.0,0.0]
for i in range(4):
    rootFreqs[i]=cumulativeBases[-1][i]/float(lRef)
    rootFreqsLog[i]=log(rootFreqs[i])
refIndeces=[]
for i in range(lRef):
    refIndeces.append(allelesLow[ref[i]])
if model=="JC":
    rootFreqs=[0.25,0.25,0.25,0.25]
    rootFreqsLog=[log(0.25),log(0.25),log(0.25),log(0.25)]
oneMutBLen=1.0/lRef
minBLenSensitivity*=oneMutBLen
# minBLen=minBLen*oneMutBLen
# maxBLen=maxBLen*oneMutBLen
# minBLenForMidNode=minBLenForMidNode*oneMutBLen

range4=range(4)
#stricter thresholds than the input one for probabilities.
thresholdProb2=thresholdProb*thresholdProb
thresholdProb4=thresholdProb2*thresholdProb2

#IMPORTANT definitions of genome list entries structure.
#The first element of a genome list entry represnts its type: 0="A", 1="C", 2="G", 3="T", 4="R", 5="N", 6="O"
#the second element represents the last position of the stretch of genome that they represent;
# a value of 1 refers to the first position of the genome.
#For entries of type "O"/6 the last element is always a normalized vector of likelihoods (they always have to sum to 1.0).
#If present, another element after the position one represents the evolutionary distance since the considered type was observed (this is always missing with type "N"/5);
#If the element is not not present, this distance is assumed to be 0.0.
#If the distance (branch length) value is present, another distance value can also be present for entries of type <5 ;
#this is to account for the fact that the observation might have occurred on the other side of the phylogeny with respect to the root;
#when this additional branch length is also present, for example if the entry is (1, 234, 0.0001, 0.0002) then this means that a "C" was observed at genome position 234, and the observation
#is separated from the root by a distance of 0.0001, while a distance of 0.0002 separates the root from the current node (or position along a branch) considered.



#if probability mass is concentrated in one nucleotide, simplify the entry from "O" to another type.
def simplfy(vec,refA):
    maxP=0.0
    maxI=0
    numA=0
    for i in range4:
        if vec[i]>maxP:
            maxP=vec[i]
            maxI=i
        if vec[i]>thresholdProb:
            numA+=1
    if maxP<thresholdProb4:
        print("Inside simplify(), all values in vector are too small - something wrong?")
        print(vec)
        raise Exception("exit")
    if numA==1:
        if maxI==refA:
            return 4
        else:
            return maxI
    else:
        return 6



#Shorten genome list by merging together R entries that are mergeable
def shorten(vec):
    entryOld=vec[0]
    index=0
    while index<len(vec)-1:
        newVec=vec[index+1]
        if newVec[0]==4 and entryOld[0]==4 and len(newVec)==len(entryOld):
            if len(newVec)==2:
                vec.pop(index)
            elif abs(newVec[2]-entryOld[2])>thresholdProb:
                index+=1
                entryOld=vec[index]
            elif len(newVec)==3:
                vec.pop(index)
            elif abs(newVec[3]-entryOld[3])<thresholdProb:
                vec.pop(index)
            else:
                index+=1
                entryOld=vec[index]
        else:
            index+=1
            entryOld=vec[index]
    return


#Like shorten(), but specific to lower likelihoods that don't need to consider root frequencies.
#UNNECESSARY, one could always use shorten() instead of it.
def shortenLower(vec):
    entryOld=vec[0]
    index=0
    while index<len(vec)-1:
        newVec=vec[index+1]
        if newVec[0]==4 and entryOld[0]==4 and len(newVec)==len(entryOld):
            if len(newVec)==2:
                vec.pop(index)
            elif abs(newVec[2]-entryOld[2])>thresholdProb:
                index+=1
                entryOld=vec[index]
            else:
                vec.pop(index)
        else:
            index+=1
            entryOld=vec[index]
    return



#define partial likelihood vector for a sample given its data
def probVectTerminalNode(diffs):
    if diffs is None:
        probVect=[(5,lRef)]
        return probVect
    pos=1
    probVect=[]
    for m in diffs:
            currPos=m[1]
            if currPos>pos:
                #region where the node with branch length bLen is identical to the ref.
                probVect.append((4,currPos-1))
                pos=currPos
            if m[0]=="n" or m[0]=="-":
                if len(m)>2:
                    length=m[2]
                else:
                    length=1
                #region with no info, store last position and length.
                probVect.append((5,currPos+length-1))
                pos=currPos+length
            elif m[0] in allelesLow:
                #position at which node allele is sure but is different from the reference.
                probVect.append((allelesLow[m[0]],currPos))
                pos=currPos+1
            else:
                # non-"n" ambiguity character; for now interpret this as ambiguity instead of as a polymorphism.
                if onlyNambiguities:
                    # if user asks to, to make things easier, interpret any ambiguity as an "n".
                    probVect.append((5,currPos))
                else:
                    #otherwise, store as an "other" scenario, where each nucleotide has its own partial likelihood.
                    probVect.append((6,currPos,ambiguities[m[0]]))
                pos=currPos+1
    if pos<=lRef:
        probVect.append((4,lRef))
    return probVect


#Update and normalize the mutation rate matrix, given new mutation counts
def updateSubMatrix(pseudoMutCounts,model,oldMutMatrix): #,mutMatrixBLen,mutMatrixBLenAdjusted
    mutMatrix=[[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0]]
    if model=="UNREST":
        for i in range4:
            tot=0.0
            for j in range4:
                if j!=i:
                    mutMatrix[i][j]=pseudoMutCounts[i][j]/rootFreqs[i]
                    tot+=mutMatrix[i][j]
            mutMatrix[i][i]=-tot
    elif model=="GTR":
        for i in range4:
            tot=0.0
            for j in range4:
                if j!=i:
                    mutMatrix[i][j]=(pseudoMutCounts[i][j]+pseudoMutCounts[j][i])/rootFreqs[i]
                    tot+=mutMatrix[i][j]
            mutMatrix[i][i]=-tot
    else:
        print("Substitution model not recognised! Exiting")
        raise Exception("exit")
    totRate=-(rootFreqs[0]*mutMatrix[0][0]+ rootFreqs[1]*mutMatrix[1][1]+ rootFreqs[2]*mutMatrix[2][2]+ rootFreqs[3]*mutMatrix[3][3] )
    for i in range4:
        for j in range4:
            mutMatrix[i][j]=mutMatrix[i][j]/totRate
    matChange=0.0
    for i in range4:
        for j in range4:
            if j!=i:
                matChange+=abs(mutMatrix[i][j]-oldMutMatrix[i][j])
    if matChange > 0.001: #the matrix has changed significantly, update it
        for i in range4:
            for j in range4:
                oldMutMatrix[i][j]=mutMatrix[i][j]
        #the matrix has changed significantly, return True so that cumulative rates can be updated as well.
        return True
    else:
        return False







#merge two partial likelihood vectors, one from above, probVect1, and one from below, probVect2
#unlike appendProb(), this function is not used on a large part of the tree at each placement, but only in a small neighbourhood;
def mergeVectorsUpDown(probVect1,bLenUp,probVect2,bLenDown,mutMatrix,useRateVariation=False,mutMatrices=None, node1isleaf=False, node2isleaf=False):
    indexEntry1, indexEntry2, pos = 0, 0, 0
    probVect=[]
    entry1=probVect1[indexEntry1]
    entry2=probVect2[indexEntry2]
    while True:
        if entry1[0]==5:
            if entry2[0]==5:
                pos=min(entry1[1],entry2[1])
                probVect.append((5,pos))

            elif entry2[0]<5:
                pos=min(entry1[1],entry2[1])
                if len(entry2)==3:
                    if bLenDown:
                        probVect.append((entry2[0],pos,entry2[2]+bLenDown,0.0))
                    else:
                        probVect.append((entry2[0],pos,entry2[2],0.0))
                else:
                    if bLenDown:
                        probVect.append((entry2[0],pos,bLenDown,0.0))
                    else:
                        probVect.append((entry2[0],pos))
            else: # entry2 case "O", entry 1 is "N"
                if useRateVariation:
                    mutMatrix=mutMatrices[pos]
                pos+=1
                if len(entry2)==4:
                    totBLen=entry2[2]
                    if bLenDown:
                        totBLen+=bLenDown
                else:
                    totBLen=bLenDown
                newVect=[]
                if totBLen:
                    for i in range4:
                        tot=0.0
                        for j in range4:
                            tot+=mutMatrix[i][j]*entry2[-1][j]
                        tot*=totBLen
                        tot+=entry2[-1][i]
                        newVect.append(tot*rootFreqs[i])
                    totSum=sum(newVect)
                    for i in range4:
                        newVect[i]/=totSum
                    probVect.append((6,pos,newVect))
                else:
                    for i in range4:
                        newVect.append(entry2[-1][i]*rootFreqs[i])
                    totSum=sum(newVect)
                    for i in range4:
                        newVect[i]/=totSum
                    probVect.append((6,pos,newVect))
        elif entry2[0]==5: #entry2 is N
            if entry1[0]<5:
                pos=min(entry1[1],entry2[1])
                if len(entry1)==2:
                    if bLenUp:
                        probVect.append((entry1[0],pos,bLenUp))
                    else:
                        probVect.append((entry1[0],pos))
                elif len(entry1)==3:
                    if bLenUp:
                        probVect.append((entry1[0],pos,entry1[2]+bLenUp))
                    else:
                        probVect.append((entry1[0],pos,entry1[2]))
                else:
                    if bLenUp:
                        probVect.append((entry1[0],pos,entry1[2],entry1[3]+bLenUp))
                    else:
                        probVect.append((entry1[0],pos,entry1[2],entry1[3]))

            else: #entry1 is "O", entry 2 is "N"
                if useRateVariation:
                    mutMatrix=mutMatrices[pos]
                pos+=1
                if len(entry1)==4:
                    totBLen=entry1[2]
                    if bLenUp:
                        totBLen+=bLenUp
                elif bLenUp:
                    totBLen=bLenUp
                else:
                    totBLen=False
                newVect=[]
                if totBLen:
                    for i in range4:
                        tot=0.0
                        for j in range4:
                            tot+=entry1[-1][j]*mutMatrix[j][i]
                        tot*=totBLen
                        tot+=entry1[-1][i]
                        newVect.append(tot)
                    totSum=sum(newVect)
                    for i in range4:
                        newVect[i]/=totSum
                    probVect.append((6,pos,newVect))
                else:
                    probVect.append((6,pos,entry1[-1]))
        elif entry2[0]==entry1[0] and entry1[0]<5:
            pos=min(entry1[1],entry2[1])
            probVect.append((entry2[0],pos))
        else: #cases where the new genome list entry will likely be of type "O"
            if entry1[0]<5:
                if len(entry1)==2:
                    totLen1=bLenUp
                else:
                    totLen1=entry1[2]
                    if bLenUp:
                        totLen1+=bLenUp
                    if len(entry1)==4:
                        totLen1+=entry1[3]
            else:
                if len(entry1)==3:
                    totLen1=bLenUp
                else:
                    totLen1=entry1[2]
                    if bLenUp:
                        totLen1+=bLenUp

            if entry2[0]<5:
                if len(entry2)==2:
                    totLen2=bLenDown
                else:
                    totLen2=entry2[2]
                    if bLenDown:
                        totLen2+=bLenDown
            else:
                if len(entry2)==3:
                    totLen2=bLenDown
                else:
                    totLen2=entry2[2]
                    if bLenDown:
                        totLen2+=bLenDown
            if entry2[0]<5 and (not totLen2): #due to 0 distance, the entry will be of same type as entry2
                if (not totLen1) and entry1[0]<5:
                    #print("mergeVectorsUpDown() returning None 1")
                    #print(entry1)
                    #print(entry2)
                    #print(totLen1)
                    #print(totLen2)
                    return None
                pos=min(entry1[1],entry2[1])
                probVect.append((entry2[0],pos))
            elif entry1[0]<5 and (not totLen1): #due to 0 distance, the entry will be of same type as entry1
                pos=min(entry1[1],entry2[1])
                probVect.append((entry1[0],pos))
            elif entry1[0]<5:
                if useRateVariation:
                    mutMatrix=mutMatrices[pos]
                if entry1[0]==4:
                    i1=refIndeces[pos]
                else:
                    i1=entry1[0]
                newVec=[]
                if len(entry1)==4:
                    rootVec=list(rootFreqs)
                    for i in range4:
                        if i==i1:
                            rootVec[i]*=(1.0+mutMatrix[i1][i1]*(entry1[2]))
                        else:
                            rootVec[i]*=mutMatrix[i][i1]*(entry1[2])
                    if bLenUp:
                        lenToRoot=entry1[3]+bLenUp
                    else:
                        lenToRoot=entry1[3]
                    for j in range4:
                        tot=0.0
                        for i in range4:
                            tot+=mutMatrix[i][j]*rootVec[i]
                        tot*=lenToRoot
                        tot+=rootVec[j]
                        newVec.append(tot)
                else:
                    if totLen1:
                        for i in range4:
                            if i==i1:
                                newVec.append(1.0+mutMatrix[i][i]*totLen1)
                            else:
                                newVec.append(mutMatrix[i1][i]*totLen1)
                    else:
                        for i in range4:
                            if i==i1:
                                newVec.append(1.0)
                            else:
                                newVec.append(0.0)
                if entry2[0]==6: #entry 2 is "O" and entry1 is a nucleotide
                    for j in range4:
                        tot=0.0
                        if totLen2:
                            for i in range4:
                                tot+=mutMatrix[j][i]*entry2[-1][i]
                            tot*=totLen2
                        tot+=entry2[-1][j]
                        newVec[j]*=tot
                    sumV=sum(newVec)
                    for i in range4:
                        newVec[i]=newVec[i]/sumV
                    state =simplfy(newVec,refIndeces[pos])
                    pos+=1
                    if state==6:
                        probVect.append((6,pos,newVec))
                    else:
                        probVect.append((state,pos))
                else: #entry 1 and entry 2 are different nucleotides (possibly reference)
                    if entry2[0]==4:
                        i2=refIndeces[pos]
                    else:
                        i2=entry2[0]
                    if totLen2:
                        for i in range4:
                            if i==i2:
                                newVec[i]*=1.0+mutMatrix[i][i]*totLen2
                            else:
                                newVec[i]*=mutMatrix[i][i2]*totLen2
                    else:
                        for i in range4:
                            if i!=i2:
                                newVec[i]=0
                    sumV=sum(newVec)
                    if not sumV:
                        print("situation")
                        print(entry1)
                        print(entry2)
                        print(bLenUp)
                        print(bLenDown)
                        print(totLen1)
                        print(totLen2)
                        print(newVec)
                    for i in range4:
                        newVec[i]=newVec[i]/sumV
                    pos+=1
                    probVect.append((6,pos,newVec))
            else: #entry1[0]==6:
                if useRateVariation:
                    mutMatrix=mutMatrices[pos]
                if totLen1:
                    newVec=[]
                    for i in range4:
                        tot=0.0
                        for j in range4:
                            tot+=mutMatrix[j][i]*entry1[-1][j]
                        tot*=totLen1
                        tot+=entry1[-1][i]
                        newVec.append(tot)
                else:
                    newVec=list(entry1[-1])

                if entry2[0]==6:
                    if totLen2:
                        for i in range4:
                            tot=0.0
                            for j in range4:
                                tot+=mutMatrix[i][j]*entry2[-1][j]
                            tot*=totLen2
                            tot+=entry2[-1][i]
                            newVec[i]*=tot
                    else:
                        for i in range4:
                            newVec[i]*=entry2[-1][i]
                else: #entry1 is "O" and entry2 is a nucleotide
                    if entry2[0]==4:
                        i2=refIndeces[pos]
                    else:
                        i2=entry2[0]
                    if totLen2:
                        for i in range4:
                            if i==i2:
                                newVec[i]*=(1.0+mutMatrix[i][i]*totLen2)
                            else:
                                newVec[i]*=mutMatrix[i][i2]*totLen2
                    else:
                        for i in range4:
                            if i!=i2:
                                newVec[i]=0.0
                sumV=sum(newVec)
                if not sumV:
                    print("mergeVectorsUpDown() returning None 3")
                    print(entry1)
                    print(entry2)
                    print(sumV)
                    print(newVec)
                    return None
                for i in range4:
                    newVec[i]=newVec[i]/sumV
                state =simplfy(newVec,refIndeces[pos])
                pos+=1
                if state==6:
                    probVect.append((6,pos,newVec))
                else:
                    probVect.append((state,pos))

        if pos==lRef:
            break
        if pos==entry1[1]:
            indexEntry1+=1
            entry1=probVect1[indexEntry1]
        if pos==entry2[1]:
            indexEntry2+=1
            entry2=probVect2[indexEntry2]

    if verbose:
        print("Merged up-down ")
        print(probVect)
    #check if the final  probVect can be simplified by merging consecutive entries
    shorten(probVect)
    if verbose:
        print("Shortened up-down merging ")
        print(probVect)
    return probVect







#merge two lower child partial likelihood vectors to create a new one
#(and also calculate the logLk of the merging if necessary, which is currently only useful for the root but could also be useful to calculate the overall total likelihood of the tree).
def mergeVectors(probVect1,bLen1,probVect2,bLen2,mutMatrix,returnLK=False,useRateVariation=False,mutMatrices=None, node1isleaf=False, node2isleaf=False):
    indexEntry1, indexEntry2, pos = 0, 0, 0
    probVect=[]
    cumulPartLk=0.0
    entry1=probVect1[indexEntry1]
    entry2=probVect2[indexEntry2]
    end=min(entry1[1],entry2[1])
    while True:
        if entry1[0]==5:
            if entry2[0]==5:
                pos=min(entry1[1],entry2[1])
                probVect.append((5,pos))
            elif entry2[0]<5:
                pos=min(entry1[1],entry2[1])
                if len(entry2)==2:
                    if bLen2:
                        probVect.append((entry2[0],pos,bLen2))
                    else:
                        probVect.append((entry2[0],pos))
                else:
                    if bLen2:
                        probVect.append((entry2[0],pos,entry2[2]+bLen2))
                    else:
                        probVect.append((entry2[0],pos,entry2[2]))
            else: # case entry2 is "O" and entry1 is "N"
                pos+=1
                if len(entry2)==3:
                    if bLen2:
                        probVect.append((6,pos,bLen2,entry2[-1]))
                    else:
                        probVect.append((6,pos,entry2[-1]))
                else:
                    if bLen2:
                        probVect.append((6,pos,entry2[2]+bLen2,entry2[-1]))
                    else:
                        probVect.append((6,pos,entry2[2],entry2[-1]))
        elif entry2[0]==5: #entry2 is N
            if entry1[0]<5:
                pos=min(entry1[1],entry2[1])
                if len(entry1)==2:
                    if bLen1:
                        probVect.append((entry1[0],pos,bLen1))
                    else:
                        probVect.append((entry1[0],pos))
                else:
                    if bLen1:
                        probVect.append((entry1[0],pos,entry1[2]+bLen1))
                    else:
                        probVect.append((entry1[0],pos,entry1[2]))
            else: #entry1 is "O" and entry2 is "N"
                pos+=1
                if len(entry1)==3:
                    if bLen1:
                        probVect.append((6,pos,bLen1,entry1[-1]))
                    else:
                        probVect.append((6,pos,entry1[-1]))
                else:
                    if bLen1:
                        probVect.append((6,pos,entry1[2]+bLen1,entry1[-1]))
                    else:
                        probVect.append((6,pos,entry1[2],entry1[-1]))

        else: #entry1 and entry2 are not "N"
            if entry1[0]<5:
                if len(entry1)==2:
                    totLen1=bLen1
                else:
                    totLen1=entry1[2]
                    if bLen1:
                        totLen1+=bLen1
            else:
                if len(entry1)==3:
                    totLen1=bLen1
                else:
                    totLen1=entry1[2]
                    if bLen1:
                        totLen1+=bLen1

            if entry2[0]<5:
                if len(entry2)==2:
                    totLen2=bLen2
                else:
                    totLen2=entry2[2]
                    if bLen2:
                        totLen2+=bLen2
            else:
                if len(entry2)==3:
                    totLen2=bLen2
                else:
                    totLen2=entry2[2]
                    if bLen2:
                        totLen2+=bLen2

            if entry2[0]==entry1[0] and entry2[0]<5: #entry1 and entry2 are two identical nucleotides
                end=min(entry1[1],entry2[1])
                probVect.append((entry2[0],end))
                if returnLK:
                    if entry2[0]==4:
                        cumulPartLk+=(totLen1+totLen2)*(cumulativeRate[end]-cumulativeRate[pos])
                    else:
                        if useRateVariation:
                            cumulPartLk+=mutMatrices[pos][entry1[0]][entry1[0]]*(totLen1+totLen2)
                        else:
                            cumulPartLk+=nonMutRates[entry1[0]]*(totLen1+totLen2)
                pos=end
            elif (not totLen1) and (not totLen2) and entry1[0]<5 and entry2[0]<5: #0 distance between different nucleotides: merge is not possible
                if returnLK:
                    return None, float("-inf")
                else:
                    return None
            elif entry1[0]<5: #entry1 is a nucleotide
                if useRateVariation:
                    mutMatrix=mutMatrices[pos]
                if entry1[0]==4:
                    i1=refIndeces[pos]
                else:
                    i1=entry1[0]
                newVec=[]
                if totLen1:
                    for i in range4:
                        if i==i1:
                            newVec.append(1.0+mutMatrix[i][i]*totLen1)
                        else:
                            newVec.append(mutMatrix[i][i1]*totLen1)
                else:
                    newVec=[0.0,0.0,0.0,0.0]
                    newVec[i1]=1.0

                if entry2[0]==6: #entry1 is a nucleotide and entry2 is "O"
                    if totLen2:
                        for j in range4:
                            tot=0.0
                            for i in range4:
                                tot+=mutMatrix[j][i]*entry2[-1][i]
                            tot*=totLen2
                            tot+=entry2[-1][j]
                            newVec[j]*=tot
                    else:
                        for j in range4:
                            newVec[j]*=entry2[-1][j]
                    sumV=sum(newVec)
                    if not sumV:
                        if returnLK:
                            return None, float("-inf")
                        else:
                            return None
                    for i in range4:
                        newVec[i]=newVec[i]/sumV
                    state =simplfy(newVec,refIndeces[pos])
                    pos+=1
                    if state==6:
                        probVect.append((6,pos,newVec))
                    else:
                        probVect.append((state,pos))
                    if returnLK:
                        cumulPartLk+=log(sumV)
                else: #entry1 and entry2 are nucleotides
                    if entry2[0]==4:
                        i2=refIndeces[pos]
                    else:
                        i2=entry2[0]
                    if totLen2:
                        for i in range4:
                            if i==i2:
                                newVec[i]*=1.0+mutMatrix[i][i]*totLen2
                            else:
                                newVec[i]*=mutMatrix[i][i2]*totLen2
                        sumV=sum(newVec)
                        for i in range4:
                            newVec[i]=newVec[i]/sumV
                        state =simplfy(newVec,refIndeces[pos])
                        pos+=1
                        if state==6:
                            probVect.append((6,pos,newVec))
                        else:
                            probVect.append((state,pos))
                        #probVect.append((6,pos,newVec))
                        if returnLK:
                            cumulPartLk+=log(sumV)
                    else:
                        pos+=1
                        probVect.append((entry2[0],pos))
                        if returnLK:
                            cumulPartLk+=log(newVec[i2])

            else: #entry1[0]==6:
                if useRateVariation:
                    mutMatrix=mutMatrices[pos]
                if totLen1:
                    newVec=[]
                    for i in range4:
                        tot=0.0
                        for j in range4:
                            tot+=mutMatrix[i][j]*entry1[-1][j]
                        tot*=totLen1
                        tot+=entry1[-1][i]
                        newVec.append(tot)
                else:
                    newVec=list(entry1[-1])

                if entry2[0]==6:
                    if totLen2:
                        for i in range4:
                            tot=0.0
                            for j in range4:
                                tot+=mutMatrix[i][j]*entry2[-1][j]
                            tot*=totLen2
                            tot+=entry2[-1][i]
                            newVec[i]*=tot
                    else:
                        for i in range4:
                            newVec[i]*=entry2[-1][i]
                    sumV=sum(newVec)
                    if not sumV:
                        if returnLK:
                            return None, float("-inf")
                        else:
                            return None
                    for i in range4:
                        newVec[i]=newVec[i]/sumV
                    state =simplfy(newVec,refIndeces[pos])
                    pos+=1
                    if state==6:
                        probVect.append((6,pos,newVec))
                    else:
                        probVect.append((state,pos))
                    if returnLK:
                        cumulPartLk+=log(sumV)
                else: #entry2 is a nucleotide and entry1 is "O"
                    if entry2[0]==4:
                        i2=refIndeces[pos]
                    else:
                        i2=entry2[0]
                    if totLen2:
                        for i in range4:
                            if i==i2:
                                newVec[i]*=(1.0+mutMatrix[i][i]*totLen2)
                            else:
                                newVec[i]*=mutMatrix[i][i2]*totLen2
                        sumV=sum(newVec)
                        for i in range4:
                            newVec[i]=newVec[i]/sumV
                        state =simplfy(newVec,refIndeces[pos])
                        pos+=1
                        if state==6:
                            probVect.append((6,pos,newVec))
                        else:
                            probVect.append((state,pos))
                        if returnLK:
                            cumulPartLk+=log(sumV)
                    else:
                        if not newVec[i2]:
                            if returnLK:
                                return None, float("-inf")
                            else:
                                return None
                        pos+=1
                        probVect.append((entry2[0],pos))
                        if returnLK:
                            cumulPartLk+=log(newVec[i2])

        if pos==lRef:
            break
        if pos==entry1[1]:
            indexEntry1+=1
            entry1=probVect1[indexEntry1]
        if pos==entry2[1]:
            indexEntry2+=1
            entry2=probVect2[indexEntry2]

    if verbose:
        print("Merged vector at root ")
        print(probVect)
    #check if the final  probVect can be simplified by merging consecutive entries
    shorten(probVect)
    if verbose:
        print("Shortened root merging ")
        print(probVect)
    if returnLK:
        return probVect, cumulPartLk
    else:
        return probVect








#calculate the probability that results from combining a lower likelihood genome list of the root with root frequencies.
#Could maybe be combined with function rootVector() ?
def findProbRoot(probVect):
    logLK=0.0
    logFactor=1.0
    pos=0
    for entry in probVect:
        if entry[0]==4:
            for i in range4:
                logLK+=rootFreqsLog[i]*(cumulativeBases[entry[1]][i]-cumulativeBases[pos][i])
        elif entry[0]<4:
            logLK+=rootFreqsLog[entry[0]]
        elif entry[0]==6:
            tot=0.0
            for i in range4:
                tot+=rootFreqs[i]*entry[-1][i]
            logFactor*=tot
        pos=entry[1]
    logLK+=log(logFactor)
    return logLK











#for the root, take lower likelihood genome list probVect, and create an overall likelihood (or upper right or upper left) genome list by multiplying likelihoods by root frequencies.
def rootVector(probVect,bLen,mutMatrix,
               useRateVariation=False,mutMatrices=None,isLeaf=False):
    newProbVect=[]
    for entry in probVect:
        if entry[0]==5:
            newProbVect.append(entry)
        elif entry[0]==6:
            if len(entry)==4:
                totBLen=entry[2]
                if bLen:
                    totBLen+=bLen
            else:
                totBLen=bLen
            newVect=[]
            if totBLen:
                if useRateVariation:
                    mutMatrix=mutMatrices[entry[1]-1]
                for i in range4:
                    tot=0.0
                    for j in range4:
                        tot+=mutMatrix[i][j]*entry[-1][j]
                    tot*=totBLen
                    tot+=entry[-1][i]
                    newVect.append(tot*rootFreqs[i])
                totSum=sum(newVect)
                for i in range4:
                    newVect[i]/=totSum
                newProbVect.append((6,entry[1],newVect))
            else:
                for i in range4:
                    newVect.append(entry[-1][i]*rootFreqs[i])
                totSum=sum(newVect)
                for i in range4:
                    newVect[i]/=totSum
                newProbVect.append((6,entry[1],newVect))
        else:
            if len(entry)==3:
                if bLen:
                    newProbVect.append((entry[0],entry[1],entry[2]+bLen,0.0))
                else:
                    newProbVect.append((entry[0],entry[1],entry[2],0.0))
            else:
                if bLen:
                    newProbVect.append((entry[0],entry[1],bLen,0.0))
                else:
                    newProbVect.append((entry[0],entry[1]))
    return newProbVect







#function to add new mutation events in new sample to the pre-exisitng pseudocounts so to improve the estimate of the substitution rates.
# probVect1 is the genome list at the node where the appending happens; probVect2 is the genome list for the new sample.
def updatePesudoCounts(probVect1,probVect2,pseudoMutCounts):
    if model!="JC":
        indexEntry1, indexEntry2, pos = 0, 0, 0
        entry1=probVect1[indexEntry1]
        entry2=probVect2[indexEntry2]
        while True:
            if entry1[0]!=entry2[0] and entry1[0]<5 and entry2[0]<5:
                if entry1[0]==4:
                    pseudoMutCounts[refIndeces[pos]][entry2[0]]+=1
                elif entry2[0]==4:
                    pseudoMutCounts[entry1[0]][refIndeces[pos]]+=1
                else:
                    pseudoMutCounts[entry1[0]][entry2[0]]+=1
                pos+=1
            else:
                pos=min(entry1[1],entry2[1])

            if pos==lRef:
                break
            if pos==entry1[1]:
                indexEntry1+=1
                entry1=probVect1[indexEntry1]
            if pos==entry2[1]:
                indexEntry2+=1
                entry2=probVect2[indexEntry2]









numNodes=[0,0,0,0,0]

#Given a tree, and a final substitution rate matrix, re-calculate all genome lists within the tree according to this matrix.
# this is useful once the matrix estimation has finished, to make sure all genome lists replect this matrix.
def reCalculateAllGenomeLists(root,mutMatrix, checkExistingAreCorrect=False,countNodes=False,countPseudoCounts=False,pseudoMutCounts=None,data=None,useRateVariation=False,mutMatrices=None):
    #first pass to update all lower likelihoods.
    node=root
    #direction 0 means from parent, direction 1 means from a child
    lastNode=None
    direction=0
    while node!=None:
        if direction==0:
            if node.children:
                node=node.children[0]
            else:
                if data!=None:
                    if node.name in data:
                        node.probVect=probVectTerminalNode(data[node.name])
                        del data[node.name]
                    else:
                        print("Error: sample name "+node.name+" not found in the input sequence data - all samples in the input tree need to have a corresponding sequence entry.")
                        raise Exception("exit")
                if countNodes:
                    numNodes[0]+=1
                    for entry in node.probVect:
                        if entry[0]<4:
                            numNodes[1]+=1
                        elif entry[0]==4:
                            numNodes[2]+=1
                        elif entry[0]==5:
                            numNodes[3]+=1
                        else:
                            numNodes[4]+=1
                lastNode=node
                node=node.up
                direction=1
        else :
            if lastNode==node.children[0]:
                node=node.children[1]
                direction=0
            else:
                newLower=mergeVectors(node.children[0].probVect,node.children[0].dist,node.children[1].probVect,
                                      node.children[1].dist,mutMatrix,returnLK=False,useRateVariation=useRateVariation,
                                      mutMatrices=mutMatrices)#, node1isleaf= node1isleaf,node2isleaf=node2isleaf )
                if checkExistingAreCorrect:
                    if areVectorsDifferentDebugging(newLower,node.probVect):
                        print("Inside reCalculateAllGenomeLists(), new lower at node is different from the old one, and it shouldn't be.")
                        print(newLower)
                        print(node.probVect)
                        raise Exception("exit")
                if newLower==None:
                    if not node.children[0].dist:
                        nodeList=[]
                        updateBLen(nodeList,node,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                        updatePartials(nodeList,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                    elif not node.children[1].dist:
                        nodeList=[]
                        updateBLen(nodeList,node.children[1],mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                        updatePartials(nodeList,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                    else:
                        print("Strange, distances>0 but inconsistent lower genome list creation in reCalculateAllGenomeLists()")
                        raise Exception("exit")
                else:
                    node.probVect=newLower
                if countNodes:
                    numNodes[0]+=1
                    for entry in node.probVect:
                        if entry[0]<4:
                            numNodes[1]+=1
                        elif entry[0]==4:
                            numNodes[2]+=1
                        elif entry[0]==5:
                            numNodes[3]+=1
                        else:
                            numNodes[4]+=1

                lastNode=node
                node=node.up
                direction=1

    #now update the other genome lists for the root
    node=root
    # newVect=rootVector(node.probVect,False,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
    # if checkExistingAreCorrect:
    # 	if areVectorsDifferentDebugging(newVect,node.probVectTot):
    # 		print("new tot at root is different from the old one, and it shouldn't be.")
    # 		print(newVect)
    # 		print(node.probVectTot)
    # 		raise Exception("exit")
    # node.probVectTot=newVect
    if node.children:
        newVect=rootVector(node.children[1].probVect,node.children[1].dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
        if checkExistingAreCorrect:
            if areVectorsDifferentDebugging(newVect,node.probVectUpRight):
                print("new probVectUpRight at root is different from the old one, and it shouldn't be.")
                print(newVect)
                print(node.probVectUpRight)
                raise Exception("exit")
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
                raise Exception("exit")
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
                    #newVect=mergeVectorsUpDown(vectUp,node.dist,node.probVect,False,mutMatrix)
                    # if newVect==None:
                    # 	print("Strange, node.dist>0 but inconsistent total genome list creation in reCalculateAllGenomeLists()")
                    # 	raise Exception("exit")
                    # else:
                    # if checkExistingAreCorrect:
                    # 	if areVectorsDifferentDebugging(newVect,node.probVectTot):
                    # 		print("new probVectTot at node is different from the old one, and it shouldn't be.")
                    # 		print(newVect)
                    # 		print(node.probVectTot)
                    # 		print(node.children)
                    # 		print("VectUp and probVect:")
                    # 		print(vectUp)
                    # 		print(node.dist)
                    # 		print(node.probVect)
                    # 		print(node.dist*mutMatrix[1][3]*0.2466/0.7533)
                    # 		raise Exception("exit")
                    # node.probVectTot=newVect
                    newVect=mergeVectorsUpDown(vectUp,node.dist/2,node.probVect,node.dist/2,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                    if checkExistingAreCorrect:
                        if areVectorsDifferentDebugging(newVect,node.probVectTotUp):
                            print("new probVectTotUp at node is different from the old one, and it shouldn't be.")
                            print(newVect)
                            print(node.probVectTotUp)
                            raise Exception("exit")
                    node.probVectTotUp=newVect
                    # if node.dist>=2*minBLenForMidNode:
                    # 	createFurtherMidNodes(node,vectUp,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                else:
                    if hasattr(node, 'probVectTotUp'):# Changed this 20/9
                        if node.probVectTotUp: # debugging
                            print('node.probvecttotup should have been none')
                        node.probVectTotUp = None # Changed this 20/9
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
                            raise Exception("exit")
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
                                raise Exception("exit")
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
                            raise Exception("exit")
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
                                raise Exception("exit")
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






#Initialize the mutation rate matrix
#If a fixed rate matrix is needed for SARS-CoV-2, this is the nucleotide mutation rate matrix from De Maio et al 2021
#mutMatrix=[[0.0,0.039,0.310,0.123],[0.140,0.0,0.022,3.028],[0.747,0.113,0.0,2.953],[0.056,0.261,0.036,0.0]]
pseudoMutCounts=[[0.0,1.0,5.0,2.0],[2.0,0.0,1.0,40.0],[5.0,2.0,0.0,20.0],[2.0,3.0,1.0,0.0]]
if model=="JC":
    mutMatrix=[[-1.0,1.0/3,1.0/3,1.0/3],[1.0/3,-1.0,1.0/3,1.0/3],[1.0/3,1.0/3,-1.0,1.0/3],[1.0/3,1.0/3,1.0/3,-1.0]]
else:
    mutMatrix=[[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0]]
    updateSubMatrix(pseudoMutCounts,model,mutMatrix)

nonMutRates=[0,0,0,0]
for i in range4:
    nonMutRates[i]=mutMatrix[i][i]

cumulativeRate=[0.0]
for i in range(lRef):
    ind=refIndeces[i]
    cumulativeRate.append(cumulativeRate[-1]+nonMutRates[ind])



#In case an input tree is given, calculate all genome lists, pseudocounts and rates, and then recalculate genome lists again according to the new rates.
if inputTree!="":
    reCalculateAllGenomeLists(tree1,mutMatrix,countPseudoCounts=True,pseudoMutCounts=pseudoMutCounts,data=data)
    if updateSubMatrix(pseudoMutCounts,model,mutMatrix):
        for i in range(lRef):
            cumulativeRate[i+1]=cumulativeRate[i]+nonMutRates[refIndeces[i]]
        for i in range4:
            nonMutRates[i]=mutMatrix[i][i]
    reCalculateAllGenomeLists(tree1,mutMatrix)
    print("Genome list for initial tree and initial pseudocounts calculated.")



#Sort samples based on distance from reference, but punishing more isolated N's and ambiguity characters.
#more ambiguous sequences are placed last this way - this is useful since ambiguous sequences are harder to place and are more likely to be less informative (and so be removed from the analysis altogether)
#def distancesFromRefPunishNs(data,samples):
def distancesFromRefPunishNs(data,samples=None):
    sampleDistances=[]
    if samples==None:
        rangeInd=range(len(data))
    else:
        rangeInd=samples
    #for sample in samples:
    for diffIndex in rangeInd:
        diffs=data[diffIndex]
        pos=1
        comparisons=0
        diffNum=0
        for m in diffs:
            currPos=m[1]
            if currPos>pos:
                #region identical to the ref.
                comparisons+=currPos-pos
                pos=currPos
            if m[0]=="n" or m[0]=="-":
                if len(m)>2:
                    pos=currPos+m[2]
                else:
                    pos=currPos+1
                diffNum+=1
            elif m[0] in allelesLow:
                comparisons+=1
                diffNum+=1
                pos=currPos+1
            else:
                pos=currPos+1
                diffNum+=1
        if pos<=lRef:
            comparisons+=lRef+1-pos
        #sampleDistances.append((diffNum*1000+lRef-comparisons,sample))
        sampleDistances.append((diffNum*1000+lRef-comparisons,diffIndex))

    from operator import itemgetter
    print("Now doing sorting")
    sampleDistances.sort(reverse=True,key=itemgetter(0))
    return sampleDistances



#function to check is one sequence is less informative than another;
#returns 0 when the 2 sequences are not comparable, otherwise returns 1 if the first is more informative or if they are identical, and 2 otherwise.
def isMinorSequence(probVect1,probVect2):
    indexEntry1, indexEntry2, pos = 0, 0, 0
    entry1=probVect1[indexEntry1]
    entry2=probVect2[indexEntry2]
    found1bigger=False
    found2bigger=False
    while True:
        if entry1[0]!=entry2[0]:
            if entry1[0]==5:
                pos=min(entry1[1],entry2[1])
                found2bigger=True
            elif entry2[0]==5:
                pos=min(entry1[1],entry2[1])
                found1bigger=True
            elif entry1[0]==6:
                if entry2[0]==4:
                    i2=refIndeces[pos]
                else:
                    i2=entry2[0]
                if entry1[-1][i2]>0.1:
                    found2bigger=True
                else:
                    return 0
                pos+=1
            elif entry2[0]==6:
                if entry1[0]==4:
                    i1=refIndeces[pos]
                else:
                    i1=entry1[0]
                if entry2[-1][i1]>0.1:
                    found1bigger=True
                else:
                    return 0
                pos+=1
            else:
                return 0
        elif entry1[0]==6:
            for j in range4:
                if entry2[-1][j]>0.1 and entry1[-1][j]<0.1:
                    found1bigger=True
                elif entry1[-1][j]>0.1 and entry2[-1][j]<0.1:
                    found2bigger=True
            pos+=1
        else:
            pos=min(entry1[1],entry2[1])
        if found1bigger and found2bigger:
            return 0
        if pos==lRef:
            break
        if pos==entry1[1]:
            indexEntry1+=1
            entry1=probVect1[indexEntry1]
        if pos==entry2[1]:
            indexEntry2+=1
            entry2=probVect2[indexEntry2]

    if found1bigger:
        if found2bigger:
            return 0
        else:
            return 1
    else:
        if found2bigger:
            return 2
        else:
            return 1






#function to calculate likelihood cost of appending sample to parent node
#UNNECESSARY this can be considered an optimized version of function appendProbNode(), although appendProb() is likely marginally faster.
#as such, when rewriting the code, one may get rid of this function altogether and use just appendProbNode() instead to reduce amount of code required.
def appendProb(probVectP,probVectC,bLen,mutMatrix):
    if not bLen:
        bLen=0.0
    Lkcost, indexEntry1, indexEntry2, totalFactor, pos = 0.0, 0, 0, 1.0, 0
    entry1=probVectP[indexEntry1]
    entry2=probVectC[indexEntry2]
    end=min(entry1[1],entry2[1])
    contribLength=bLen
    while True:
        if entry2[0]==5: # case N
            pos=min(entry1[1],entry2[1])
        elif entry1[0]==5: # case N
            #if parent node is type "N", in theory we might have to calculate the contribution of root nucleotides;
            # however, if this node is "N" then every other node in the current tree is "N", so we can ignore this since this contribution cancels out in relative terms.
            pos=min(entry1[1],entry2[1])
        elif entry1[0]==4: # case entry1 is R
            if entry2[0]==4:
                end=min(entry1[1],entry2[1])
                if len(entry1)==2:
                    Lkcost+=bLen*(cumulativeRate[end]-cumulativeRate[pos])
                else:
                    contribLength=bLen+entry1[2]
                    if len(entry1)==3:
                        Lkcost+=contribLength*(cumulativeRate[end]-cumulativeRate[pos])
                    else:
                        #here contribution from root frequency gets added and subtracted so it's ignored
                        Lkcost+=(contribLength+entry1[3])*(cumulativeRate[end]-cumulativeRate[pos])
                pos=end
            elif entry2[0]==6:
                i1=refIndeces[pos]
                if len(entry1)==4:
                    contribLength=bLen+entry1[3]
                    if entry2[2][i1]>0.1:
                        contribLength+=entry1[2]
                        #here contribution from root frequency can also be also ignored
                        Lkcost+=nonMutRates[i1]*contribLength
                    else:
                        tot=0.0
                        for i in range4:
                            if i1==i:
                                tot2=rootFreqs[i]*(1.0+nonMutRates[i]*entry1[2])
                            else:
                                tot2=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]
                            tot3=0.0
                            for j in range4:
                                if entry2[2][j]>0.1:
                                    tot3+=mutMatrix[i][j]
                            tot3*=contribLength
                            if entry2[2][i]>0.1:
                                tot3+=1.0
                            tot+=tot2*tot3
                        totalFactor*=(tot/rootFreqs[i1])
                else:
                    if entry2[2][i1]>0.1:
                        if len(entry1)==3:
                            Lkcost+=nonMutRates[i1]*(bLen+entry1[2])
                        else:
                            Lkcost+=nonMutRates[i1]*bLen
                    else:
                        tot=0.0
                        for j in range4:
                            if entry2[2][j]>0.1:
                                tot+=mutMatrix[i1][j]
                        if len(entry1)==3:
                            totalFactor*=tot*(bLen+entry1[2])
                        else:
                            totalFactor*=tot*bLen
                pos+=1

            else: #entry1 is reference and entry2 is a different but single nucleotide
                if len(entry1)==2:
                    totalFactor*=mutMatrix[refIndeces[pos]][entry2[0]]*bLen
                elif len(entry1)==3:
                    totalFactor*=mutMatrix[refIndeces[pos]][entry2[0]]*(bLen+entry1[2])
                else:
                    i1=refIndeces[pos]
                    i2=entry2[0]
                    totalFactor*=((rootFreqs[i1]*mutMatrix[i1][i2]*(bLen+entry1[3])*(1.0+nonMutRates[i1]*entry1[2])+rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]*(1.0+nonMutRates[i2]*(bLen+entry1[3])))/rootFreqs[i1])
                pos+=1

        # entry1 is of type "O"
        elif entry1[0]==6:
            if len(entry1)==3:
                bLen13=bLen
            else:
                bLen13=bLen+entry1[2]
            if entry2[0]==6:
                tot=0.0
                for j in range4:
                    tot2=0.0
                    for j2 in range4:
                        if entry2[2][j2]>0.1:
                            tot2+=mutMatrix[j][j2]
                    tot2*=bLen13
                    if entry2[2][j]>0.1:
                        tot2+=1.0
                    tot+=tot2*entry1[-1][j]
                totalFactor*=tot
            else:
                if entry2[0]==4:
                    i2=refIndeces[pos]
                else:
                    i2=entry2[0]
                totalFactor*=(entry1[-1][i2]+bLen13*(entry1[-1][0]*mutMatrix[0][i2]+entry1[-1][1]*mutMatrix[1][i2]+entry1[-1][2]*mutMatrix[2][i2]+entry1[-1][3]*mutMatrix[3][i2]))
            pos+=1

        else: #entry1 is a non-ref nuc
            i1=entry1[0]
            if entry2[0]==i1:
                if len(entry1)==2:
                    Lkcost+=nonMutRates[i1]*bLen
                elif len(entry1)==3:
                    Lkcost+=nonMutRates[i1]*(bLen+entry1[2])
                else:
                    Lkcost+=nonMutRates[i1]*(bLen+entry1[2]+entry1[3])
            else: #entry1 and entry2 are of different types
                if entry2[0]==6:
                    if len(entry1)==4:
                        bLen15=bLen+entry1[3]
                        if entry2[2][i1]>0.1:
                            Lkcost+=nonMutRates[i1]*(bLen15+entry1[2])
                        else:
                            tot=0.0
                            for i in range4:
                                if i1==i:
                                    tot2=rootFreqs[i]*(1.0+nonMutRates[i1]*entry1[2])
                                else:
                                    tot2=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]
                                tot3=0.0
                                for j in range4:
                                    if entry2[2][j]>0.1:
                                        tot3+=mutMatrix[i][j]
                                if entry2[2][i]>0.1:
                                    tot+=tot2*(1.0+bLen15*tot3)
                                else:
                                    tot+=tot2*bLen15*tot3
                            totalFactor*=(tot/rootFreqs[i1])
                    else:
                        if entry2[2][i1]>0.1:
                            if len(entry1)==2:
                                Lkcost+=nonMutRates[i1]*bLen
                            else:
                                Lkcost+=nonMutRates[i1]*(bLen+entry1[2])
                        else:
                            tot=0.0
                            for j in range4:
                                if entry2[2][j]>0.1:
                                    tot+=mutMatrix[i1][j]
                            if len(entry1)==2:
                                totalFactor*=tot*bLen
                            else:
                                totalFactor*=tot*(bLen+entry1[2])
                #entry2 is a nucleotide type (like entry1)
                else:
                    if entry2[0]==4:
                        i2=refIndeces[pos]
                    else:
                        i2=entry2[0]
                    if len(entry1)==2:
                        totalFactor*=mutMatrix[i1][i2]*bLen
                    elif len(entry1)==3:
                        totalFactor*=mutMatrix[i1][i2]*(bLen+entry1[2])
                    else:
                        #here we ignore contribution of non-parsimonious mutational histories
                        totalFactor*=((rootFreqs[i1]*mutMatrix[i1][i2]*(bLen+entry1[3])*(1.0+nonMutRates[i1]*entry1[2]) + rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]*(1.0+nonMutRates[i2]*(bLen+entry1[3])) )/rootFreqs[i1])
            pos+=1

        if totalFactor<=minimumCarryOver:
            if totalFactor<sys.float_info.min:
                return float("-inf")
            Lkcost+=log(totalFactor)
            totalFactor=1.0

        if pos==lRef:
            break
        if pos==entry1[1]:
            indexEntry1+=1
            entry1=probVectP[indexEntry1]
        if pos==entry2[1]:
            indexEntry2+=1
            entry2=probVectC[indexEntry2]
    return Lkcost+log(totalFactor)





#number of samples that could have been placed as major of another sample but weren't due to sample placement order
totalMissedMinors=[0]








#function traversing the tree to find the best node in the tree where to re-append the given subtree (rooted at node.children[child]) to improve the topology of the current tree.
# bestLKdiff is the best likelihood cost found for the current placement (after optimizing the branch length).
# removedBLen is such branch length that optimizes the current placement - it will be used to place the subtree attached at other nodes of the tree.

def findBestParentTopology(node,child,bestLKdiff,removedBLen,mutMatrix,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,useRateVariation=False,mutMatrices=None):
    #currentLK = bestLKdiff
    bestNode=node
    bestNodes=[]
    nodesToVisit=[]
    removedPartials=node.children[child].probVect
    removedPartialsIsLeaf = (node.children[child].children==[])
    originalLK = bestLKdiff
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
        originalBLens = (node.dist, node.children[1 - child].dist, removedBLen) #v0.1.9
        originalPlacement = node.children[1 - child]
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
            originalPlacement = node.children[1 - child].children[0] #v0.1.9
            originalBLens = (0.0, node.children[1 - child].children[0].dist, removedBLen) #v0.1.9
        else:
            originalPlacement = node.children[1 - child] #v0.1.9
            originalBLens = (0.0, node.children[1 - child].dist, removedBLen) #v0.1.9

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
                midProb=appendProbNode(midTot,removedPartials,removedBLen,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,
                                       node2isleaf=removedPartialsIsLeaf) # because removedPartials=node.children[child].probVect
                #notably, this has an improved value for the first loop over NodesToVisit, using node.children[1-child].probVect
                if midProb>bestLKdiff:
                    bestLKdiff=midProb
                    bestNode=t1
                    failedPasses=0
                if midProb>bestLKdiff-thresholdLogLKtopology/factorOptimizePlacementLKvsSearchLK:
                    #if needsUpdating, then add to the tuple also the information on the up and down genome lists to use to recalculate intermediate genome lists at varying branch lengths
                    if needsUpdating:
                        bestNodes.append((t1,midProb,passedPartials,t1.probVect,distance,midTot))
                    else:
                        bestNodes.append((t1,midProb))
                if midProb<(lastLK-thresholdLogLKconsecutivePlacement): #placement at current node is considered failed if placement likelihood is not improved by a certain margin compared to best placement so far for the nodes above it.
                    failedPasses+=1
            else:
                if hasattr(t1, 'probVectTotUp'): # Changed this 20/9
                    t1.probVectTotUp = None # Changed this 20/9
                midProb=lastLK

            #keep crawling down into children nodes unless the stop criteria for the traversal are satisfied.
            traverseChildren=False
            if strictTopologyStopRules:
                if failedPasses<=allowedFailsTopology and midProb>(bestLKdiff-thresholdLogLKtopology) and t1.children:
                    traverseChildren=True
            else:
                if failedPasses<=allowedFailsTopology or midProb>(bestLKdiff-thresholdLogLKtopology):
                    if t1.children:
                        traverseChildren=True
            if traverseChildren:
                #if (failedPasses<=allowedFailsTopology or nodeProb>(bestLKdiff-thresholdLogLKtopology) ) and len(t1.children)==2:
                child=t1.children[0]
                otherChild=t1.children[1]
                if needsUpdating:
                    vectUpRight=mergeVectorsUpDown(passedPartials,distance,otherChild.probVect,otherChild.dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,
                                                   node2isleaf=otherChild.children==[] )
                else:
                    vectUpRight=t1.probVectUpRight
                if vectUpRight!=None:
                    nodesToVisit.append((child,0,vectUpRight,child.dist,needsUpdating,midProb,failedPasses, False))
                child=t1.children[1]
                otherChild=t1.children[0]
                if needsUpdating:
                    vectUpLeft=mergeVectorsUpDown(passedPartials,distance,otherChild.probVect,otherChild.dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,
                                                  node2isleaf=otherChild.children == [])
                else:
                    vectUpLeft=t1.probVectUpLeft
                if vectUpLeft!=None:
                    nodesToVisit.append((child,0,vectUpLeft,child.dist,needsUpdating,midProb,failedPasses, False)) #no leaf

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
                midProb=appendProbNode(midTot,removedPartials,removedBLen,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices, node2isleaf=removedPartialsIsLeaf)
                if midProb>bestLKdiff:
                    bestLKdiff=midProb
                    bestNode=t1
                    failedPasses=0
                if midProb>=(bestLKdiff-thresholdLogLKtopology/factorOptimizePlacementLKvsSearchLK):
                    #if needsUpdating, then add to the tuple also the information on the up and down genome lists to use to recalculate intermediate genome lists at varying branch lengths
                    if needsUpdating:
                        bestNodes.append((t1,midProb,vectUp,midBottom,t1.dist,midTot)) #VectUp can never be a leaf
                    else:
                        bestNodes.append((t1,midProb))
                if midProb<(lastLK-thresholdLogLKconsecutivePlacement): #placement at current node is considered failed if placement likelihood is not improved by a certain margin compared to best placement so far for the nodes above it.
                    failedPasses+=1
            else:
                midProb=lastLK
                if hasattr(t1, 'probVectTotUp'): # Changed this 20/9
                    t1.probVectTotUp = None # Changed this 20/9

            #testing for the stop rule of the traversal process
            keepTraversing=False
            if strictTopologyStopRules:
                if failedPasses<=allowedFailsTopology and midProb>(bestLKdiff-thresholdLogLKtopology):
                    keepTraversing=True
            else:
                if failedPasses<=allowedFailsTopology or midProb>(bestLKdiff-thresholdLogLKtopology):
                    keepTraversing=True
            if keepTraversing:
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
                        nodesToVisit.append((otherChild,0,vectUp,otherChild.dist,needsUpdating,midProb,failedPasses, False))
                    #now pass the crawling up to the parent node
                    if needsUpdating:
                        if midBottom==None:
                            midBottom=mergeVectors(otherChild.probVect,otherChild.dist,passedPartials,distance,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,
                                                   node1isleaf=otherChild.children==[], node2isleaf=passedPartialsIsLeaf )#, node1isleaf= node1isleaf,node2isleaf=node2isleaf )
                            if midBottom==None:
                                continue
                    else:
                        midBottom=t1.probVect
                    nodesToVisit.append((t1.up,upChild+1,midBottom,t1.dist,needsUpdating,midProb,failedPasses, False))
                #now consider case of root node
                else:
                    if needsUpdating:
                        vectUp=rootVector(passedPartials,distance,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices, isLeaf=passedPartialsIsLeaf)
                    else:
                        if direction==1:
                            vectUp=t1.probVectUpLeft
                        else:
                            vectUp=t1.probVectUpRight
                    nodesToVisit.append((otherChild,0,vectUp,otherChild.dist,needsUpdating,midProb,failedPasses, False))

    #Initial exploration is finished.
    #Now, for each branch within threshold likelihood distance from the best found, optimize branch lengths.
    #Use optimized scores to select final best branch
    bestBranchLengths = originalBLens #version 0.1.9

    #bestBranchLengths=(None,None,None)
    bestScore= bestLKdiff #currentLK -1 #float('inf') #
    compensanteForBranchLengthChange=True #true by default, also add the LK cost for changing the Blen in the initial tree.
    secondBranchLengthOptimizationRound=False #false by default
    if not bestNodes: # v0.1.9
        return originalPlacement, originalLK, originalBLens
    BLenHaveBeenOptimized = False # v0.1.9
    for nodePair in bestNodes: # it is strange that bestScore is not replaced in the following part?
        score=nodePair[1]
        if score>=bestLKdiff-thresholdLogLKtopology/factorOptimizePlacementLKvsSearchLK:
            t1=nodePair[0]
            if len(nodePair)==2:
                #optimize branch lengths of appendage
                if t1==t1.up.children[0]:
                    upVect=t1.up.probVectUpRight
                else:
                    upVect=t1.up.probVectUpLeft
                downVect=t1.probVect
                distance=t1.dist
                midTot=t1.probVectTotUp
                downVectIsLeaf = (t1.children==[])
            else:
                upVect=nodePair[2]
                downVect=nodePair[3]
                distance=nodePair[4]
                midTot=nodePair[5]
                downVectIsLeaf = (nodePair[0].children==[])

            bestAppendingLength=estimateBranchLengthWithDerivative(midTot,removedPartials,mutMatrix, node2isleaf=removedPartialsIsLeaf) #midTot is no leaf


            #now optimize appending location
            midLowerVector=mergeVectors(downVect,distance/2,removedPartials,bestAppendingLength,mutMatrix,
                                        node1isleaf=downVectIsLeaf,node2isleaf=removedPartialsIsLeaf) # When changing node1isleaf or node2isleaf, it remains 0.
            bestTopLength=estimateBranchLengthWithDerivative(upVect,midLowerVector,mutMatrix, node2isleaf=False) #midLowerVector is no leaf

            midTopVector=mergeVectorsUpDown(upVect,bestTopLength,removedPartials,bestAppendingLength,mutMatrix,
                                            node2isleaf=removedPartialsIsLeaf)
            bestBottomLength=estimateBranchLengthWithDerivative(midTopVector,downVect,mutMatrix,
                                                                node2isleaf=downVectIsLeaf) #with the original midTopVect it estimates 0.000709709809765853
            newMidVector=mergeVectorsUpDown(upVect,bestTopLength,downVect,bestBottomLength,mutMatrix,
                                            node2isleaf=downVectIsLeaf)
            appendingCost=appendProbNode(newMidVector,removedPartials,bestAppendingLength,mutMatrix, node2isleaf=removedPartialsIsLeaf) #shouldn't rate variation be used here and the original one?

            if compensanteForBranchLengthChange: #if wanted, do a more thorough examination of the appending cost, taking into account a possible change in branch length of the branch on which to be appended.
                initialCost=appendProbNode(upVect,downVect,distance,mutMatrix, node2isleaf=downVectIsLeaf)
                newPartialCost=appendProbNode(upVect,downVect,bestBottomLength+bestTopLength,mutMatrix, node2isleaf=downVectIsLeaf)
                optimizedScore=appendingCost+newPartialCost-initialCost
            else:
                optimizedScore=appendingCost
            if optimizedScore>=bestScore: # if this statement is never true, bestBranchLengths remains (None,none, none)
                BLenHaveBeenOptimized = True #v0.1.9
                bestNode=t1
                bestScore=optimizedScore
                bestBranchLengths=(bestTopLength,bestBottomLength,bestAppendingLength)
    if not BLenHaveBeenOptimized: #v0.1.9
        bestBranchLengths = (bestNode.dist / 2, bestNode.dist / 2, removedBLen)
    return bestNode, bestScore, bestBranchLengths



#function to find the best node in the tree where to append the new sample; traverses the tree and tries to append the sample at each node and mid-branch nodes,
# but stops traversing when certain criteria are met.
def findBestParentForNewSample(root,diffs,sample,mutMatrix):
    bestNodes=[]
    bestNode=root
    bestBranchLengths=(False,False,oneMutBLen)
    if not root.children: #check if the new leaf is strictly less informative than already placed leaf
        comparison=isMinorSequence(root.probVect,diffs)
        if comparison==1:
            root.minorSequences.append(sample)
            return root, 1.0, None
        elif comparison==2:
            totalMissedMinors[0]+=1
    rootVect=rootVector(root.probVect,False,mutMatrix)
    bestLKdiff=appendProb(rootVect,diffs,oneMutBLen,mutMatrix)
    nodesToVisit=[]
    for child in root.children:
        nodesToVisit.append((child,bestLKdiff,0))
    while nodesToVisit:
        t1,parentLK,failedPasses=nodesToVisit.pop()
        if not t1.children: #check if the new leaf is strictly less informative than already placed leaf
            comparison=isMinorSequence(t1.probVect,diffs)
            if comparison==1:
                t1.minorSequences.append(sample)
                return t1, 1.0, None
            elif comparison==2:
                totalMissedMinors[0]+=1

        if t1.dist and t1.up!=None: # try first placing as a descendant of the mid-branch point of the branch above the current node.
            LKdiff=appendProb(t1.probVectTotUp,diffs,oneMutBLen,mutMatrix)
            if LKdiff>=bestLKdiff:
                bestLKdiff=LKdiff
                bestNode=t1
                failedPasses=0
                bestNodes.append((t1,LKdiff))
            elif LKdiff>bestLKdiff-thresholdLogLK/factorOptimizePlacementLKvsSearchLK:
                bestNodes.append((t1,LKdiff))
            if LKdiff<(parentLK-thresholdLogLKconsecutivePlacement): #placement at current node is considered failed if placement likelihood is not improved by a certain margin compared to best placement so far for the nodes above it.
                failedPasses+=1
        else:
            LKdiff=parentLK
        #keep trying to place at children nodes, unless placement has failed too many times already or the likelihood cost suboptimality is above a certain threshold
        if strictInitialStopRules:
            if failedPasses<=allowedFails and LKdiff>(bestLKdiff-thresholdLogLK):
                for c in t1.children:
                    nodesToVisit.append((c,LKdiff,failedPasses))
        else:
            if failedPasses<=allowedFails or LKdiff>(bestLKdiff-thresholdLogLK):
                for c in t1.children:
                    nodesToVisit.append((c,LKdiff,failedPasses))
    #Initial exploration is finished.
    #Now, for each branch within threshold likelihood distance from the best found, optimize branch lengths.
    #Use optimized scores to select final best branch
    if bestNode!=root:
        bestBranchLengths=(bestNode.dist/2,bestNode.dist/2,oneMutBLen)
    bestScore=bestLKdiff
    compensanteForBranchLengthChange=True
    secondBranchLengthOptimizationRound=False
    for nodePair in bestNodes:
        score=nodePair[1]
        if score>=bestLKdiff-thresholdLogLK/factorOptimizePlacementLKvsSearchLK:
            node=nodePair[0]
            #optimize branch lengths of appendage
            if node==node.up.children[0]:
                upVect=node.up.probVectUpRight
            else:
                upVect=node.up.probVectUpLeft
            bestAppendingLength=estimateBranchLengthWithDerivative(node.probVectTotUp,diffs,mutMatrix)
            #now optimize appending location
            midLowerVector=mergeVectors(node.probVect,node.dist/2,diffs,bestAppendingLength,mutMatrix)
            bestTopLength=estimateBranchLengthWithDerivative(upVect,midLowerVector,mutMatrix)
            midTopVector=mergeVectorsUpDown(upVect,bestTopLength,diffs,bestAppendingLength,mutMatrix)
            bestBottomLength=estimateBranchLengthWithDerivative(midTopVector,node.probVect,mutMatrix)
            newMidVector=mergeVectorsUpDown(upVect,bestTopLength,node.probVect,bestBottomLength,mutMatrix)
            if secondBranchLengthOptimizationRound: #if wanted, do a second round of branch length optimization
                bestAppendingLength=estimateBranchLengthWithDerivative(newMidVector,diffs,mutMatrix)
                midLowerVector=mergeVectors(node.probVect,bestBottomLength,diffs,bestAppendingLength,mutMatrix)
                bestTopLength=estimateBranchLengthWithDerivative(upVect,midLowerVector,mutMatrix)
                midTopVector=mergeVectorsUpDown(upVect,bestTopLength,diffs,bestAppendingLength,mutMatrix)
                bestBottomLength=estimateBranchLengthWithDerivative(midTopVector,node.probVect,mutMatrix)
                newMidVector=mergeVectorsUpDown(upVect,bestTopLength,node.probVect,bestBottomLength,mutMatrix)
            appendingCost=appendProb(newMidVector,diffs,bestAppendingLength,mutMatrix)
            if compensanteForBranchLengthChange: #if wante, do a more thorough examination of the appending cost, taking into account a possible change in branch length of the branch on which to be appended.
                initialCost=appendProbNode(upVect,node.probVect,node.dist,mutMatrix)
                newPartialCost=appendProbNode(upVect,node.probVect,bestBottomLength+bestTopLength,mutMatrix)
                optimizedScore=appendingCost+newPartialCost-initialCost
            else:
                optimizedScore=appendingCost
            if optimizedScore>=bestScore:
                bestNode=node
                bestScore=optimizedScore
                bestBranchLengths=(bestTopLength,bestBottomLength,bestAppendingLength)

    return bestNode, bestScore, bestBranchLengths







#Check if two genome lists represent the same partial likelihoods or not.
#This is used to traverse the tree and only update genome lists if they changed after a change in the tree.
def areVectorsDifferent(probVect1,probVect2):
    if probVect2==None:
        return True
    indexEntry1, indexEntry2, pos = 0, 0, 0
    entry1=probVect1[indexEntry1]
    entry2=probVect2[indexEntry2]
    while True:
        if entry1[0]!=entry2[0]:
            return True
        if len(entry1)!=len(entry2):
            return True
        if entry1[0]<5:
            if len(entry1)>2:
                if abs(entry1[2] - entry2[2])>thresholdProb:
                    return True
                if len(entry1)==4:
                    if abs(entry1[3] - entry2[3])>thresholdProb:
                        return True
        if entry1[0]==6:
            if len(entry1)==4:
                if abs(entry1[2] - entry2[2])>thresholdProb:
                    return True
            for i in range4:
                diffVal=abs(entry1[-1][i] - entry2[-1][i])
                #example thresholdDiffForUpdate=0.01
                #example thresholdFoldChangeUpdate=1.5
                if diffVal:
                    if (not entry1[-1][i]) or (not entry2[-1][i]):
                        return True
                    if diffVal>thresholdDiffForUpdate or (diffVal>thresholdProb and ( (diffVal/entry1[-1][i]>thresholdFoldChangeUpdate)  or  (diffVal/entry2[-1][i]>thresholdFoldChangeUpdate))):
                        return True
        #update pos, end, etc
        pos=min(entry1[1],entry2[1])
        if pos==lRef:
            break
        if pos==entry1[1]:
            indexEntry1+=1
            entry1=probVect1[indexEntry1]
        if pos==entry2[1]:
            indexEntry2+=1
            entry2=probVect2[indexEntry2]
    return False



#Check if two genome lists represent the same partial likelihoods or not.
#This is a less strict version used for debugging purposes.
def areVectorsDifferentDebugging(probVect1,probVect2,threshold=0.00001):
    if probVect2==None:
        print('none case')
        return True
    indexEntry1, indexEntry2, pos = 0, 0, 0
    entry1=probVect1[indexEntry1]
    entry2=probVect2[indexEntry2]
    while True:
        if entry1[0]!=6 and entry2[0]!=6:
            if entry1[0]!=entry2[0]:
                print('diff type case')
                return True
            if len(entry1)!=len(entry2):
                if len(entry2) > len(entry1) and type(entry2[-1])==bool: #1 of them is with error rates
                    entry2 = entry2[:(len(entry1)-len(entry2))]
                elif len(entry2) < len(entry1) and type(entry1[-1]) == bool:
                    entry1 = entry1[:(len(entry2)-len(entry1))]
                else:
                    print('diff len case')
                    return True
            if entry1[0]<5:
                if len(entry1)>2:
                    if abs(entry1[2] - entry2[2])>thresholdProb:
                        print('diff bLen case')
                        return True
                    if len(entry1)==4:
                        if abs(entry1[3] - entry2[3])>thresholdProb:
                            print('diff bLen2 case')
                            return True
        elif entry1[0]==6 and entry2[0]==6 :
            if len(entry1)==4 and len(entry2)==4 :
                if abs(entry1[2] - entry2[2])>thresholdProb:
                    print("66 bLen case")
                    return True
            elif len(entry1)!=len(entry2):
                print("66 entry len")
                return True
            for i in range4:
                diffVal=abs(entry1[-1][i] - entry2[-1][i])
                if diffVal:
                    if (not entry1[-1][i]) or (not entry2[-1][i]):
                        print("66 0 like case")
                        return True
                    if diffVal>0.01 or (diffVal>threshold and ( (diffVal/entry1[-1][i]>thresholdFoldChangeUpdate)  or  (diffVal/entry2[-1][i]>thresholdFoldChangeUpdate))):
                        print("66 like diff case")
                        print(entry1)
                        print(entry2)
                        return True
        else:
            if not (entry1[0]==5 and entry2[0]==5):
                if (entry1[0]==5 or entry2[0]==5):
                    print("N case")
                    return True
                if entry1[0]<5:
                    if entry1[0]==4:
                        i1=refIndeces[pos]
                    else:
                        i1=entry1[0]
                    if entry2[-1][i1]+threshold<1.0:
                        print("i1 case")
                        print(entry1)
                        print(entry2)
                        return True
                elif entry2[0]<5:
                    if entry2[0]==4:
                        i2=refIndeces[pos]
                    else:
                        i2=entry2[0]
                    if entry1[-1][i2]+threshold<1.0:
                        print("i2 case")
                        print(entry1)
                        print(entry2)
                        return True
        #update pos, end, etc
        pos=min(entry1[1],entry2[1])
        if pos==lRef:
            break
        if pos==entry1[1]:
            indexEntry1+=1
            entry1=probVect1[indexEntry1]
        if pos==entry2[1]:
            indexEntry2+=1
            entry2=probVect2[indexEntry2]
    return False








#if updating genome lists in updatePartials() creates an inconsistency, this function can increase the length of a 0-length branch to resolve the inconsistency.
#In doing so, it updates, the input list of nodes to visit and update.
def updateBLen(nodeList,node,mutMatrix,useRateVariation=False,mutMatrices=None):
    cNode=node
    node=node.up
    if cNode==node.children[0]:
        vectUp=node.probVectUpRight
        cNum=0
    else:
        vectUp=node.probVectUpLeft
        cNum=1
    bestLength=estimateBranchLengthWithDerivative(vectUp,cNode.probVect,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices, node2isleaf=(cNode.children==[]))
    if bestLength: #cNode.dist != bestLength:
        cNode.dist=bestLength
        node.dirty=True
        cNode.dirty=True
        nodeList.append((cNode,2))
        nodeList.append((node,cNum))





# update the partials iteratively starting from the nodes in nodeList
#each entry in nodeList contains the node it refers to, and the direction where the update comes from (0 is left child, 1 is right child, 2 is parent)
def updatePartials(nodeList,mutMatrix,useRateVariation=False,mutMatrices=None):
    while nodeList:
        updatedBLen=False # if there has been an inconsistency, function updateBLen() has been called, and so there is no point continuing with some updates.
        node, direction = nodeList.pop()
        node.dirty=True
        if node.up != None:
            if node==node.up.children[0]:
                childNumUp=0
                vectUpUp=node.up.probVectUpRight
            else:
                childNumUp=1
                vectUpUp=node.up.probVectUpLeft
        #change in likelihoods is coming from parent node
        if direction==2:
            if node.dist : #if necessary, update the total probabilities at the mid node.
                newTot=mergeVectorsUpDown(vectUpUp,node.dist/2,node.probVect,node.dist/2,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices, node2isleaf=node.children==[])
                if newTot==None:
                    if node.dist>1e-100:
                        print("inside updatePartials(), from parent: should not have happened since node.dist>0")
                    updateBLen(nodeList,node,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                    updatedBLen=True
                else:
                    node.probVectTotUp=newTot
            else:  # Changed this 20/9
                if hasattr(node, 'probVectTotUp'):  # Changed this 20/9
                    node.probVectTotUp = None  # Changed this 20/9
            if len(node.children)>0 and (not updatedBLen): #at valid internal node, update upLeft and upRight, and if necessary add children to nodeList.
                child0Vect=node.children[0].probVect
                child1Vect=node.children[1].probVect
                dist0=node.children[0].dist
                dist1=node.children[1].dist
                newUpRight=mergeVectorsUpDown(vectUpUp,node.dist,child1Vect,dist1,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices, node2isleaf=node.children[1].children==[])
                if newUpRight==None:
                    if (not node.dist) and (not dist1):
                        updateBLen(nodeList,node,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                    else:
                        print("Strange: None vector from non-zero distances in updatePartials() from parent direction.")
                        raise Exception("exit")
                    updatedBLen=True
                if not updatedBLen:
                    newUpLeft=mergeVectorsUpDown(vectUpUp,node.dist,child0Vect,dist0,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,node2isleaf=node.children[0].children==[])
                    if newUpLeft==None:
                        if (not node.dist) and (not dist0) :
                            updateBLen(nodeList,node,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                        else:
                            print("Strange: None vector from non-zero distances in updatePartials() from parent direction, child0.")
                            raise Exception("exit")
                        updatedBLen=True
                if not updatedBLen:
                    if areVectorsDifferent(node.probVectUpRight,newUpRight):
                        node.probVectUpRight=newUpRight
                        nodeList.append((node.children[0],2))
                    if areVectorsDifferent(node.probVectUpLeft,newUpLeft):
                        node.probVectUpLeft=newUpLeft
                        nodeList.append((node.children[1],2))

        else: #change in likelihoods is coming from child number "direction".
            childNum=direction
            otherChildNum=1-childNum
            childDist=node.children[childNum].dist
            otherChildDist=node.children[otherChildNum].dist
            otherChildVect=node.children[otherChildNum].probVect
            probVectDown=node.children[childNum].probVect
            if childNum:
                otherVectUp=node.probVectUpRight
            else:
                otherVectUp=node.probVectUpLeft

            #update lower likelihoods
            newVect=mergeVectors(otherChildVect,otherChildDist,probVectDown,childDist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,
                                 node1isleaf=node.children[otherChildNum].children==[],
                                 node2isleaf=node.children[childNum].children==[])
            if newVect==None:
                if (not childDist) and (not otherChildDist):
                    updateBLen(nodeList,node.children[childNum],mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                    updatedBLen=True
                else:
                    print("Strange: None vector from non-zero distances in updatePartials() from child direction.")
                    print(childDist)
                    print(probVectDown)
                    print(otherChildDist)
                    print(otherChildVect)
                    raise Exception("exit")
            else:
                oldProbVect=node.probVect
                node.probVect=newVect

            #update total mid-branches likelihood
            if not updatedBLen:
                if node.dist and node.up!=None:
                    newTot=mergeVectorsUpDown(vectUpUp,node.dist/2,node.probVect,node.dist/2,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,
                                              node2isleaf=node.children==[]) #.children[childNum]
                    if newTot==None:
                        updateBLen(nodeList,node,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                        updatedBLen=True
                        print("inside updatePartials(), from child: should not have happened since node.dist>0")
                    else:
                        node.probVectTotUp=newTot
                elif node.up!=None: # Changed this 20/9
                    if hasattr(node, 'probVectTotUp'): # Changed this 20/9
                        node.probVectTotUp = None # Changed this 20/9

            if not updatedBLen:
                #update likelihoods at parent node
                if areVectorsDifferent(node.probVect,oldProbVect):
                    if node.up != None:
                        nodeList.append((node.up,childNumUp))
                    if verbose:
                        print("Old upRight and UpLeft:")
                        print(node.probVectUpRight)
                        print(node.probVectUpLeft)

                #update likelihoods at sibling node
                if node.up != None:
                    if verbose:
                        print("Trying to move update to sibling while updating partials")
                    newUpVect=mergeVectorsUpDown(vectUpUp,node.dist,probVectDown,childDist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,
                                                 node2isleaf=node.children[childNum].children==[])
                else:
                    if verbose:
                        print("reached root while moving up updating likelihoods and trying to move to sibling")
                    newUpVect=rootVector(probVectDown,childDist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,
                                         isLeaf=node.children[childNum].children==[])
                if newUpVect==None:
                    if (not node.dist) and (not childDist):
                        updateBLen(nodeList,node,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                        updatedBLen=True
                    else:
                        print("Strange: None vector from non-zero distances in updatePartials() from child direction, newUpVect.")
                        raise Exception("exit")
                else:
                    if areVectorsDifferent(otherVectUp,newUpVect):
                        if verbose:
                            print("Moving update to sibling; newUpVect and otherVectUp")
                            print(newUpVect)
                            print(otherVectUp)
                        if childNum:
                            node.probVectUpRight=newUpVect
                        else:
                            node.probVectUpLeft=newUpVect
                        nodeList.append((node.children[otherChildNum],2))
        if verbose:
            print("New down, upRight, UpLeft and tot:")
            print(node.probVect)
            print(node.probVectUpRight)
            print(node.probVectUpLeft)








#we know that sample "sample", with partials "newPartials", is best placed near a node resulting in logLK contribution of newChildLK
# now explore exactly which position near node is best for placement (direct child, sibling or sibling of existing child), and which are the best branch lengths,
#then add the sample at that position of the tree, and update all the internal probability vectors.
#UNNECESSARY? Could probably just be replaced by the more general placeSubtreeOnTree().
def placeSampleOnTree(node,newPartials,sample,newChildLK, bestUpLength, bestDownLength, bestAppendingLength,mutMatrix,pseudoMutCounts):
    tryNewRoot=False
    if node.up==None:
        tryNewRoot=True
        totRoot=rootVector(node.probVect,False,mutMatrix)
        bestAppendingLength=estimateBranchLengthWithDerivative(totRoot,newPartials,mutMatrix)
        root=node
        newChildLK=appendProb(totRoot,newPartials,bestAppendingLength,mutMatrix)
    # v0.1.9
    else:
        if node.up.children[0]==node:
            child=0
            vectUp=node.up.probVectUpRight
        else:
            child=1
            vectUp=node.up.probVectUpLeft
        if not bestUpLength:
            pNode=node.up
            while (not pNode.dist) and (pNode.up!=None):
                pNode=pNode.up
            if pNode.up==None:
                root=pNode
                tryNewRoot=True
                if (not bestDownLength) or (bestDownLength>1.01*node.dist) or (bestDownLength<0.99*node.dist):
                    node.dist=bestDownLength
                    nodeList=[(node,2),(node.up,child)]
                    updatePartials(nodeList,mutMatrix)

    #in case of best placement as a descendant appended exactly at the root node, attempt also to create new root
    if tryNewRoot:
        node=root
        probOldRoot = findProbRoot(node.probVect)
        rootUpLeft=rootVector(node.probVect,bestAppendingLength/2,mutMatrix)
        bestRightLength=estimateBranchLengthWithDerivative(rootUpLeft,newPartials,mutMatrix)
        rootUpRight=rootVector(newPartials,bestRightLength,mutMatrix)
        bestLeftLength=estimateBranchLengthWithDerivative(rootUpRight,node.probVect,mutMatrix)
        secondBranchLengthOptimizationRound=True
        if secondBranchLengthOptimizationRound: #if wanted, do a second round of branch length optimization
            rootUpLeft=rootVector(node.probVect,bestLeftLength,mutMatrix)
            bestRightLength=estimateBranchLengthWithDerivative(rootUpLeft,newPartials,mutMatrix)
            rootUpRight=rootVector(newPartials,bestRightLength,mutMatrix)
            bestLeftLength=estimateBranchLengthWithDerivative(rootUpRight,node.probVect,mutMatrix)
        probVectRoot,probRoot = mergeVectors(node.probVect,bestLeftLength,newPartials,bestRightLength,mutMatrix,returnLK=True)
        probRoot+= findProbRoot(probVectRoot)
        parentLKdiff=probRoot-probOldRoot
        if parentLKdiff<=newChildLK: #best is just placing as descendant of the root
            bestRightLength=bestAppendingLength
            bestLeftLength=False
            probVectRoot=mergeVectors(node.probVect,bestLeftLength,newPartials,bestRightLength,mutMatrix)
            rootUpRight=rootVector(newPartials,bestRightLength,mutMatrix)
        #now add new root to the tree
        newRoot=Tree()
        newRoot.probVect=probVectRoot
        newRoot.probVectUpRight=rootUpRight
        newRoot.probVectUpLeft=rootVector(node.probVect,bestLeftLength,mutMatrix)
        node.up=newRoot
        node.dist=bestLeftLength
        newRoot.add_child(node)
        newNode=Tree(name=sample,dist=bestRightLength)
        newNode.minorSequences=[]
        newNode.up=newRoot
        newRoot.add_child(newNode)
        newNode.probVect=newPartials
        if bestRightLength:
            newNode.probVectTotUp=mergeVectorsUpDown(newRoot.probVectUpLeft,bestRightLength/2,newPartials,bestRightLength/2,mutMatrix)
        if verbose:
            print("new root added to tree")
            print(newRoot.probVect)
            print(newRoot.children[0].probVect)
            print(newNode.probVect)
        #updatePartialsFromTop(node,newRoot.probVectUpRight,mutMatrix)
        nodeList=[(node,2)]
        updatePartials(nodeList,mutMatrix)
        return newRoot

    #print("adding internal node")
    #in all other cases (not attempting to add a new root) create a new internal node in the tree and add sample as a descendant.
    if node.up.children[0]==node:
        child=0
        vectUp=node.up.probVectUpRight
    else:
        child=1
        vectUp=node.up.probVectUpLeft
    newInternalNode=Tree()
    node.up.children[child]=newInternalNode
    newInternalNode.up=node.up
    newInternalNode.add_child(node)
    node.up=newInternalNode
    node.dist=bestDownLength
    newNode=Tree(name=sample,dist=bestAppendingLength)
    newNode.minorSequences=[]
    newNode.up=newInternalNode
    newInternalNode.add_child(newNode)
    newInternalNode.dist=bestUpLength
    newNode.probVect=newPartials
    newInternalNode.probVect=mergeVectors(node.probVect,bestDownLength,newPartials,bestAppendingLength,mutMatrix)
    newInternalNode.probVectUpRight=mergeVectorsUpDown(vectUp,bestUpLength,newPartials,bestAppendingLength,mutMatrix)
    newInternalNode.probVectUpLeft=mergeVectorsUpDown(vectUp,bestUpLength,node.probVect,bestDownLength,mutMatrix)
    if bestUpLength:
        newInternalNode.probVectTotUp=mergeVectorsUpDown(vectUp,bestUpLength/2,newInternalNode.probVect,bestUpLength/2,mutMatrix)
    if bestAppendingLength:
        newNode.probVectTotUp=mergeVectorsUpDown(newInternalNode.probVectUpLeft,bestAppendingLength/2,newPartials,bestAppendingLength/2,mutMatrix)
        updatePesudoCounts(newInternalNode.probVectUpLeft,newPartials,pseudoMutCounts)
    if not bestDownLength:
        node.probVectTotUp=None
    if verbose:
        print("new internal node added to tree")
        print(newInternalNode.probVect)
    nodeList=[(node,2),(newInternalNode.up,child)]
    updatePartials(nodeList,mutMatrix)

    return None




#set all descendant nodes to dirty.
#So far thses flags are used to prevent traversing the same part of the tree multiple times.
def setAllDirty(node):
    nextLeaves=[node]
    while nextLeaves:
        nextNode=nextLeaves.pop()
        nextNode.dirty=True
        for c in nextNode.children:
            nextLeaves.append(c)




#function to calculate likelihood cost of appending node to parent node
#differently from appendProb, this allows the bottom node to be internal, not just a sample.
def appendProbNode(probVectP,probVectC,bLen,mutMatrix,useRateVariation=False,mutMatrices=None, node1isleaf=False, node2isleaf=False):
    Lkcost, indexEntry1, indexEntry2, totalFactor, pos = 0.0, 0, 0, 1.0, 0
    entry1=probVectP[indexEntry1] #parent
    entry2=probVectC[indexEntry2] #child
    end=min(entry1[1],entry2[1])
    contribLength=bLen
    while True:
        if entry2[0]==5: # case entry1 is N
            pos=min(entry1[1],entry2[1])
            pass
        elif entry1[0]==5: # case entry2 is N
            #if parent node is type "N", in theory we might have to calculate the contribution of root nucleotides;
            # however, if this node is "N" then every other node in the current tree is "N", so we can ignore this since this contribution cancels out in relative terms.
            pos=min(entry1[1],entry2[1])
            pass
        else:
            #contribLength will be here the total length from the root or from the upper node, down to the down node.
            if entry1[0]<5:
                if len(entry1)==2:
                    contribLength=bLen
                elif len(entry1)==3:
                    contribLength=entry1[2]
                    if bLen:
                        contribLength+=bLen
                else:
                    contribLength=entry1[3]
                    if bLen:
                        contribLength+=bLen
            else:
                if len(entry1)==3:
                    contribLength=bLen
                else:
                    contribLength=entry1[2]
                    if bLen:
                        contribLength+=bLen
            if entry2[0]<5:
                if len(entry2)==3:
                    if contribLength:
                        contribLength+=entry2[2]
                    else:
                        contribLength=entry2[2]
            else:
                if len(entry2)==4:
                    if contribLength:
                        contribLength+=entry2[2]
                    else:
                        contribLength=entry2[2]
            #assert(contribLength == contribLength2)
            if entry1[0]==4: # case entry1 is R
                if entry2[0]==4:
                    if len(entry1)==4:
                        end=min(entry1[1],entry2[1])
                        contribLength+=entry1[2]
                        Lkcost+=contribLength*(cumulativeRate[end]-cumulativeRate[pos])
                        pos=end
                    else:
                        if contribLength:
                            end=min(entry1[1],entry2[1])
                            Lkcost+=contribLength*(cumulativeRate[end]-cumulativeRate[pos])
                            pos=end
                        else:
                            pos=min(entry1[1],entry2[1])

                #entry1 is reference and entry2 is of type "O"
                elif entry2[0]==6:
                    if useRateVariation:
                        mutMatrix=mutMatrices[pos]
                    i1=refIndeces[pos]
                    if len(entry1)==4:
                        tot=0.0
                        for i in range4:
                            if i1==i:
                                tot2=rootFreqs[i]*(1.0+mutMatrix[i][i]*entry1[2])
                            else:
                                tot2=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]
                            if contribLength:
                                tot3=0.0
                                for j in range4:
                                    tot3+=mutMatrix[i][j]*entry2[-1][j]
                                tot+=tot2*(entry2[-1][i] + contribLength*tot3)
                            else:
                                tot+=tot2*entry2[-1][i]
                        tot/=rootFreqs[i1]
                    else:
                        if contribLength:
                            tot=0.0
                            for j in range4:
                                tot+=mutMatrix[i1][j]*entry2[-1][j]
                            tot*=contribLength
                            tot+=entry2[-1][i1]
                        else:
                            tot=entry2[-1][i1]
                    totalFactor*=tot
                    pos+=1

                else: #entry1 is R and entry2 is a different but single nucleotide
                    if useRateVariation:
                        mutMatrix=mutMatrices[pos]
                    if len(entry1)==4:
                        i1=refIndeces[pos]
                        i2=entry2[0]
                        if contribLength:
                            totalFactor*=((rootFreqs[i1]*mutMatrix[i1][i2]*contribLength*(1.0+mutMatrix[i1][i1]*entry1[2])+rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]*(1.0+mutMatrix[i2][i2]*contribLength))/rootFreqs[i1])
                        else:
                            totalFactor*=((rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2])/rootFreqs[i1])
                    else:
                        if contribLength:
                            totalFactor*=mutMatrix[refIndeces[pos]][entry2[0]]*contribLength
                        else:
                            return float("-inf")
                    pos+=1

            # entry1 is of type "O"
            elif entry1[0]==6:
                if useRateVariation:
                    mutMatrix=mutMatrices[pos]
                if entry2[0]==6:
                    if contribLength:
                        tot=0.0
                        for j in range4:
                            tot+=entry1[-1][j]*(entry2[-1][j] + contribLength*(mutMatrix[j][0]*entry2[-1][0]+mutMatrix[j][1]*entry2[-1][1]+mutMatrix[j][2]*entry2[-1][2]+mutMatrix[j][3]*entry2[-1][3]) )
                        totalFactor*=tot
                    else:
                        tot=0.0
                        for j in range4:
                            tot+=entry1[-1][j]*entry2[-1][j]
                        totalFactor*=tot
                else: #entry1 is "O" and entry2 is a nucleotide
                    if entry2[0]==4:
                        i2=refIndeces[pos]
                    else:
                        i2=entry2[0]
                    if contribLength:
                        totalFactor*=(entry1[-1][i2]+contribLength*(entry1[-1][0]*mutMatrix[0][i2]+entry1[-1][1]*mutMatrix[1][i2]+entry1[-1][2]*mutMatrix[2][i2]+entry1[-1][3]*mutMatrix[3][i2]))
                    else:
                        totalFactor*=entry1[-1][i2]
                pos+=1

            else: #entry1 is a non-ref nuc
                if useRateVariation:
                    mutMatrix=mutMatrices[pos]
                if entry2[0]==entry1[0]:
                    if len(entry1)==4:
                        contribLength+=entry1[2]
                    if contribLength:
                        Lkcost+=mutMatrix[entry1[0]][entry1[0]]*contribLength
                else: #entry1 is a nucleotide and entry2 is not the same as entry1
                    i1=entry1[0]
                    if entry2[0]<5: #entry2 is a nucleotide
                        if entry2[0]==4:
                            i2=refIndeces[pos]
                        else:
                            i2=entry2[0]
                        if len(entry1)==4:
                            if contribLength:
                                totalFactor*=(( rootFreqs[i1]*mutMatrix[i1][i2]*contribLength*(1.0+mutMatrix[i1][i1]*entry1[2]) + rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]*(1.0+mutMatrix[i2][i2]*contribLength) )/rootFreqs[i1])
                            else:
                                totalFactor*=((rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2] )/rootFreqs[i1])
                        else:
                            if contribLength:
                                totalFactor*=mutMatrix[i1][i2]*contribLength
                            else:
                                return float("-inf")
                    else: #entry1 is a nucleotide and entry2 is of type "O"
                        if len(entry1)==4:
                            tot=0.0
                            for i in range4:
                                if i1==i:
                                    tot2=rootFreqs[i]*(1.0+mutMatrix[i][i]*entry1[2])
                                else:
                                    tot2=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]
                                tot3=0.0
                                for j in range4:
                                    tot3+=mutMatrix[i][j]*entry2[-1][j]
                                tot+=tot2*(entry2[-1][i] + contribLength*tot3)
                            totalFactor*=(tot/rootFreqs[i1])
                        else:
                            tot=0.0
                            for j in range4:
                                tot+=mutMatrix[i1][j]*entry2[-1][j]
                            tot*=contribLength
                            tot+=entry2[-1][i1]
                            totalFactor*=tot
                pos+=1

        if totalFactor<=minimumCarryOver:
            if totalFactor<sys.float_info.min:
                return float("-inf")
            Lkcost+=log(totalFactor)
            totalFactor=1.0
        if pos==lRef:
            break
        if pos==entry1[1]:
            indexEntry1+=1
            entry1=probVectP[indexEntry1]
        if pos==entry2[1]:
            indexEntry2+=1
            entry2=probVectC[indexEntry2]
    return Lkcost+log(totalFactor)


#calculate derivative starting from coefficients.
def calculateDerivative(ais,t):
    derivative=0.0
    for ai in ais:
        derivative+=1.0/(ai+t)
    return derivative




#function to optimize branch lengths.
#calculate features of the derivative of the likelihood cost function wrt the branch length, then finds branch length that minimizes likelihood cost.
def estimateBranchLengthWithDerivative(probVectP,probVectC,mutMatrix,useRateVariation=False,mutMatrices=None, node2isleaf=False):
    c1=0.0
    ais=[]
    indexEntry1, indexEntry2, pos = 0, 0, 0
    entry1=probVectP[indexEntry1]
    entry2=probVectC[indexEntry2]
    end=min(entry1[1],entry2[1])
    while True:
        if entry2[0]==5: # case entry1 is N
            pos=min(entry1[1],entry2[1])
            pass
        elif entry1[0]==5: # case entry2 is N
            #if parent node is type "N", in theory we might have to calculate the contribution of root nucleotides;
            # however, if this node is "N" then every other node in the current tree is "N", so we can ignore this since this contribution cancels out in relative terms.
            pos=min(entry1[1],entry2[1])
            pass
        else:
            #contribLength will be here the total length from the root or from the upper node, down to the down node.
            if entry1[0]<5:
                if len(entry1)==2:
                    contribLength=False
                elif len(entry1)==3:
                    contribLength=entry1[2]
                else:
                    contribLength=entry1[3]
            else:
                if len(entry1)==3:
                    contribLength=False
                else:
                    contribLength=entry1[2]
            if entry2[0]<5:
                if len(entry2)==3:
                    if contribLength:
                        contribLength+=entry2[2]
                    else:
                        contribLength=entry2[2]
            else:
                if len(entry2)==4:
                    if contribLength:
                        contribLength+=entry2[2]
                    else:
                        contribLength=entry2[2]

            if entry1[0]==4: # case entry1 is R
                if entry2[0]==4:
                    end=min(entry1[1],entry2[1])
                    c1+=(cumulativeRate[end]-cumulativeRate[pos])
                    pos=end

                #entry1 is reference and entry2 is of type "O"
                elif entry2[0]==6:
                    if useRateVariation:
                        mutMatrix=mutMatrices[pos]
                        #print(mutMatrix)
                        #print(pos)
                    i1=refIndeces[pos]
                    if len(entry1)==4:
                        coeff0=rootFreqs[i1]*entry2[-1][i1]
                        coeff1=0.0
                        for i in range4:
                            coeff0+=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]*entry2[-1][i]
                            coeff1+=mutMatrix[i1][i]*entry2[-1][i]
                        coeff1*=rootFreqs[i1]
                        if contribLength:
                            coeff0+=coeff1*contribLength
                    else:
                        coeff0=entry2[-1][i1]
                        coeff1=0.0
                        for j in range4:
                            coeff1+=mutMatrix[i1][j]*entry2[-1][j]
                        if contribLength:
                            coeff0+=coeff1*contribLength
                    if coeff1<0.0:
                        c1+=coeff1/coeff0
                    elif coeff1:
                        coeff0=coeff0/coeff1
                        ais.append(coeff0)
                    pos+=1

                else: #entry1 is R and entry2 is a different but single nucleotide
                    if len(entry1)==4:
                        if useRateVariation:
                            mutMatrix=mutMatrices[pos]
                        i1=refIndeces[pos]
                        i2=entry2[0]
                        if contribLength:
                            coeff0=rootFreqs[i1]*mutMatrix[i1][i2]*contribLength+rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]
                            coeff1=rootFreqs[i1]*mutMatrix[i1][i2]
                        else:
                            coeff0=rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]
                            coeff1=rootFreqs[i1]*mutMatrix[i1][i2]
                        coeff0=coeff0/coeff1
                    else:
                        if contribLength:
                            coeff0=contribLength
                        else:
                            coeff0=0.0
                    ais.append(coeff0)
                    pos+=1

            # entry1 is of type "O"
            elif entry1[0]==6:
                if useRateVariation:
                    mutMatrix=mutMatrices[pos]
                if entry2[0]==6:
                    coeff0=entry1[-1][0]*entry2[-1][0]+entry1[-1][1]*entry2[-1][1]+entry1[-1][2]*entry2[-1][2]+entry1[-1][3]*entry2[-1][3]
                    coeff1=0.0
                    for i in range4:
                        for j in range4:
                            coeff1+=entry1[-1][i]*entry2[-1][j]*mutMatrix[i][j]
                    if contribLength:
                        coeff0+=coeff1*contribLength
                else: #entry1 is "O" and entry2 is a nucleotide
                    if entry2[0]==4:
                        i2=refIndeces[pos]
                    else:
                        i2=entry2[0]
                    coeff0=entry1[-1][i2]
                    coeff1=0.0
                    for i in range4:
                        coeff1+=entry1[-1][i]*mutMatrix[i][i2]
                    if contribLength:
                        coeff0+=coeff1*contribLength
                if coeff1<0.0:
                    c1+=coeff1/coeff0
                elif coeff1:
                    coeff0=coeff0/coeff1
                    ais.append(coeff0)
                pos+=1

            else: #entry1 is a non-ref nuc
                if useRateVariation:
                    mutMatrix=mutMatrices[pos]
                if entry2[0]==entry1[0]:
                    c1+=mutMatrix[entry1[0]][entry1[0]]
                else: #entry1 is a nucleotide and entry2 is not the same as entry1
                    i1=entry1[0]
                    if entry2[0]<5: #entry2 is a nucleotide
                        if entry2[0]==4:
                            i2=refIndeces[pos]
                        else:
                            i2=entry2[0]

                        if len(entry1)==4:
                            if contribLength:
                                coeff0=rootFreqs[i1]*mutMatrix[i1][i2]*contribLength+rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]
                                coeff1=rootFreqs[i1]*mutMatrix[i1][i2]
                            else:
                                coeff0=rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]
                                coeff1=rootFreqs[i1]*mutMatrix[i1][i2]
                            coeff0=coeff0/coeff1
                        else:
                            if contribLength:
                                coeff0=contribLength
                            else:
                                coeff0=0.0
                        ais.append(coeff0)

                    else: #entry1 is a nucleotide and entry2 is of type "O"
                        if len(entry1)==4:
                            coeff0=rootFreqs[i1]*entry2[-1][i1]
                            coeff1=0.0
                            for i in range4:
                                coeff0+=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]*entry2[-1][i]
                                coeff1+=mutMatrix[i1][i]*entry2[-1][i]
                            coeff1*=rootFreqs[i1]
                            if contribLength:
                                coeff0+=coeff1*contribLength
                        else:
                            coeff0=entry2[-1][i1]
                            coeff1=0.0
                            for j in range4:
                                coeff1+=mutMatrix[i1][j]*entry2[-1][j]
                            if contribLength:
                                coeff0+=coeff1*contribLength
                        if coeff1<0.0:
                            c1+=coeff1/coeff0
                        elif coeff1:
                            coeff0=coeff0/coeff1
                            ais.append(coeff0)
                pos+=1
        #print(c1, ais, pos)
        if pos==lRef:
            break
        if pos==entry1[1]:
            indexEntry1+=1
            entry1=probVectP[indexEntry1]
        if pos==entry2[1]:
            indexEntry2+=1
            entry2=probVectC[indexEntry2]

    #now optimized branch length based on coefficients
    #minBLenSensitivity=minBLen/100000
    c1=-c1
    n=len(ais)
    if n==0:
        return False
    else:
        tDown = n / c1 - min(ais) #v0.1.9
        if tDown <= 0.0:
            return 0.0
        # vDown=calculateDerivative(ais,tDown)
        vDown = 0.0
        for ai in ais:
            vDown += 1.0 / (ai + tDown)
        tUp = n / c1 - max(ais)
        if tUp <= minBLenSensitivity:
            if min(ais):
                tUp = 0.0
            else:
                tUp = minBLenSensitivity
        vUp = 0.0
        for ai in ais:
            vUp += 1.0 / (ai + tUp)
    if vDown>c1+minBLenSensitivity or vUp<c1-minBLenSensitivity:
        if vUp<c1-minBLenSensitivity and (not tUp):
            return 0.0
        print(c1,ais)
        print("Initial values")
        print(tDown,vDown,tUp,vUp)
        print("Initial border parameters don't fit expectations:")
        print(c1,ais)
        print(probVectP)
        print(probVectC)
        print(mutMatrix)
        #raise Exception("exit")

    while tDown-tUp>minBLenSensitivity:
        tMiddle=(tUp+tDown)/2
        vMiddle=calculateDerivative(ais,tMiddle)
        #print(tUp,tMiddle,tDown,vMiddle,c1)
        if vMiddle>c1:
            tUp=tMiddle
        else:
            tDown=tMiddle
    return tUp





#traverse the tree to optimize (only) the length of the branches using the derivative approach.
def traverseTreeToOptimizeBranchLengths(root,mutMatrix,testing=False,useRateVariation=False,mutMatrices=None):
    totalLKimprovementBL=0.0
    updates=0
    nodesTraversed=0
    if root.children:
        nodesToTraverse=[root.children[0],root.children[1]]
    else:
        return 0
    while nodesToTraverse:
        nodesTraversed+=1
        node=nodesToTraverse.pop()
        if node==node.up.children[0]:
            upVect=node.up.probVectUpRight
            child=0
        else:
            upVect=node.up.probVectUpLeft
            child=1
        if node.dirty:
            bestLength=estimateBranchLengthWithDerivative(upVect,node.probVect,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices, node2isleaf=(node.children==[]))


            if bestLength or node.dist:
                if testing:
                    currentCost=appendProbNode(upVect,node.probVect,node.dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,  node2isleaf=(node.children==[]))
                    newCost=appendProbNode(upVect,node.probVect,bestLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices, node2isleaf=(node.children==[]))
                    if newCost<currentCost-thresholdLogLKconsecutivePlacement:
                        print("New branch length")
                        print(currentCost,node.dist,newCost,bestLength)
                        print("Worse length?")
                        print("Improvement LK "+str(newCost-currentCost))
                        print(upVect)
                        print(node.probVect)
                    totalLKimprovementBL+=newCost-currentCost
                if (not bestLength) or (not node.dist) or node.dist/bestLength>1.01 or node.dist/bestLength<0.99:
                    node.dist=bestLength
                    updates+=1
                    nodeList=[(node,2),(node.up,child)]
                    updatePartials(nodeList,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
        for child in node.children:
            nodesToTraverse.append(child)
    #print(nodesTraversed)
    if testing:
        return totalLKimprovementBL
    else:
        return updates


#
# upVect,node.probVect,node.dist
# ([(4, 2, 0.16470633438995388, 0.0),
#   (1, 3, 0.16470633438995388, 0.0),
#   (4, 11, 0.16470633438995388, 0.0),
#   (1, 12, 0.16470633438995388, 0.0),
#   (0, 13, 0.16470633438995388, 0.0),
#   (4, 16, 0.16470633438995388, 0.0)],
#  [(4, 2), (1, 3), (4, 11), (1, 12), (4, 14), (0, 15), (4, 16)],
#  False)






#we know that subtree "appendedNode", with partials "newPartials", is best placed as child of "node" resulting in logLK contribution of newChildLK
# now explore exactly which position near node is best for placement (direct child, sibling or sibling of existing child), and which are the best branch lengths,
#then add the subtree at that position of the tree, and update all the internal probability vectors.
def placeSubtreeOnTree(node,newPartials,appendedNode,newChildLK,bestBranchLengths,mutMatrix,useRateVariation=False,mutMatrices=None, isLeaf=False):
    #isLeaf says whether newPartials is from a leaf, which comes from appendedNode.probvect. #isLeaf = appendedNode.children ==[]
    bestAppendingLength=bestBranchLengths[2]
    bestUpLength=bestBranchLengths[0]
    bestDownLength=bestBranchLengths[1]
    tryNewRoot=False
    # if node.up==None:
    # 	tryNewRoot=True
    # 	totRoot=rootVector(node.probVect,False,mutMatrix)
    # 	bestAppendingLength=estimateBranchLengthWithDerivative(totRoot,newPartials,mutMatrix)
    # 	root=node
    # 	newChildLK=appendProbNode(totRoot,newPartials,bestAppendingLength,mutMatrix)
    #elif not bestUpLength:
    if node.up.children[0] == node:
        child = 0
        vectUp = node.up.probVectUpRight
    else:
        child = 1
        vectUp = node.up.probVectUpLeft

    if not bestUpLength:
        pNode=node.up
        while (not pNode.dist) and (pNode.up!=None):
            pNode=pNode.up
        if pNode.up==None:
            root=pNode
            tryNewRoot=True
            if (not bestDownLength) or (bestDownLength > 1.01 * node.dist) or (bestDownLength < 0.99 * node.dist):
                node.dist = bestDownLength
                nodeList = [(node, 2), (node.up, child)]
                updatePartials(nodeList, mutMatrix, useRateVariation=useRateVariation, mutMatrices=mutMatrices)
    #in case of best placement as a descendant appended exactly at the root node, attempt also to create new root
    if tryNewRoot:
        node=root
        probOldRoot = findProbRoot(node.probVect)
        rootUpLeft=rootVector(node.probVect,bestAppendingLength/2,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
        bestRightLength=estimateBranchLengthWithDerivative(rootUpLeft,newPartials,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,node2isleaf=(node.children==[]))
        rootUpRight=rootVector(newPartials,bestRightLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
        bestLeftLength=estimateBranchLengthWithDerivative(rootUpRight,node.probVect,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,node2isleaf=(node.children==[]))
        secondBranchLengthOptimizationRound=True
        if secondBranchLengthOptimizationRound: #if wanted, do a second round of branch length optimization
            rootUpLeft=rootVector(node.probVect,bestLeftLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
            bestRightLength=estimateBranchLengthWithDerivative(rootUpLeft,newPartials,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,node2isleaf=(node.children==[]))
            rootUpRight=rootVector(newPartials,bestRightLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
            bestLeftLength=estimateBranchLengthWithDerivative(rootUpRight,node.probVect,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,node2isleaf=(node.children==[]))
        probVectRoot,probRoot = mergeVectors(node.probVect,bestLeftLength,newPartials,
                                             bestRightLength,mutMatrix,returnLK=True,
                                             useRateVariation=useRateVariation,mutMatrices=mutMatrices,
                                             node1isleaf=node.children==[],
                                             node2isleaf=isLeaf)
        probRoot+= findProbRoot(probVectRoot)
        parentLKdiff=probRoot-probOldRoot
        if parentLKdiff<=newChildLK: #best is just placing as descendant of the root
            bestRightLength=bestAppendingLength
            bestLeftLength=False
            probVectRoot=mergeVectors(node.probVect,bestLeftLength,newPartials,bestRightLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,
                                             node1isleaf=node.children==[],
                                             node2isleaf=isLeaf)
            rootUpRight=rootVector(newPartials,bestRightLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices, isLeaf=isLeaf)
        #now add new root to the tree
        newRoot=Tree()
        newRoot.dirty=True
        newRoot.probVect=probVectRoot
        newRoot.probVectUpRight=rootUpRight
        newRoot.probVectUpLeft=rootVector(node.probVect,bestLeftLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices, isLeaf=node.children==[])
        node.up=newRoot
        node.dist=bestLeftLength
        newRoot.add_child(node)
        #newNode=Tree(dist=bestRightLength)
        appendedNode.up=newRoot
        newRoot.add_child(appendedNode)
        appendedNode.dist=bestRightLength
        #newNode.probVect=newPartials
        if verbose:
            print("new root added to tree")
            print(newRoot.probVect)
            print(newRoot.children[0].probVect)
            print(appendedNode.probVect)
        #updatePartialsFromTop(node,newRoot.probVectUpRight,mutMatrix)
        nodeList=[(node,2),(appendedNode,2)]
        updatePartials(nodeList,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
        return newRoot

    #print("adding internal node")
    #in all other cases (not attempting to add a new root) create a new internal node in the tree and add sample as a descendant.
    if node.up.children[0]==node:
        child=0
        vectUp=node.up.probVectUpRight
    else:
        child=1
        vectUp=node.up.probVectUpLeft
    newInternalNode=Tree()
    newInternalNode.dirty=True
    node.up.children[child]=newInternalNode
    newInternalNode.up=node.up
    newInternalNode.add_child(node)
    node.up=newInternalNode
    node.dist=bestDownLength
    #newNode=Tree(name=sample,dist=bestAppendingLength)
    #newNode.minorSequences=[]
    appendedNode.up=newInternalNode
    appendedNode.dist=bestAppendingLength
    newInternalNode.add_child(appendedNode)
    newInternalNode.dist=bestUpLength
    #newNode.probVect=newPartials
    newInternalNode.probVect=mergeVectors(node.probVect,bestDownLength,newPartials,bestAppendingLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,
                                             node1isleaf=(node.children==[]),
                                             node2isleaf=isLeaf) #isLeaf, says whether the appendedNode has children or not
    newInternalNode.probVectUpRight=mergeVectorsUpDown(vectUp,bestUpLength,newPartials,bestAppendingLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,
                                             node1isleaf=False,
                                             node2isleaf=isLeaf)
    newInternalNode.probVectUpLeft=mergeVectorsUpDown(vectUp,bestUpLength,node.probVect,bestDownLength,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,
                                             node1isleaf=False,
                                             node2isleaf=node.children==[])
    if bestUpLength:
        newInternalNode.probVectTotUp=mergeVectorsUpDown(vectUp,bestUpLength/2,newInternalNode.probVect,bestUpLength/2,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,
                                                         node1isleaf=False, node2isleaf=False)  #newinternal node is no leaf.
    if not bestDownLength:
        node.probVectTotUp=None
    if verbose:
        print("new internal node added to tree")
        print(newInternalNode.probVect)
    nodeList=[(node,2),(newInternalNode.up,child),(appendedNode,2)]
    updatePartials(nodeList,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)

    return None





#remove node from the current position in the tree and re-attach it at a new given place new bestNode.
# First remove node from the tree, then update the genome lists;
# then find the exact best reattachment of node and update the genome lists again using function placeSubtreeOnTree().

def cutAndPasteNode(node,bestNode,bestBranchLengths,bestLK,mutMatrix,useRateVariation=False,mutMatrices=None):
    #remove node from the tree
    #verbose=True
    if verbose or debugging:
        print("In cutAndPasteNode() removing subtree from the tree, subtree root partials: ")
        print(node.probVect)
        print("likelihoods to which it is attached:")
        if node==node.up.children[0]:
            print(node.up.probVectUpRight)
        else:
            print(node.up.probVectUpLeft)
    parentNode=node.up
    if node==parentNode.children[0]:
        sibling=parentNode.children[1]
    else:
        sibling=parentNode.children[0]
    if parentNode.up!=None:
        if parentNode==parentNode.up.children[0]:
            childP=0
        else:
            childP=1

    if parentNode.up!=None:
        parentNode.up.children[childP]=sibling
    sibling.up=parentNode.up
    if sibling.dist:
        if parentNode.dist:
            sibling.dist+=parentNode.dist
    else:
        sibling.dist=parentNode.dist
    #update likelihood lists after node removal
    if sibling.up==None:
        sibling.dist=1.0
        if verbose:
            print("cutAndPasteNode(), sibling node is root")
        #sibling.probVectTot=rootVector(sibling.probVect,False,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
        if sibling.children:
            sibling.probVectUpRight=rootVector(sibling.children[1].probVect,sibling.children[1].dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices, isLeaf=sibling.children[1].children==[])
            sibling.probVectUpLeft=rootVector(sibling.children[0].probVect,sibling.children[0].dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices, isLeaf=sibling.children[0].children==[])
            if debugging:
                print("sibling node is root, children distances: "+str(sibling.children[0].dist)+" "+str(sibling.children[1].dist)+" new upLeft:")
                print(sibling.probVectUpLeft)
            nodeList=[(sibling.children[0],2),(sibling.children[1],2)]
            updatePartials(nodeList,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
            if debugging:
                print("sibling node is root, after updatePartials, children distances: "+str(sibling.children[0].dist)+" "+str(sibling.children[1].dist)+" new upLeft:")
                print(sibling.probVectUpLeft)
                print(sibling.children)
            #updatePartialsFromTop(sibling.children[0],sibling.probVectUpRight,mutMatrix)
            #updatePartialsFromTop(sibling.children[1],sibling.probVectUpLeft,mutMatrix)
    else:
        nodeList=[(sibling,2),(sibling.up,childP)]
        updatePartials(nodeList,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
        if verbose:
            if sibling.dist:
                print("cutAndPasteNode(), sibling node is not root. New overall likelihoods at previously sibling node: ")
                print(sibling.probVectTot)
        #updatePartialsFromTop(sibling,vectUp,mutMatrix)
        #updatePartialsFromBottom(sibling.up,sibling.probVect,childP,sibling,mutMatrix)
    #re-place the node and re-update the vector lists
    newRoot = placeSubtreeOnTree(bestNode,node.probVect,node,bestLK,bestBranchLengths, mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices, isLeaf=node.children==[])
    # if debugging:
    # 	root=node
    # 	while root.up!=None:
    # 		root=root.up
        #print("after re-attaching subtree to the tree: "+createBinaryNewick(root))
        #print(sibling.probVectUpLeft)
        #print(sibling.children)
        #print(sibling.up); #print("And subtree "+createBinaryNewick(sibling))
    #if the root of the tree has changed, return the new root

    if sibling.up==None:
        return sibling
    else:
        return newRoot


# try to find a re-placement of a dirty node of the tree to improve the topology.
# Cut out the subtree at this node, and look for somewhere else in the tree where to attach it (an SPR move).
# To find the best location of the new re-attachment, we traverse the tree starting at the current attachment node, and for each node visited we evaluate the reattachment,
# as done by the findBestParentTopology() function.
# To avoid traversing the whole tree for each SPR move, we use stopping conditions similar to those used in the initial sample placement process.
# After we find the best SPR move for the given node, we execute it using the cutAndPasteNode() function.

def traverseTreeForTopologyUpdate(node,mutMatrix,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,thresholdTopologyPlacement=thresholdTopologyPlacement,useRateVariation=False,mutMatrices=None):
    #track if the root has changed, so that the new root node can be returned.
    newRoot=None
    # has the branch length been updated so far? If we change the branch length above the current node, even if we don't perform any SPR, we still need to update the genome lists in the tree.
    bLenChanged=False
    totalImprovement=0.0
    # we avoid the root node since it cannot be re-placed with SPR moves
    if node.up!=None:
        #evaluate current placement
        parentNode=node.up
        if parentNode.children[0]==node:
            child=0
            vectUp=parentNode.probVectUpRight
        else:
            child=1
            vectUp=parentNode.probVectUpLeft
        #score of current tree
        bestOriginalBLen=node.dist
        bestCurrenBLen=node.dist
        originalLK=appendProbNode(vectUp,node.probVect,bestCurrenBLen,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices, node2isleaf=(node.children==[]))
        bestCurrentLK=originalLK
        if bestCurrentLK<thresholdTopologyPlacement:
            bestCurrenBLen=estimateBranchLengthWithDerivative(vectUp,node.probVect,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices, node2isleaf=(node.children==[]))
            if bestCurrenBLen or node.dist:
                bestCurrentLK=appendProbNode(vectUp,node.probVect,bestCurrenBLen,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices, node2isleaf=(node.children==[]))
                if (not bestCurrenBLen) or (not node.dist) or node.dist/bestCurrenBLen>1.01 or node.dist/bestCurrenBLen<0.99:
                    bLenChanged=True; totalImprovement = bestCurrentLK - originalLK
        topologyUpdated=False
        if bestCurrentLK<thresholdTopologyPlacement:
            #now find the best place on the tree where to re-attach the subtree rooted ar "node"
            #but to do that we need to consider new vector probabilities after removing the node that we want to replace
            # this is done using findBestParentTopology().
            if debugging:
                if bLenChanged:
                    print("Better branch length found at node "+str(node)+", from "+str(node.dist)+" (LK "+str(originalLK)+") to "+str(bestCurrenBLen)+" (LK "+str(bestCurrentLK)+"). Children")
                    print(node.children)

            #bestNodeSoFar , bestLKdiff , bestIsMidNode = findBestParentTopologyOld(parentNode,child,bestCurrentLK,bestCurrenBLen,mutMatrix)
            #bestNodeSoFar , bestLKdiff , bestIsMidNode = findBestParent(parentNode,child,bestCurrentLK,bestCurrenBLen,mutMatrix,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
            bestNodeSoFar , bestLKdiff , bestBranchLengths = findBestParentTopology(parentNode,child,bestCurrentLK,bestCurrenBLen,mutMatrix,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
            # if bestBranchLengths == (None, None, None): # remove , only for debugging
            #     bestNodeSoFar , bestLKdiff , bestBranchLengths = findBestParentTopology(parentNode,child,bestCurrentLK,bestCurrenBLen,mutMatrix,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,useRateVariation=useRateVariation,mutMatrices=mutMatrices)

            # if bestLKdiff>thresholdProb2:
            # 	print("Strange, LK cost is positive")
            # 	print(bestLKdiff)
            # 	print(bestBranchLengths)
            # 	print(bestNodeSoFar.probVectTotUp)
            # 	print(node.probVect)
            # 	raise Exception("exit")
            #elif bestLKdiff<-1e50:
            if bestLKdiff<-1e50:
                print("Error: found likelihood cost is very heavy, this might mean that the reference used is not the same used to generate the input diff file: ", str(bestLKdiff))
                raise Exception("exit")
            #if debugging:
            #	print("Best SPR move found for this node has cost "+str(bestLKdiff)+" affecting node with children")
            #	print(bestNodeSoFar.children)
            #	print("With appending distance "+str(bestCurrenBLen))
            if bestLKdiff+thresholdTopologyPlacement>bestCurrentLK:
                topologyUpdated=True
                topNode=node.up
                if bestNodeSoFar==topNode:
                    topologyUpdated=False
                while (not topNode.dist) and (topNode.up!=None):
                    topNode=topNode.up
                if bestNode==topNode and (not bestBranchLengths[1]):
                    topologyUpdated=False
                parentNode=node.up
                if node==parentNode.children[0]:
                    sibling=parentNode.children[1]
                else:
                    sibling=parentNode.children[0]
                if bestNode==sibling:
                    topologyUpdated=False
                if bestNodeSoFar.up == sibling and (not bestBranchLengths[0]):
                    topologyUpdated = False

                if topologyUpdated:
                    totalImprovement=(bestLKdiff-originalLK)
                    if verbose:
                        print("\n\n In traverseTreeForTopologyUpdate() found SPR move with improvement "+str(totalImprovement))
                    if debugging:
                        print("\n\n In traverseTreeForTopologyUpdate() found SPR move with improvement "+str(totalImprovement))
                        print("Performing SPR move, detaching node "+str(node)+" with children")
                        print(node.children)
                        print("and reattaching it around node "+str(bestNodeSoFar)+" with children")
                        print(bestNodeSoFar.children)
                        root=node
                        while root.up!=None:
                            root=root.up
                        #print("Tree before cutAndPasteNode(): "+createBinaryNewick(root))
                    if hasattr(bestNodeSoFar, "name"):
                        if bestNodeSoFar.name =='S121': findBestParentTopology(parentNode, child, bestCurrentLK, bestCurrenBLen, mutMatrix, strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,    useRateVariation=useRateVariation, mutMatrices=mutMatrices)
                    newRoot = cutAndPasteNode(node,bestNodeSoFar,bestBranchLengths,bestLKdiff,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                    bLenChanged=False
        if (not topologyUpdated) and bLenChanged:
            if debugging:
                print("Changing branch length (pos 3) from "+str(node.dist)+" to "+str(bestCurrenBLen)+" at node "+str(node)+" with children ")
                print(node.children)
            node.dist=bestCurrenBLen
            # if debugging:
            # 	root=node
            # 	while root.up!=None:
            # 		root=root.up
            # 	print("New tree "+createBinaryNewick(root))
            #updatePartialsFromTop(node,vectUp,mutMatrix)
            #updatePartialsFromBottom(node.up,node.probVect,child,node,mutMatrix)
            nodeList=[(node,2),(node.up,child)]
            updatePartials(nodeList,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
    #if debugging:
    #	print("End of traverseTreeForTopologyUpdate, tree after cutAndPasteNode(): "+createBinaryNewick(newRoot))
    return newRoot,totalImprovement


#traverse the tree (here the input "node" will usually be the root), and for each dirty node ancountered, call traverseTreeForTopologyUpdate()
# to attempt an SPR move by cutting the subtree rooted at this dirty node and trying to re-append it elsewhere.
def startTopologyUpdates(node,mutMatrix,checkEachSPR=False,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,thresholdTopologyPlacement=thresholdTopologyPlacement,useRateVariation=False,mutMatrices=None):
    nodesToVisit=[node]
    totalImprovement=0.0
    newRoot=None
    numNodes=0
    while nodesToVisit:
        newNode=nodesToVisit.pop()
        for c in newNode.children:
            #if c.dirty:
            nodesToVisit.append(c)
        if newNode.dirty:
            newNode.dirty=False
            if checkEachSPR:
                root=newNode
                while root.up!=None:
                    root=root.up
                reCalculateAllGenomeLists(root,mutMatrix, checkExistingAreCorrect=True,useRateVariation=useRateVariation,mutMatrices=mutMatrices)# remove!                #print("Pre-SPR tree: "+createBinaryNewick(root))
                oldTreeLK=calculateTreeLikelihood(root,mutMatrix,checkCorrectness=True,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                #print("Pre-SPR tree likelihood: "+str(oldTreeLK))
                reCalculateAllGenomeLists(root,mutMatrix, checkExistingAreCorrect=True,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
            newRoot2,improvement=traverseTreeForTopologyUpdate(newNode,mutMatrix,strictTopologyStopRules=strictTopologyStopRules,allowedFailsTopology=allowedFailsTopology,thresholdLogLKtopology=thresholdLogLKtopology,thresholdTopologyPlacement=thresholdTopologyPlacement,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
            if checkEachSPR:
                #print(" apparent improvement "+str(improvement))
                root=newNode
                while root.up!=None:
                    root=root.up
                newTreeLK=calculateTreeLikelihood(root,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                reCalculateAllGenomeLists(root,mutMatrix, checkExistingAreCorrect=True,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                #print("Post-SPR tree likelihood: "+str(newTreeLK))
                print("In startTopologyUpdates, LK score of improvement "+str(newTreeLK)+" - "+str(oldTreeLK)+" = "+str(newTreeLK-oldTreeLK)+", was supposed to be "+str(improvement))
                if newTreeLK-oldTreeLK < improvement-0.1:#1.0: #changed
                    print("In startTopologyUpdates, LK score of improvement "+str(newTreeLK)+" - "+str(oldTreeLK)+" = "+str(newTreeLK-oldTreeLK)+" is less than what is supposed to be "+str(improvement))
                    #raise Exception("exit")
                    # In startTopologyUpdates, LK score of improvement -4029.2451376863723 - -4028.3851970278783 = -0.859940658493997 is less than what is supposd to be 1.1674419238238682, after SPR move to S121
            totalImprovement+=improvement
            if newRoot2!=None:
                newRoot=newRoot2
            numNodes+=1
            if (numNodes%1000)==0:
                print("Processed topology for "+str(numNodes)+" nodes.")
    return newRoot,totalImprovement




#create newick string of a given tree (input node is assumed to be the root).
#This function can create multifurcations, which are not welcome by some software - use createBinaryNewick() if you want only binary trees
def createNewick(node):
    nextNode=node
    stringList=[]
    direction=0
    while nextNode!=None:
        if nextNode.children:
            if direction==0:
                stringList.append("(")
                nextNode=nextNode.children[0]
            elif direction==1:
                stringList.append(",")
                nextNode=nextNode.children[1]
                direction=0
            else:
                if nextNode.dist:
                    stringList.append("):"+str(nextNode.dist))
                else:
                    stringList.append("):0.0")
                if nextNode.up!=None:
                    if nextNode.up.children[0]==nextNode:
                        direction=1
                    else:
                        direction=2
                nextNode=nextNode.up
        else:
            if len(nextNode.minorSequences)>0:
                stringList.append("("+nextNode.name+":0.0")
                for s2 in nextNode.minorSequences:
                    stringList.append(","+s2+":0.0")
                if nextNode.dist:
                    stringList.append("):"+str(nextNode.dist))
                else:
                    stringList.append("):0.0")
            else:
                if nextNode.dist:
                    stringList.append(nextNode.name+":"+str(nextNode.dist))
                else:
                    stringList.append(nextNode.name+":0.0")
            if nextNode.up!=None:
                if nextNode.up.children[0]==nextNode:
                    direction=1
                else:
                    direction=2
            nextNode=nextNode.up
    stringList.append(";")
    return "".join(stringList)


def createBinaryNewick(node):
    nextNode=node
    stringList=[]
    direction=0
    while nextNode!=None:
        if nextNode.children:
            if direction==0:
                stringList.append("(")
                nextNode=nextNode.children[0]
            elif direction==1:
                stringList.append(",")
                nextNode=nextNode.children[1]
                direction=0
            else:
                if nextNode.dist:
                    stringList.append("):"+str(nextNode.dist))
                else:
                    stringList.append("):"+str(0.0))
                if nextNode.up!=None:
                    if nextNode.up.children[0]==nextNode:
                        direction=1
                    else:
                        direction=2
                nextNode=nextNode.up
        else:
            if len(nextNode.minorSequences)>0:
                for i in nextNode.minorSequences:
                    stringList.append("(")
                stringList.append(nextNode.name+":")
                for s2 in nextNode.minorSequences:
                    stringList.append(str(0.0)+","+s2+":"+str(0.0)+"):")
                if nextNode.dist:
                    stringList.append(str(nextNode.dist))
                else:
                    stringList.append(str(0.0))
            else:
                if nextNode.dist:
                    stringList.append(nextNode.name+":"+str(nextNode.dist))
                else:
                    stringList.append(nextNode.name+":"+str(0.0))
            if nextNode.up!=None:
                if nextNode.up.children[0]==nextNode:
                    direction=1
                else:
                    direction=2
            nextNode=nextNode.up
    stringList.append(";")
    return "".join(stringList)






#Given a tree, and a final substitution rate matrix, calculate the likelihood of the tree
def calculateTreeLikelihood(root,mutMatrix,checkCorrectness=False,useRateVariation=False,mutMatrices=None):
    node=root
    #direction 0 means from parent, direction 1 means from a child
    lastNode=None
    direction=0
    totalLK=0.0
    while node!=None:
        if direction==0:
            if node.children:
                node=node.children[0]
            else:
                lastNode=node
                node=node.up
                direction=1
                #we can ignore the likelihood contribution from the normalization of the likelihoods of the ambiguity characters at terminal nodes:
                # these are anyway the same for all trees and substitution models!
        else :
            if lastNode==node.children[0]:
                node=node.children[1]
                direction=0
            else:
                newLower, Lkcontribution=mergeVectors(node.children[0].probVect,node.children[0].dist,
                                                      node.children[1].probVect,node.children[1].dist,mutMatrix,
                                                      returnLK=True,useRateVariation=useRateVariation,mutMatrices=mutMatrices,
                                                      node1isleaf=node.children[0].children==[], node2isleaf=node.children[1].children==[])

                totalLK+=Lkcontribution
                if newLower==None:
                    print("Strange, inconsistent lower genome list creation in calculateTreeLikelihood(); "
                          "old list, and children lists")
                    print(node.probVect)
                    print(node.children[0].probVect)
                    print(node.children[1].probVect)
                    raise Exception("exit")
                elif numTopologyImprovements and checkCorrectness and areVectorsDifferentDebugging(node.probVect,newLower):
                    print("Strange, while calculating tree likelihood encountered non-updated lower likelihood genome list at node "+str(node)+" with children ")
                    print(node.children)
                    print("new list")
                    print(newLower)
                    print("old list")
                    print(node.probVect)
                    print("child 0 list")
                    print(node.children[0].probVect)
                    print("child 1 list")
                    print(node.children[1].probVect)
                    raise Exception("exit")
                lastNode=node
                node=node.up
                direction=1
    #now add contribution from the root
    totalLK+=findProbRoot(root.probVect)
    # print(findProbRoot(root.probVect))
    return totalLK





#Given a tree and its genome lists, calculate mutation counts and waiting times for expectation maximization estimation of substitution rates
def expectationMaximizationCalculationRates(root,useRateVariation=False):
    node=root
    #direction 0 means from parent, direction 1 means from a child
    lastNode=None
    direction=0
    counts=[[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0]]
    waitingTimes=[0.0,0.0,0.0,0.0]
    if useRateVariation:
        totTreeLength=0.0
        waitingTimesSites=[]
        countsSites=[]
        #list used to track missed waiting times across the genome due to Ns
        trackingNs=[]
        for i in range(lRef):
            waitingTimesSites.append([0.0,0.0,0.0,0.0])
            countsSites.append(0.0)
            trackingNs.append(0.0)
        trackingNs.append(0.0)

    while node!=None:
        if direction==0:
            if node.dist and node.up!=None:
                if useRateVariation:
                    totTreeLength+=node.dist
                #update counts and waiting times
                if node==node.up.children[0]:
                    probVectP=node.up.probVectUpRight
                else:
                    probVectP=node.up.probVectUpLeft
                probVectC=node.probVect
                indexEntry1, indexEntry2, pos = 0, 0, 0
                entry1=probVectP[indexEntry1]
                entry2=probVectC[indexEntry2]
                end=min(entry1[1],entry2[1])
                contribLength=node.dist

                while True:
                    if entry2[0]==5: # case N
                        if useRateVariation:
                            trackingNs[pos]-=node.dist
                        pos=min(entry1[1],entry2[1])
                        if useRateVariation:
                            trackingNs[pos]+=node.dist
                    elif entry1[0]==5: # case N
                        if useRateVariation:
                            trackingNs[pos]-=node.dist
                        pos=min(entry1[1],entry2[1])
                        if useRateVariation:
                            trackingNs[pos]+=node.dist
                    else:
                        if entry1[0]<5:
                            if len(entry1)==2:
                                totLen1=node.dist
                            elif len(entry1)==3:
                                totLen1=entry1[2]
                                if node.dist:
                                    totLen1+=node.dist
                            else:
                                totLen1=entry1[3]
                                if node.dist:
                                    totLen1+=node.dist
                                #consider contribution across the root twice, each time only considering the bit on the own side of the root,
                        else:
                            if len(entry1)==3:
                                totLen1=node.dist
                            else:
                                totLen1=entry1[2]
                                if node.dist:
                                    totLen1+=node.dist

                        if entry2[0]<5:
                            if len(entry2)==2:
                                totLen2=False
                            else:
                                totLen2=entry2[2]
                        else:
                            if len(entry2)==3:
                                totLen2=False
                            else:
                                totLen2=entry2[2]
                        if entry1[0]==4: # case entry1 is R
                            if entry2[0]==4:
                                end=min(entry1[1],entry2[1])
                                if not totLen2:
                                    for i in range4:
                                        waitingTimes[i]+=totLen1*(cumulativeBases[end][i]-cumulativeBases[pos][i])
                                pos=end

                            elif entry2[0]==6:
                                if not totLen2:
                                    i1=refIndeces[pos]
                                    normalization=0.0
                                    if len(entry1)==4:
                                        contribLength=node.dist+entry1[3]
                                        if useRateVariation:
                                            waitingTimesSites[pos][i1]-=contribLength
                                        for i in range4:
                                            if i1==i:
                                                prob=rootFreqs[i]*(1.0+nonMutRates[i]*entry1[2])
                                                tot3=0.0
                                                for j in range4:
                                                    tot3+=mutMatrix[i][j]*entry2[-1][j]
                                                tot3*=contribLength
                                                tot3+=entry2[-1][i]
                                                normalization+=prob*tot3
                                            else:
                                                prob=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]*(1.0+nonMutRates[i]*contribLength)*entry2[-1][i]
                                                normalization+=prob
                                        for i in range4:
                                            if i1==i:
                                                prob=rootFreqs[i]*(1.0+nonMutRates[i]*entry1[2])
                                                for j in range4:
                                                    if j==i:
                                                        tot3=prob*(1.0+nonMutRates[i]*contribLength)*entry2[-1][j]
                                                        waitingTimes[i]+=contribLength*tot3/normalization
                                                        if useRateVariation:
                                                            waitingTimesSites[pos][i]+=contribLength*tot3/normalization
                                                    else:
                                                        tot3=prob*mutMatrix[i][j]*contribLength*entry2[-1][j]/normalization
                                                        waitingTimes[i]+=(contribLength/2)*tot3
                                                        waitingTimes[j]+=(contribLength/2)*tot3
                                                        counts[i][j]+=tot3
                                                        if useRateVariation:
                                                            waitingTimesSites[pos][i]+=(contribLength/2)*tot3
                                                            waitingTimesSites[pos][j]+=(contribLength/2)*tot3
                                                            countsSites[pos]+=tot3
                                            else:
                                                prob=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]*(1.0+nonMutRates[i]*contribLength)*entry2[-1][i]
                                                waitingTimes[i]+=contribLength*prob/normalization
                                                if useRateVariation:
                                                    waitingTimesSites[pos][i]+=contribLength*prob/normalization

                                    else:
                                        if useRateVariation:
                                            waitingTimesSites[pos][i1]-=totLen1
                                        for i in range4:
                                            if i1==i:
                                                normalization+=(1.0+nonMutRates[i]*totLen1)*entry2[-1][i]
                                            else:
                                                normalization+=mutMatrix[i1][i]*totLen1*entry2[-1][i]
                                        for i in range4:
                                            if i1==i:
                                                prob=(1.0+nonMutRates[i]*totLen1)*entry2[-1][i]
                                                waitingTimes[i]+=totLen1*prob/normalization
                                                if useRateVariation:
                                                    waitingTimesSites[pos][i]+=totLen1*prob/normalization
                                            else:
                                                prob=mutMatrix[i1][i]*totLen1*entry2[-1][i]/normalization
                                                waitingTimes[i1]+=(totLen1/2)*prob
                                                waitingTimes[i]+=(totLen1/2)*prob
                                                counts[i1][i]+=prob
                                                if useRateVariation:
                                                    waitingTimesSites[pos][i1]+=(totLen1/2)*prob
                                                    waitingTimesSites[pos][i]+=(totLen1/2)*prob
                                                    countsSites[pos]+=prob
                                pos+=1

                            else: #entry1 is reference and entry2 is a different but single nucleotide
                                if not totLen2:
                                    i1=refIndeces[pos]
                                    i2=entry2[0]
                                    if len(entry1)<4:
                                        if useRateVariation:
                                            waitingTimesSites[pos][i1]-=totLen1/2
                                            waitingTimesSites[pos][i2]+=totLen1/2
                                            countsSites[pos]+=1
                                        waitingTimes[i1]+=(totLen1/2)
                                        waitingTimes[i2]+=(totLen1/2)
                                        counts[i1][i2]+=1
                                    else:
                                        contribLength=node.dist+entry1[3]
                                        prob1=rootFreqs[i1]*mutMatrix[i1][i2]*contribLength*(1.0+nonMutRates[i1]*entry1[2])
                                        prob2=rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]*(1.0+nonMutRates[i2]*contribLength)
                                        normalization=prob1+prob2
                                        prob1=prob1/normalization
                                        prob2=prob2/normalization
                                        waitingTimes[i1]+=(contribLength/2)*prob1
                                        waitingTimes[i2]+=(contribLength/2)*prob1
                                        counts[i1][i2]+=prob1
                                        waitingTimes[i2]+=contribLength*prob2
                                        if useRateVariation:
                                            waitingTimesSites[pos][i1]-=contribLength
                                            waitingTimesSites[pos][i1]+=(contribLength/2)*prob1
                                            waitingTimesSites[pos][i2]+=(contribLength/2)*prob1
                                            waitingTimesSites[pos][i2]+=contribLength*prob2
                                            countsSites[pos]+=prob1
                                pos+=1

                        # entry1 is of type "O"
                        elif entry1[0]==6:
                            if not totLen2:
                                #bLen13=totLen1
                                normalization=0.0
                                if useRateVariation:
                                    waitingTimesSites[pos][refIndeces[pos]]-=totLen1
                                if entry2[0]==6:
                                    for i in range4:
                                        for j in range4:
                                            if i==j:
                                                normalization+=entry1[-1][i]*(1.0+nonMutRates[i]*totLen1)*entry2[-1][j]
                                            else:
                                                normalization+=entry1[-1][i]*mutMatrix[i][j]*totLen1*entry2[-1][j]
                                    for i in range4:
                                        for j in range4:
                                            if i==j:
                                                prob=entry1[-1][i]*(1.0+nonMutRates[i]*totLen1)*entry2[-1][j]
                                                waitingTimes[i]+=totLen1*prob/normalization
                                                if useRateVariation:
                                                    waitingTimesSites[pos][i]+=totLen1*prob/normalization
                                            else:
                                                prob=entry1[-1][i]*mutMatrix[i][j]*totLen1*entry2[-1][j]/normalization
                                                waitingTimes[i]+=(totLen1/2)*prob
                                                waitingTimes[j]+=(totLen1/2)*prob
                                                counts[i][j]+=prob
                                                if useRateVariation:
                                                    waitingTimesSites[pos][i]+=(totLen1/2)*prob
                                                    waitingTimesSites[pos][j]+=(totLen1/2)*prob
                                                    countsSites[pos]+=prob
                                else:
                                    if entry2[0]==4:
                                        i2=refIndeces[pos]
                                    else:
                                        i2=entry2[0]
                                    for i in range4:
                                        if i==i2:
                                            normalization+=entry1[-1][i]*(1.0+nonMutRates[i]*totLen1)
                                        else:
                                            normalization+=entry1[-1][i]*mutMatrix[i][i2]*totLen1
                                    for i in range4:
                                        if i==i2:
                                            prob=entry1[-1][i]*(1.0+nonMutRates[i]*totLen1)
                                            waitingTimes[i]+=totLen1*prob/normalization
                                            if useRateVariation:
                                                waitingTimesSites[pos][i]+=totLen1*prob/normalization
                                        else:
                                            prob=entry1[-1][i]*mutMatrix[i][i2]*totLen1/normalization
                                            waitingTimes[i]+=(totLen1/2)*prob
                                            waitingTimes[i2]+=(totLen1/2)*prob
                                            counts[i][i2]+=prob
                                            if useRateVariation:
                                                waitingTimesSites[pos][i]+=(totLen1/2)*prob
                                                waitingTimesSites[pos][i2]+=(totLen1/2)*prob
                                                countsSites[pos]+=prob
                            pos+=1

                        else: #entry1 is a non-ref nuc
                            i1=entry1[0]
                            if entry2[0]==i1:
                                if not totLen2:
                                    waitingTimes[i1]+=totLen1
                                    if useRateVariation:
                                        waitingTimesSites[pos][i1]+=totLen1
                                        waitingTimesSites[pos][refIndeces[pos]]-=totLen1
                            else: #entry1 and entry2 are of different types
                                if entry2[0]==6:
                                    if not totLen2:
                                        normalization=0.0
                                        if len(entry1)==4:
                                            contribLength=node.dist+entry1[3]
                                            if useRateVariation:
                                                waitingTimesSites[pos][refIndeces[pos]]-=contribLength
                                            for i in range4:
                                                if i1==i:
                                                    prob=rootFreqs[i]*(1.0+nonMutRates[i]*entry1[2])
                                                    tot3=0.0
                                                    for j in range4:
                                                        tot3+=mutMatrix[i][j]*entry2[-1][j]
                                                    tot3*=contribLength
                                                    tot3+=entry2[-1][i]
                                                    normalization+=prob*tot3
                                                else:
                                                    prob=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]*(1.0+nonMutRates[i]*contribLength)*entry2[-1][i]
                                                    normalization+=prob
                                            for i in range4:
                                                if i1==i:
                                                    prob=rootFreqs[i]*(1.0+nonMutRates[i]*entry1[2])
                                                    for j in range4:
                                                        if j==i:
                                                            tot3=prob*(1.0+nonMutRates[i]*contribLength)*entry2[-1][j]
                                                            waitingTimes[i]+=contribLength*tot3/normalization
                                                            if useRateVariation:
                                                                waitingTimesSites[pos][i]+=contribLength*tot3/normalization
                                                        else:
                                                            tot3=prob*mutMatrix[i][j]*contribLength*entry2[-1][j]/normalization
                                                            waitingTimes[i]+=(contribLength/2)*tot3
                                                            waitingTimes[j]+=(contribLength/2)*tot3
                                                            counts[i][j]+=tot3
                                                            if useRateVariation:
                                                                waitingTimesSites[pos][i]+=(contribLength/2)*tot3
                                                                waitingTimesSites[pos][j]+=(contribLength/2)*tot3
                                                                countsSites[pos]+=tot3
                                                else:
                                                    prob=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]*(1.0+nonMutRates[i]*contribLength)*entry2[-1][i]
                                                    waitingTimes[i]+=contribLength*prob/normalization
                                                    if useRateVariation:
                                                        waitingTimesSites[pos][i]+=contribLength*prob/normalization
                                        else:
                                            if useRateVariation:
                                                waitingTimesSites[pos][refIndeces[pos]]-=totLen1
                                            for i in range4:
                                                if i1==i:
                                                    normalization+=(1.0+nonMutRates[i]*totLen1)*entry2[-1][i]
                                                else:
                                                    normalization+=mutMatrix[i1][i]*totLen1*entry2[-1][i]
                                            for i in range4:
                                                if i1==i:
                                                    prob=(1.0+nonMutRates[i]*totLen1)*entry2[-1][i]
                                                    waitingTimes[i]+=totLen1*prob/normalization
                                                    if useRateVariation:
                                                        waitingTimesSites[pos][i]+=totLen1*prob/normalization
                                                else:
                                                    prob=mutMatrix[i1][i]*totLen1*entry2[-1][i]/normalization
                                                    waitingTimes[i1]+=(totLen1/2)*prob
                                                    waitingTimes[i]+=(totLen1/2)*prob
                                                    counts[i1][i]+=prob
                                                    if useRateVariation:
                                                        waitingTimesSites[pos][i1]+=(totLen1/2)*prob
                                                        waitingTimesSites[pos][i]+=(totLen1/2)*prob
                                                        countsSites[pos]+=prob
                                #entry2 is a nucleotide type (like entry1)
                                else:
                                    if not totLen2:
                                        if entry2[0]==4:
                                            i2=refIndeces[pos]
                                        else:
                                            i2=entry2[0]
                                        if len(entry1)<4:
                                            if useRateVariation:
                                                waitingTimesSites[pos][refIndeces[pos]]-=totLen1
                                                waitingTimesSites[pos][i1]+=(totLen1/2)
                                                waitingTimesSites[pos][i2]+=(totLen1/2)
                                                countsSites[pos]+=1
                                            waitingTimes[i1]+=(totLen1/2)
                                            waitingTimes[i2]+=(totLen1/2)
                                            counts[i1][i2]+=1
                                        else:
                                            contribLength=node.dist+entry1[3]
                                            prob1=rootFreqs[i1]*mutMatrix[i1][i2]*contribLength*(1.0+nonMutRates[i1]*entry1[2])
                                            prob2=rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]*(1.0+nonMutRates[i2]*contribLength)
                                            normalization=prob1+prob2
                                            prob1=prob1/normalization
                                            prob2=prob2/normalization
                                            waitingTimes[i1]+=(contribLength/2)*prob1
                                            waitingTimes[i2]+=(contribLength/2)*prob1
                                            counts[i1][i2]+=prob1
                                            waitingTimes[i2]+=contribLength*prob2
                                            if useRateVariation:
                                                waitingTimesSites[pos][refIndeces[pos]]-=contribLength
                                                waitingTimesSites[pos][i1]+=(contribLength/2)*prob1
                                                waitingTimesSites[pos][i2]+=(contribLength/2)*prob1
                                                countsSites[pos]+=prob1
                                                waitingTimesSites[pos][i2]+=contribLength*prob2

                            pos+=1

                    if pos==lRef:
                        break
                    if pos==entry1[1]:
                        indexEntry1+=1
                        entry1=probVectP[indexEntry1]
                    if pos==entry2[1]:
                        indexEntry2+=1
                        entry2=probVectC[indexEntry2]

            if node.children:
                node=node.children[0]
            else:
                lastNode=node
                node=node.up
                direction=1

        else :
            if lastNode==node.children[0]:
                node=node.children[1]
                direction=0
            else:
                lastNode=node
                node=node.up
                direction=1

    if model=="UNREST":
        for i in range4:
            if not waitingTimes[i]:
                for j in range4:
                    counts[i][j]=0.0
            else:
                for j in range4:
                    if i!=j:
                        counts[i][j]/=waitingTimes[i]
                counts[i][i]=-sum(counts[i])
    elif model=="GTR":
        newRates=[[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0]]
        for i in range4:
            if not waitingTimes[i]:
                for j in range4:
                    newRates[i][j]=0.0
            else:
                for j in range4:
                    if i!=j:
                        newRates[i][j]=(counts[i][j]+counts[j][i])/waitingTimes[i]
                newRates[i][i]=-sum(newRates[i])
        counts=newRates
    else:
        print("Expectation Maximization for given model "+model+" not implemented yet")
        raise Exception("exit")
    totRate=-(rootFreqs[0]*counts[0][0]+ rootFreqs[1]*counts[1][1]+ rootFreqs[2]*counts[2][2]+ rootFreqs[3]*counts[3][3] )
    if totRate:
        for i in range4:
            for j in range4:
                counts[i][j]=counts[i][j]/totRate

    if useRateVariation:
        siteRates=[]
        totRate=0.0
        for i in range(lRef):
            waitingTimesSites[i][refIndeces[i]]+=totTreeLength+trackingNs[i]
            totExpected=0.0
            for j in range4:
                totExpected-=waitingTimesSites[i][j]*counts[j][j]
            if not totExpected:
                siteRates.append(1.0)
            else:
                siteRates.append(countsSites[i]/totExpected)
            totRate+=siteRates[-1]
        totRate=totRate/lRef
        for i in range(lRef):
            siteRates[i]=siteRates[i]/totRate
        #print(siteRates)
        #print(totRate)

        return counts, siteRates
    else:
        return counts, None





startTime = time()
if runOnlyExample:
    data={}
    data["refSeq"]=[]
    data["Seq1"]=[('t',1),('-',10,1)]
    data["Seq1bis"]=[('t',1),('-',11,1)]
    data["Seq2"]=[('t',1),('a',3)]
    data["Seq3"]=[('t',1),('a',3),('t',5)]
    data["Seq4"]=[('a',2)]
    data["Seq5"]=[('a',2),('t',4)]
    data["Seq6"]=[('a',2),('t',6)]
    print(data)
    print(ref[:10])
    samples=data.keys()


if inputTree=="":
    print("distances with no sample names")
    distances=distancesFromRefPunishNs(data)
else:
    print("distances with sample names")
    distances=distancesFromRefPunishNs(data,samples=data.keys())
print("Distances from the reference calculated")

if inputTree=="": #initialize tree to just the initial root sample
    #extract root genome among those closest to the reference but not empty
    firstSample=distances.pop()
    t1=Tree(name=firstSample[1])
    #t1.probVect=probVectTerminalNode(data.pop(firstSample[1]))
    t1.probVect=probVectTerminalNode(data[firstSample[1]])
    data[firstSample[1]]=None
    t1.probVectTot=rootVector(t1.probVect,False,mutMatrix)
    t1.minorSequences=[]
else:
    t1=tree1

timeFinding=0.0
timePlacing=0.0
numSamples=0
while distances:
    d=distances.pop()
    numSamples+=1
    sample=d[1]
    newPartials=probVectTerminalNode(data[sample])
    data[sample]=None
    if (numSamples%updateSubstMatrixEveryThisSamples)==0:
        if model!="JC":
            if updateSubMatrix(pseudoMutCounts,model,mutMatrix):
                for i in range4:
                    nonMutRates[i]=mutMatrix[i][i]
                for i in range(lRef):
                    cumulativeRate[i+1]=cumulativeRate[i]+nonMutRates[refIndeces[i]]
    if (numSamples%1000)==0:
        print("Sample num "+str(numSamples))
    start=time()
    bestNode , bestScore, bestBranchLengths = findBestParentForNewSample(t1,newPartials,sample,mutMatrix)
    # node , bestNewLK, isMidNode, bestUpLK, bestDownLK, bestDownNode, adjustBLen=findBestParentOld(t1,newPartials,sample,mutMatrix)
    # if bestNewLK>bestScore+0.5:
    # 	print("Suboptimal placement "+str(bestScore)+" vs "+str(bestNewLK))
    # 	bestNode , bestScore, bestBranchLengths, nodesVisited = findBestParentForNewSample(t1,newPartials,sample,mutMatrix,debug=True)
    # 	node , bestNewLK, isMidNode, bestUpLK, bestDownLK, bestDownNode, adjustBLen=findBestParentOld(t1,newPartials,sample,mutMatrix,debug=True)
    # 	raise Exception("exit")
    timeFinding+=(time()-start)
    if bestBranchLengths != None: #v0.1.9
        start=time()
        newRoot=placeSampleOnTree(bestNode,newPartials,sample,bestScore, bestBranchLengths[0], bestBranchLengths[1], bestBranchLengths[2],mutMatrix,pseudoMutCounts)
        if newRoot!=None:
            t1=newRoot
        timePlacing+=(time()-start)

#print("number of plausible branches per placed sample:")
#print(float(sum(numBestNodes))/len(numBestNodes))

if runOnlyExample:
    print("Tree after initial placement:")
    newickString=createNewick(t1)
    print(newickString)
    print("Now making change to the tree to create imperfection. Moving")
    nodeToReplace=t1.children[1].children[0].children[0].children[1]
    destination=t1.children[0].children[1].children[1]
    newickString=createNewick(nodeToReplace)
    print(newickString)
    print("to")
    newickString=createNewick(destination)
    print(newickString)
    cutAndPasteNode(nodeToReplace,destination,False,0.0001,-100.0,mutMatrix)
    newickString=createNewick(t1)
    print(newickString)

data.clear()


#put sample names in the tree
if debugging and inputTree=="":
    nextLeaves=[t1]
    while nextLeaves:
        node=nextLeaves.pop()
        if not node.children:
            node.name="S"+str(node.name)
            print(node.name)
            print(node.probVect)
            for m in range(len(node.minorSequences)):
                node.minorSequences[m]="S"+str(node.minorSequences[m])
        else:
            for c in node.children:
                nextLeaves.append(c)


#recalculate all genome lists according to the final substitution model
if inputTree=="" or largeUpdate or rateVariation:
    start=time()
    reCalculateAllGenomeLists(t1,mutMatrix,countNodes=True)
    if model!="JC" or rateVariation:
        mutMatrix , siteRates = expectationMaximizationCalculationRates(t1,useRateVariation=rateVariation)
        if rateVariation:
            for i in range4:
                nonMutRates[i]=mutMatrix[i][i]
            for i in range(lRef):
                cumulativeRate[i+1]=cumulativeRate[i]+nonMutRates[refIndeces[i]]*siteRates[i]
            mutMatrices=[]
            for i in range(lRef):
                mutMatrices.append([])
                for j in range4:
                    mutMatrices[i].append(list(mutMatrix[j]))
                    for k in range4:
                        mutMatrices[i][j][k]*=siteRates[i]
        else:
            for i in range4:
                nonMutRates[i]=mutMatrix[i][i]
            for i in range(lRef):
                cumulativeRate[i+1]=cumulativeRate[i]+nonMutRates[refIndeces[i]]
            mutMatrices=None
        reCalculateAllGenomeLists(t1,mutMatrix,countNodes=False,useRateVariation=rateVariation,mutMatrices=mutMatrices)
    timeRecalculation=time()-start
    print("Time to recalculate all genome lists: "+str(timeRecalculation))

    print("Number of nodes: "+str(numNodes[0]))
    print("Os per node: "+str(float(numNodes[4])/numNodes[0]))
    print("Nucs per node: "+str(float(numNodes[1])/numNodes[0]))
    print("Ns per node: "+str(float(numNodes[3])/numNodes[0]))
    #print("R per node: "+str(float(numNodes[2])/numNodes[0]))
    #print("Non-O per node: "+str(float(numNodes[1]+numNodes[2]+numNodes[3])/numNodes[0]))







#----------------------------------------------------------
"""
ERROR RATE FUNCTIONS
Below are my error rate functions, that fall in three main categories: 
1) Functions, specifically used when the error rate is incorporated. 
2) Functions from the original MAPLE code
   that I have adapted, e.g. in terms of the length of the entries that is checked for, and additional calculations. 
3) Functions to test the error rate functions

Other functions outside of this scope have also undergone some small changes to pass on whether a genomelist arrived from a leaf or not, 
which is crucial to extract the flags for each of its entries.  Most of these changes can be found by searching for 'isLeaf' or for 'errorRate'. 

Later we should also write error functions for appendProb, for appending single nodes to the tree to enable online inference with error rates.
"""

def getErrorRatesSiteSpecific(file):
    fileI = open(file, "r")
    line = fileI.readline()
    fileI.close()
    return [float(errorRate) for errorRate in line.split(", ")] #CodeImprovement: check that these values are not nan and that the length equals lRef

if errorRateSiteSpecific:
        #check if file exists
    if not os.path.isfile(errorRateSiteSpecific):
        print("errorRateSiteSpecific file " + errorRateSiteSpecific + " not found, quitting MAPLE tree inference. Use option --input to specify a valid input file.")
        raise Exception("exit")
    errorRates = getErrorRatesSiteSpecific(file=errorRateSiteSpecific)
    cumulativeErrorRate = [0]*(len(errorRates)+1) #creating cumulative error rates.
    cumulativeErrorRate[0] = errorRates[0]
    for i in range(0,len(errorRates)):
        cumulativeErrorRate[i+1] = cumulativeErrorRate[i] + errorRates[i] # CodeImprovement: to make estimations more accurate, one could scale these log(1-errorRate) here, which would be of little extra computational demand.

def getPartialVec(i12, flag, totLen, errorRate, upNode=False):
    """
    When merging a CATGR-entry with an O entry, or for the case that T1=/=T2, it is useful to
    first create the partial vectors of the entries seperately, instead of working with one 'newVec'.
    With this function we create the partial likelihood vector for the considered entry:
    (if needed) we account for the error rates, and (if needed) we propagate the vector along its branch.
    :param i12: i1 or i2;  the nucleotide that was stored in the original entry
    :param flag: whether we should take into account the error rate
    :param totLen: branch length along which the vector should be 'mutated'
    :param upNode: should be True if the nucleotide was observed above the branch.
    :return: the new likelihood vector for this entry.
    """
    if flag:
        partialVec1 = [flag * errorRate / 3] * 4  # without error rate [0.0, 0.0, 0.0, 0.0]
        partialVec1[i12] = 1.0 - errorRate * flag  #Approximation 1, without error rate: 1.0
        if totLen:
            mutatedPartialVec1 = [0]*4 #Approximation 1
            for j in range4:
                tot = 0.0
                for i in range4:
                    if upNode:
                        tot += mutMatrix[i][j] * partialVec1[i]
                    else:
                        tot += mutMatrix[j][i] * partialVec1[i]
                tot *= totLen
                tot += partialVec1[j]
                mutatedPartialVec1[j] = tot
            partialVec1 = mutatedPartialVec1
    else:  # no flag 1
        if totLen:
            partialVec1 = []
            for i in range4:
                if i == i12:
                    partialVec1.append(1.0 + mutMatrix[i][i] * totLen)
                else:
                    if upNode:
                        partialVec1.append(mutMatrix[i12][i] * totLen)
                    else:
                        partialVec1.append(mutMatrix[i][i12] * totLen)
        else:
            partialVec1 = [0.0, 0.0, 0.0, 0.0]
            partialVec1[i12] = 1.0
    return partialVec1


def addErrorTerminalNode(node):
    """
    This function should only be run the first time the errors are incorporated in the tree.
    It takes into account the error rates in the partial likelihood vectors of O-type entries, .
    which arrive from ambiguity characters other than N, and we should derive the new probVector using the errorRate.
    For the other entry types, nothing needs to be done: Since these are leaf nodes, no flags need to be set.

    :param node: a leaf node
    """
    # When no site specific errors are used, this function could be sped up by
    # simply changing the ambiguities list to incorporate error rates right away. When the genome list entry has a pointer to the amiguity list,
    # it will be updated inmediately.
    for entry in node.probVect: # loop over the lower likelihood entries in the genome list
        #for position specific errors one could decide to check if errorRate[i] / 3 > thresholdProb: in such case the error rate doesnt have to be considered for the position
        # for position specific errors; check if errorRate[i] > thresholdHighError: store as O vector.
        if entry[0] == 6: # we are dealing with ambiguity characters other than N, and we should derive the new probVector using the errorRate.. We should not encounter this if onlyNambiguities==True:
            sumUnnormalizedVector =  sum([bool(item) for item in entry[-1]])
            if errorRateSiteSpecific: errorRate = errorRates[entry[1]] #errorRates[pos], or should I use errorRates[pos-1]?
            if sumUnnormalizedVector ==2:
                for i in range4:             # M, instead of [0.5, 0.5, 0, 0] we will now get [0.5- ⅓ε, 0.5- ⅓ε,  ⅓ε,  ⅓ε]
                    if entry[-1][i]==0:
                        entry[-1][i] = errorRate/3
                    else: #if entry[-1][i]==0.5:
                        entry[-1][i] -= errorRate/3
            elif sumUnnormalizedVector == 3:
                for i in range4:             # for V instead of [⅓, ⅓, ⅓, 0] we get, [ ⅓ - ε/9,  ⅓ -ε/9,  ⅓ -ε/9,  ⅓ε]
                    if entry[-1][i] == 0:
                        entry[-1][i] = errorRate/3
                    else:  # if entry[-1][i]==⅓:
                        entry[-1][i] -= errorRate/9 # ⅓ - ε/9
            print(entry[-1])

def getFlag(entry, isLeaf=False):
    """
    Given an entry, and whether this entry is from a genomeList that arrived from a leaf node,
    return whether the flag is True or False
    """
    if entry[0] >= 5:  # cases of O or N
        flag = False
    elif len(entry) >= 4: #length 4 or 5 and no O or N entry
        flag = entry[-1]
    elif isLeaf: #length 2 and a leaf
        flag = True
    else:  # other cases of length 2, which should have a 0 branch length element:
        flag = False
    assert (type(flag) == bool)
    return flag


def findProbRootError(probVect):
    """
    The errorrate version of 'findProbRootError'.
    """
    global errorRate
    logLK=0.0
    logFactor=1.0
    pos=0
    for entry in probVect:
        flag = getFlag(entry, isLeaf=False)
        if errorRateSiteSpecific: errorRate = errorRates[pos]
        if entry[0]==4:
            for i in range4:
                logLK+=rootFreqsLog[i]*(cumulativeBases[entry[1]][i]-cumulativeBases[pos][i]) - errorRate*flag*(entry[1] - pos)
        elif entry[0]<4:
            logLK+=rootFreqsLog[entry[0]] - errorRate*flag #K += log(𝝿(τ)*(1-e)) = log(𝝿(τ)) + log(1-e)) ~=log(𝝿(τ) -e)
        elif entry[0]==6:
            tot=0.0
            for i in range4:
                tot+=rootFreqs[i]*entry[-1][i]
            logFactor*=tot
        pos=entry[1]
    logLK+=log(logFactor)
    return logLK



def appendProbNodeErrorRate(probVectP,probVectC,bLen,mutMatrix,useRateVariation=False,mutMatrices=None, node2isleaf=False):
    """
    The errorrate version of 'appendProbNode'.
    :param node2isleaf: True if probVectC arrives from a leaf node. (probVectP can never arrive from a leaf node)
    """
    global errorRate
    Lkcost, indexEntry1, indexEntry2, totalFactor, pos = 0.0, 0, 0, 1.0, 0
    entry1=probVectP[indexEntry1]
    entry2=probVectC[indexEntry2]
    end=min(entry1[1],entry2[1])
    while True:
        if entry2[0]==5: # case entry1 is N
            pos=min(entry1[1],entry2[1])
            pass
        elif entry1[0]==5: # case entry2 is N
            #if parent node is type "N", in theory we might have to calculate the contribution of root nucleotides;
            # however, if this node is "N" then every other node in the current tree is "N", so we can ignore this since this contribution cancels out in relative terms.
            pos=min(entry1[1],entry2[1])
            pass
        else:
            if entry1[0]<5:
                if len(entry1)==2:
                    contribLength=bLen
                elif len(entry1)==4: #errorrate ==4 instead of ==3
                    contribLength=entry1[2]
                    if bLen:
                        contribLength+=bLen
                else: #len==5
                    contribLength=entry1[3]
                    if bLen:
                        contribLength+=bLen
            else:#error rate: O /N entry stays the same.
                if len(entry1)==3:
                    contribLength=bLen
                else:
                    contribLength=entry1[2]
                    if bLen:
                        contribLength+=bLen
            if entry2[0]<5:
                if len(entry2)==4: #ErrorRate: =4 instead of ==3
                    if contribLength:
                        contribLength+=entry2[2]
                    else:
                        contribLength=entry2[2]
            else:#error rate: O /N entry stays the same.
                if len(entry2)==4:
                    if contribLength:
                        contribLength+=entry2[2]
                    else:
                        contribLength=entry2[2]
            flag1, flag2 = getFlag(entry1, isLeaf=False), getFlag(entry2, isLeaf=node2isleaf) #node1isleaf must be false bcs it is the parent node.
            if errorRateSiteSpecific: errorRate = errorRates[pos]
            if entry1[0]==4: # case entry1 is R
                if entry2[0]==4: # case entry2 is R
                    if len(entry1)==5: # == 5 instead of ==4 with errorRate
                        end=min(entry1[1],entry2[1])
                        contribLength+=entry1[2]
                        if flag1 + flag2:
                            if errorRateSiteSpecific: cumErrorRate  = cumulativeErrorRate[end]-cumulativeErrorRate[pos]
                            else: cumErrorRate = - errorRate*(end-pos)               #CodeImprovement: it may be faster to calculate a cumulative errorrate array also for non-site speficic error rates, as this would reduce the number of if statements and get rid of the multiplication.
                            Lkcost += cumErrorRate*(flag1+flag2)
                        Lkcost+=contribLength*(cumulativeRate[end]-cumulativeRate[pos])
                        pos=end
                    else:
                        end=min(entry1[1],entry2[1])
                        if flag1 + flag2:
                            if errorRateSiteSpecific: cumErrorRate  = cumulativeErrorRate[end]-cumulativeErrorRate[pos]
                            else: cumErrorRate = -errorRate * (end - pos)  # Approximation 2
                            Lkcost += cumErrorRate * (flag1 + flag2) # also if not contribLength, will the error rate be taken into account.
                        if contribLength:
                            Lkcost += contribLength * (cumulativeRate[end] - cumulativeRate[pos])
                        pos = end

                #Case τ2= O: entry1 is reference and entry2 is of type "O"
                elif entry2[0]==6: #flag2 will be false
                    if useRateVariation:
                        mutMatrix=mutMatrices[pos]
                    i1=refIndeces[pos]
                    if len(entry1)==5: #==5 instead ==4 ErrorRate:
                        tot=0.0
                        for i in range4:
                            if i1==i:
                                tot2 = rootFreqs[i] * (1.0 + mutMatrix[i][i] * entry1[2] - flag1*errorRate) #Approximation 2
                            else:
                                tot2=rootFreqs[i]*(mutMatrix[i][i1]*entry1[2] + flag1*errorRate/3) #Approximation 2
                            if contribLength:
                                tot3=0.0
                                for j in range4:
                                    tot3+=mutMatrix[i][j]*entry2[-1][j]
                                tot+=tot2*(entry2[-1][i] + contribLength*tot3)
                            else:
                                tot+=tot2*entry2[-1][i]
                        tot/=rootFreqs[i1]
                    else: #entry1 will have flag1==False. only for entries accross the root, can it inherit from a child node.
                        if contribLength:
                            tot=0.0
                            for j in range4:
                                tot+=mutMatrix[i1][j]*entry2[-1][j]
                            tot*=contribLength
                            tot+=entry2[-1][i1]
                        else:
                            tot=entry2[-1][i1]
                    totalFactor*=tot
                    pos+=1

                else: #entry1 is R and entry2 is a different but single nucleotide
                    if useRateVariation:
                        mutMatrix=mutMatrices[pos]
                    if len(entry1)==5: #errorrate ==5 instead of ==4:
                        i1=refIndeces[pos]
                        i2=entry2[0]
                        # CodeImprovement: To speed this up, one could add a seperate formula for cases where contriblength =0. For such cases, the calculations would be less involved.
                        totalFactor*=(((mutMatrix[i1][i2]*contribLength+errorRate/3*flag2)* #prob that mutation (or error) occured on entry2's side of the root
                                       (1.0+mutMatrix[i1][i1]*entry1[2]-errorRate*flag1)+ #MULTIPLIED by prob that NO mutation (or error) occured on entry1's side of the root..
                                       rootFreqs[i2]/rootFreqs[i1]*(mutMatrix[i2][i1]*entry1[2]+errorRate/3*flag1)* #PLUS  the prob that mutation (or error)  occured on entry1's side of the root
                                       (1.0+mutMatrix[i2][i2]*contribLength-errorRate*flag2)) #MULTIPLIED by the prob that NO mutation (or error) occured on entry2's side of the root
                                      )
                    else: #only flag2 can be true.
                        if contribLength or flag2:
                            totalFactor*=(mutMatrix[refIndeces[pos]][entry2[0]]*contribLength + flag2*errorRate/3) #error rate. Flag1 will be false here.
                        else:
                            return float("-inf")
                    pos+=1
            # entry1 is of type "O"
            elif entry1[0]==6:
                if useRateVariation:
                    mutMatrix=mutMatrices[pos]
                if entry2[0]==6: # entry1 and entry2 is of type "O"; no Errorrate
                    if contribLength:
                        tot=0.0
                        for j in range4:
                            tot+=entry1[-1][j]*(entry2[-1][j] + contribLength*(mutMatrix[j][0]*entry2[-1][0]+mutMatrix[j][1]*entry2[-1][1]+mutMatrix[j][2]*entry2[-1][2]+mutMatrix[j][3]*entry2[-1][3]) )
                        totalFactor*=tot
                    else:
                        tot=0.0
                        for j in range4:
                            tot+=entry1[-1][j]*entry2[-1][j]
                        totalFactor*=tot
                else: #entry1 is "O" and entry2 is a nucleotide
                    if entry2[0]==4:
                        i2=refIndeces[pos]
                    else:
                        i2=entry2[0]
                    if contribLength or flag2:
                        tot2=0
                        for i in range4:
                            tot2 += entry1[-1][i]*mutMatrix[i][i2]
                        tot1_2 = (entry1[-1][i2]*(1 - 4/3*errorRate * flag2) + 1/3*errorRate*flag2 + contribLength*tot2 )
                        totalFactor *= tot1_2
                    else:
                        totalFactor*=entry1[-1][i2]
                pos+=1
            else: #entry1 is a non-ref nuc
                if useRateVariation:
                    mutMatrix=mutMatrices[pos]
                if entry2[0]==entry1[0]: #entry1 = entry2
                    if len(entry1)==5: #errorrate ==5 instead of ==4
                        contribLength+=entry1[2]
                    if contribLength or (flag1 +flag2):
                        Lkcost+=mutMatrix[entry1[0]][entry1[0]]*contribLength + (flag1 +flag2)*log(1-errorRate) #OR: if error rate is very low, - errorRate(flag1 +flag2)
                else: #entry1 is a nucleotide and entry2 is not the same as entry1
                    i1=entry1[0]
                    if entry2[0]<5: #entry2 is a nucleotide
                        if entry2[0]==4:
                            i2=refIndeces[pos]
                        else:
                            i2=entry2[0]
                        if len(entry1)==5: #errorrate ==5 instead of ==4
                            totalFactor*=(((mutMatrix[i1][i2]*contribLength+errorRate/3*flag2)* #prob that mutation (or error) occured on entry2's side of the root
                                           (1.0+mutMatrix[i1][i1]*entry1[2]-errorRate*flag1)+ #MULTIPLIED by prob that NO mutation (or error) occured on entry1's side of the root..
                                           rootFreqs[i2]/rootFreqs[i1]*(mutMatrix[i2][i1]*entry1[2]+errorRate/3*flag1)* #PLUS  the prob that mutation (or error)  occured on entry1's side of the root
                                           (1.0+mutMatrix[i2][i2]*contribLength-errorRate*flag2)) #MULTIPLIED by the prob that NO mutation (or error) occured on entry2's side of the root
                                          )
                        else:
                            if contribLength or flag2:
                                totalFactor*=(mutMatrix[i1][i2]*contribLength + flag2*errorRate/3) #error rate. Flag1 will be false
                            else:
                                return float("-inf")
                    else: #entry1 is a nucleotide and entry2 is of type "O"
                        if len(entry1)==5: #errorrate ==5 instead of ==4
                            tot=0.0
                            for i in range4:
                                if i1==i:
                                    tot2=rootFreqs[i]*(1.0+mutMatrix[i][i]*entry1[2]-errorRate*flag1)
                                else:
                                    tot2=rootFreqs[i]*(mutMatrix[i][i1]*entry1[2] + flag1*errorRate/3)
                                tot3=0.0
                                for j in range4:
                                    tot3+=mutMatrix[i][j]*entry2[-1][j]
                                tot+=tot2*(entry2[-1][i] + contribLength*tot3)
                            totalFactor*=(tot/rootFreqs[i1])
                        else: # no error rate.
                            tot=0.0
                            for j in range4:
                                tot+=mutMatrix[i1][j]*entry2[-1][j]
                            tot*=contribLength
                            tot+=entry2[-1][i1]
                            totalFactor*=tot
                pos+=1
        if totalFactor<=minimumCarryOver:
            if totalFactor<sys.float_info.min:
                return float("-inf")
            Lkcost+=log(totalFactor)
            totalFactor=1.0
        if pos==lRef:
            break
        if pos==entry1[1]:
            indexEntry1+=1
            entry1=probVectP[indexEntry1]
        if pos==entry2[1]:
            indexEntry2+=1
            entry2=probVectC[indexEntry2]
    return Lkcost+log(totalFactor)


def mergeVectorsUpDownErrorDetection(probVect1,bLenUp,probVect2,bLenDown,mutMatrix,useRateVariation=False,mutMatrices=None, node1isleaf=False, node2isleaf=False):
    """
    This function is used to aid error detection, to calculate the overall likelihood vectors of leaf nodes.
    """
    global errorRate
    indexEntry1, indexEntry2, pos = 0, 0, 0
    probVect=[]
    entry1=probVect1[indexEntry1]
    entry2=probVect2[indexEntry2]
    while True:
        flag1, flag2 = getFlag(entry1, isLeaf=node1isleaf), getFlag(entry2, isLeaf=node2isleaf)
        if errorRateSiteSpecific: errorRate = errorRates[pos]
        if entry1[0]==5:
            if entry2[0]==5:
                pos=min(entry1[1],entry2[1])
                probVect.append((5,pos))

            elif entry2[0]<5:
                pos=min(entry1[1],entry2[1])
                if len(entry2)==4: #ErrorRate (entry2 can be of length 4 or 2)
                    if bLenDown:
                        probVect.append((entry2[0],pos,entry2[2]+bLenDown,0.0, flag2)) #ErrorRate
                    else:
                        probVect.append((entry2[0],pos,entry2[2],0.0, flag2)) #ErrorRate
                else: #len=2
                    if bLenDown or flag2:
                        probVect.append((entry2[0],pos,bLenDown,0.0, flag2)) #ErrorRate
                    else:
                        probVect.append((entry2[0],pos))
            else: # entry2 case "O", entry 1 is "N"; no flags
                if useRateVariation:
                    mutMatrix=mutMatrices[pos]
                pos+=1
                if len(entry2)==4:
                    totBLen=entry2[2]
                    if bLenDown:
                        totBLen+=bLenDown
                else:
                    totBLen=bLenDown
                newVect=[]
                if totBLen:
                    for i in range4:
                        tot=0.0
                        for j in range4:
                            tot+=mutMatrix[i][j]*entry2[-1][j]
                        tot*=totBLen
                        tot+=entry2[-1][i]
                        newVect.append(tot*rootFreqs[i])
                    totSum=sum(newVect)
                    for i in range4:
                        newVect[i]/=totSum
                    probVect.append((6,pos,newVect))
                else:
                    for i in range4:
                        newVect.append(entry2[-1][i]*rootFreqs[i])
                    totSum=sum(newVect)
                    for i in range4:
                        newVect[i]/=totSum
                    probVect.append((6,pos,newVect))
        elif entry2[0]==5: #entry2 is N
            if entry1[0]<5:
                pos=min(entry1[1],entry2[1])
                if len(entry1)==2:
                    if bLenUp or flag1:
                        probVect.append((entry1[0],pos,bLenUp, flag1))
                    else:
                        probVect.append((entry1[0],pos))
                elif len(entry1)==4: #errorRate: ==4
                    if bLenUp:
                        probVect.append((entry1[0],pos,entry1[2]+bLenUp, flag1))
                    else:
                        probVect.append((entry1[0],pos,entry1[2], flag1))
                else:
                    if bLenUp:
                        probVect.append((entry1[0],pos,entry1[2],entry1[3]+bLenUp, flag1))
                    else:
                        probVect.append((entry1[0],pos,entry1[2],entry1[3], flag1))

            else: #entry1 is "O", entry 2 is "N"; no error rate.
                if useRateVariation:
                    mutMatrix=mutMatrices[pos]
                pos+=1
                if len(entry1)==4:
                    totBLen=entry1[2]
                    if bLenUp:
                        totBLen+=bLenUp
                elif bLenUp:
                    totBLen=bLenUp
                else:
                    totBLen=False
                newVect=[]
                if totBLen:
                    for i in range4:
                        tot=0.0
                        for j in range4:
                            tot+=entry1[-1][j]*mutMatrix[j][i]
                        tot*=totBLen
                        tot+=entry1[-1][i]
                        newVect.append(tot)
                    totSum=sum(newVect)
                    for i in range4:
                        newVect[i]/=totSum
                    probVect.append((6,pos,newVect))
                else:
                    probVect.append((6,pos,entry1[-1]))
        elif entry2[0]==entry1[0] and entry1[0]<5:
            pos=min(entry1[1],entry2[1])
            probVect.append((entry2[0],pos))
        else: #cases where the new genome list entry will likely be of type "O"
            if entry1[0]<5:
                if len(entry1)==2:
                    totLen1=bLenUp
                else:
                    totLen1=entry1[2]
                    if bLenUp:
                        totLen1+=bLenUp
                    if len(entry1)==4:
                        totLen1+=entry1[3]
            else:
                if len(entry1)==3:
                    totLen1=bLenUp
                else:
                    totLen1=entry1[2]
                    if bLenUp:
                        totLen1+=bLenUp

            if entry2[0]<5:
                if len(entry2)==2:
                    totLen2=bLenDown
                else:
                    totLen2=entry2[2]
                    if bLenDown:
                        totLen2+=bLenDown
            else:
                if len(entry2)==3:
                    totLen2=bLenDown
                else:
                    totLen2=entry2[2]
                    if bLenDown:
                        totLen2+=bLenDown
            if entry1[0]<5: #entry 1 and 2 are different nucleotides
                if useRateVariation:
                    mutMatrix=mutMatrices[pos]
                if entry1[0]==4:
                    i1=refIndeces[pos]
                else:
                    i1=entry1[0]
                newVec=[]
                if len(entry1)==5: #I could also integrate as part of getpartial vec?.
                    rootVec=list(rootFreqs)
                    for i in range4:
                        if i==i1:
                            rootVec[i]*=(1.0+mutMatrix[i1][i1]*entry1[2]-errorRate*flag1)
                        else:
                            rootVec[i]*=(mutMatrix[i][i1]*entry1[2] + errorRate*flag1)
                    if bLenUp:
                        lenToRoot=entry1[3]+bLenUp
                    else:
                        lenToRoot=entry1[3]
                    for j in range4:
                        tot=0.0
                        for i in range4:
                            tot+=mutMatrix[i][j]*rootVec[i]
                        tot*=lenToRoot
                        tot+=rootVec[j]
                        newVec.append(tot)
                else:
                    newVec = getPartialVec(i1, flag1, totLen1, errorRate, upNode=True)

                if entry2[0]==6: #entry 2 is "O" and entry1 is a nucleotide
                    for j in range4:
                        tot=0.0
                        if totLen2:
                            for i in range4:
                                tot+=mutMatrix[j][i]*entry2[-1][i]
                            tot*=totLen2
                        tot+=entry2[-1][j]
                        newVec[j]*=tot
                    sumV=sum(newVec)
                    for i in range4:
                        newVec[i]=newVec[i]/sumV
                    state =simplfy(newVec,refIndeces[pos])
                    pos+=1
                    if state==6:
                        probVect.append((6,pos,newVec))
                    else:
                        probVect.append((state,pos))
                else: #entry 1 and entry 2 are different nucleotides (possibly reference)
                    if entry2[0]==4:
                        i2=refIndeces[pos]
                    else:
                        i2=entry2[0]
                    partialVec2 = getPartialVec(i2, flag2, totLen2, errorRate, upNode=False)
                    #when flag2 is true and a Blen is added, the probabilities of observing other nucleotides than i2 get way smaller.
                    for i in range4:
                        newVec[i] *= partialVec2[i]
                    sumV=sum(newVec)
                    if not sumV:
                        print("situation")
                        print(entry1)
                        print(entry2)
                        print(bLenUp)
                        print(bLenDown)
                        print(totLen1)
                        print(totLen2)
                        print(newVec)
                    for i in range4:
                        newVec[i]=newVec[i]/sumV
                    pos+=1
                    probVect.append((6,pos,newVec))
            else: #entry1[0]==6:
                if useRateVariation:
                    mutMatrix=mutMatrices[pos]
                if totLen1:
                    newVec=[]
                    for i in range4:
                        tot=0.0
                        for j in range4:
                            tot+=mutMatrix[j][i]*entry1[-1][j]
                        tot*=totLen1
                        tot+=entry1[-1][i]
                        newVec.append(tot)
                else:
                    newVec=list(entry1[-1])

                if entry2[0]==6:
                    if totLen2:
                        for i in range4:
                            tot=0.0
                            for j in range4:
                                tot+=mutMatrix[i][j]*entry2[-1][j]
                            tot*=totLen2
                            tot+=entry2[-1][i]
                            newVec[i]*=tot
                    else:
                        for i in range4:
                            newVec[i]*=entry2[-1][i]
                else: #entry1 is "O" and entry2 is a nucleotide
                    if entry2[0]==4:
                        i2=refIndeces[pos]
                    else:
                        i2=entry2[0]

                    partialVec2 = getPartialVec(i2, flag2, totLen2, errorRate, upNode=False)
                    for i in range4:
                        newVec[i] *= partialVec2[i]
                sumV=sum(newVec)
                if not sumV:
                    print("mergeVectorsUpDown() returning None 3")
                    print(entry1)
                    print(entry2)
                    print(sumV)
                    print(newVec)
                    return None
                for i in range4:
                    newVec[i]=newVec[i]/sumV
                state =simplfy(newVec,refIndeces[pos])
                pos+=1
                if state==6:
                    probVect.append((6,pos,newVec))
                else:
                    probVect.append((state,pos)) #no flags needed

        if pos==lRef:
            break
        if pos==entry1[1]:
            indexEntry1+=1
            entry1=probVect1[indexEntry1]
        if pos==entry2[1]:
            indexEntry2+=1
            entry2=probVect2[indexEntry2]

    if verbose:
        print("Merged up-down ")
        print(probVect)
    #check if the final  probVect can be simplified by merging consecutive entries
    shorten(probVect)
    if verbose:
        print("Shortened up-down merging ")
        print(probVect)
    return probVect


def mergeVectorsUpDownError(probVect1,bLenUp,probVect2,bLenDown,mutMatrix,useRateVariation=False,mutMatrices=None, node1isleaf=False, node2isleaf=False):
    """
    The errorrate version of 'mergeVectorsUpDownError'.
    :param node2isleaf: True if probVect2 arrives from a leaf node. (probVect1 can never arrive from a leaf node)
    """
    global errorRate
    indexEntry1, indexEntry2, pos = 0, 0, 0
    probVect=[]
    entry1=probVect1[indexEntry1]
    entry2=probVect2[indexEntry2]
    while True:
        flag1, flag2 = getFlag(entry1, isLeaf=False), getFlag(entry2, isLeaf=node2isleaf)
        if errorRateSiteSpecific: errorRate = errorRates[pos]
        if entry1[0]==5:
            if entry2[0]==5:
                pos=min(entry1[1],entry2[1])
                probVect.append((5,pos))

            elif entry2[0]<5:
                pos=min(entry1[1],entry2[1])
                if len(entry2)==4: #ErrorRate (entry2 can be of length 4 or 2)
                    if bLenDown:
                        probVect.append((entry2[0],pos,entry2[2]+bLenDown,0.0, flag2)) #ErrorRate
                    else:
                        probVect.append((entry2[0],pos,entry2[2],0.0, flag2)) #ErrorRate
                else: #len=2
                    if bLenDown or flag2:
                        probVect.append((entry2[0],pos,bLenDown,0.0, flag2)) #ErrorRate
                    else:
                        probVect.append((entry2[0],pos))
            else: # entry2 case "O", entry 1 is "N"; no flags
                if useRateVariation:
                    mutMatrix=mutMatrices[pos]
                pos+=1
                if len(entry2)==4:
                    totBLen=entry2[2]
                    if bLenDown:
                        totBLen+=bLenDown
                else:
                    totBLen=bLenDown
                newVect=[]
                if totBLen:
                    for i in range4:
                        tot=0.0
                        for j in range4:
                            tot+=mutMatrix[i][j]*entry2[-1][j]
                        tot*=totBLen
                        tot+=entry2[-1][i]
                        newVect.append(tot*rootFreqs[i])
                    totSum=sum(newVect)
                    for i in range4:
                        newVect[i]/=totSum
                    probVect.append((6,pos,newVect))
                else:
                    for i in range4:
                        newVect.append(entry2[-1][i]*rootFreqs[i])
                    totSum=sum(newVect)
                    for i in range4:
                        newVect[i]/=totSum
                    probVect.append((6,pos,newVect))
        elif entry2[0]==5: #entry2 is N
            if entry1[0]<5:
                pos=min(entry1[1],entry2[1])
                if len(entry1)==2:
                    if bLenUp or flag1:
                        probVect.append((entry1[0],pos,bLenUp, flag1))
                    else:
                        probVect.append((entry1[0],pos))
                elif len(entry1)==4: #errorRate: ==4
                    if bLenUp:
                        probVect.append((entry1[0],pos,entry1[2]+bLenUp, flag1))
                    else:
                        probVect.append((entry1[0],pos,entry1[2], flag1))
                else:
                    if bLenUp:
                        probVect.append((entry1[0],pos,entry1[2],entry1[3]+bLenUp, flag1))
                    else:
                        probVect.append((entry1[0],pos,entry1[2],entry1[3], flag1))

            else: #entry1 is "O", entry 2 is "N"; no error rate.
                if useRateVariation:
                    mutMatrix=mutMatrices[pos]
                pos+=1
                if len(entry1)==4:
                    totBLen=entry1[2]
                    if bLenUp:
                        totBLen+=bLenUp
                elif bLenUp:
                    totBLen=bLenUp
                else:
                    totBLen=False
                newVect=[]
                if totBLen:
                    for i in range4:
                        tot=0.0
                        for j in range4:
                            tot+=entry1[-1][j]*mutMatrix[j][i]
                        tot*=totBLen
                        tot+=entry1[-1][i]
                        newVect.append(tot)
                    totSum=sum(newVect)
                    for i in range4:
                        newVect[i]/=totSum
                    probVect.append((6,pos,newVect))
                else:
                    probVect.append((6,pos,entry1[-1]))
        elif entry2[0]==entry1[0] and entry1[0]<5:
            pos=min(entry1[1],entry2[1])
            probVect.append((entry2[0],pos))
        else: #cases where the new genome list entry will likely be of type "O"
            if entry1[0]<5:
                if len(entry1)==2:
                    totLen1=bLenUp
                else:
                    totLen1=entry1[2]
                    if bLenUp:
                        totLen1+=bLenUp
                    if len(entry1)==4:
                        totLen1+=entry1[3]
            else:
                if len(entry1)==3:
                    totLen1=bLenUp
                else:
                    totLen1=entry1[2]
                    if bLenUp:
                        totLen1+=bLenUp

            if entry2[0]<5:
                if len(entry2)==2:
                    totLen2=bLenDown
                else:
                    totLen2=entry2[2]
                    if bLenDown:
                        totLen2+=bLenDown
            else:
                if len(entry2)==3:
                    totLen2=bLenDown
                else:
                    totLen2=entry2[2]
                    if bLenDown:
                        totLen2+=bLenDown
            if entry2[0]<5 and (not totLen2) and (not flag2):  #due to 0 distance, the entry will be of same type as entry2
                if (not totLen1) and entry1[0]<5:
                    #print("mergeVectorsUpDown() returning None 1")
                    #print(entry1)
                    #print(entry2)
                    #print(totLen1)
                    #print(totLen2)
                    return None
                pos=min(entry1[1],entry2[1])
                probVect.append((entry2[0],pos))
            elif entry1[0]<5 and (not totLen1): #due to 0 distance, the entry will be of same type as entry1
                pos=min(entry1[1],entry2[1])
                probVect.append((entry1[0],pos))
            elif entry1[0]<5: #entry 1 and 2 are different nucleotides
                if useRateVariation:
                    mutMatrix=mutMatrices[pos]
                if entry1[0]==4:
                    i1=refIndeces[pos]
                else:
                    i1=entry1[0]
                newVec=[]
                if len(entry1)==5: #I could also integrate as part of getpartial vec?.
                    rootVec=list(rootFreqs)
                    for i in range4:
                        if i==i1:
                            rootVec[i]*=(1.0+mutMatrix[i1][i1]*entry1[2]-errorRate*flag1)
                        else:
                            rootVec[i]*=(mutMatrix[i][i1]*entry1[2] + errorRate*flag1)
                    if bLenUp:
                        lenToRoot=entry1[3]+bLenUp
                    else:
                        lenToRoot=entry1[3]
                    for j in range4:
                        tot=0.0
                        for i in range4:
                            tot+=mutMatrix[i][j]*rootVec[i]
                        tot*=lenToRoot
                        tot+=rootVec[j]
                        newVec.append(tot)
                else:
                    newVec = getPartialVec(i1, flag1, totLen1, errorRate, upNode=True) # check again if this is correct?
                if entry2[0]==6: #entry 2 is "O" and entry1 is a nucleotide
                    for j in range4:
                        tot=0.0
                        if totLen2:
                            for i in range4:
                                tot+=mutMatrix[j][i]*entry2[-1][i]
                            tot*=totLen2
                        tot+=entry2[-1][j]
                        newVec[j]*=tot
                    sumV=sum(newVec)
                    for i in range4:
                        newVec[i]=newVec[i]/sumV
                    state =simplfy(newVec,refIndeces[pos])
                    pos+=1
                    if state==6:
                        probVect.append((6,pos,newVec))
                    else:
                        probVect.append((state,pos))
                else: #entry 1 and entry 2 are different nucleotides (possibly reference)
                    if entry2[0]==4:
                        i2=refIndeces[pos]
                    else:
                        i2=entry2[0]
                    partialVec2 = getPartialVec(i2, flag2, totLen2, errorRate, upNode=False)
                    #when flag2 is true and a Blen is added, the probabilities of observing other nucleotides than i2 get way smaller.
                    for i in range4:
                        newVec[i] *= partialVec2[i]
                    sumV=sum(newVec)
                    if not sumV:
                        print("situation")
                        print(entry1)
                        print(entry2)
                        print(bLenUp)
                        print(bLenDown)
                        print(totLen1)
                        print(totLen2)
                        print(newVec)
                    for i in range4:
                        newVec[i]=newVec[i]/sumV
                    pos+=1
                    probVect.append((6,pos,newVec))
            else: #entry1[0]==6:
                if useRateVariation:
                    mutMatrix=mutMatrices[pos]
                if totLen1:
                    newVec=[]
                    for i in range4:
                        tot=0.0
                        for j in range4:
                            tot+=mutMatrix[j][i]*entry1[-1][j]
                        tot*=totLen1
                        tot+=entry1[-1][i]
                        newVec.append(tot)
                else:
                    newVec=list(entry1[-1])

                if entry2[0]==6:
                    if totLen2:
                        for i in range4:
                            tot=0.0
                            for j in range4:
                                tot+=mutMatrix[i][j]*entry2[-1][j]
                            tot*=totLen2
                            tot+=entry2[-1][i]
                            newVec[i]*=tot
                    else:
                        for i in range4:
                            newVec[i]*=entry2[-1][i]
                else: #entry1 is "O" and entry2 is a nucleotide
                    if entry2[0]==4:
                        i2=refIndeces[pos]
                    else:
                        i2=entry2[0]
                    partialVec2 = getPartialVec(i2, flag2, totLen2, errorRate, upNode=False)
                    for i in range4:
                        newVec[i] *= partialVec2[i]
                sumV=sum(newVec)
                if not sumV:
                    print("mergeVectorsUpDown() returning None 3")
                    print(entry1)
                    print(entry2)
                    print(sumV)
                    print(newVec)
                    return None
                for i in range4:
                    newVec[i]=newVec[i]/sumV
                state =simplfy(newVec,refIndeces[pos])
                pos+=1
                if state==6:
                    probVect.append((6,pos,newVec))
                else:
                    probVect.append((state,pos)) #no flags needed

        if pos==lRef:
            break
        if pos==entry1[1]:
            indexEntry1+=1
            entry1=probVect1[indexEntry1]
        if pos==entry2[1]:
            indexEntry2+=1
            entry2=probVect2[indexEntry2]

    if verbose:
        print("Merged up-down ")
        print(probVect)
    #check if the final  probVect can be simplified by merging consecutive entries
    shorten(probVect)
    if verbose:
        print("Shortened up-down merging ")
        print(probVect)
    return probVect



def mergeVectorsError(probVect1, bLen1, probVect2, bLen2, mutMatrix, returnLK=False, useRateVariation=False,
                 mutMatrices=None, node1isleaf=False, node2isleaf=False):
    """
    The errorrate version of 'mergeVectors'.
    :param node1isleaf: True if probVect1 arrives from a leaf node
    :param node2isleaf: True if probVect2 arrives from a leaf node
    """
    global errorRate
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
                    if bLen2 or node2isleaf: #To account for cases with 0 branch length elements, when the flag needs to be set to True.
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
                    if bLen1 or node1isleaf: #To account for cases with 0 branch length elements, when the flag needs to be set to True.
                        probVect.append((entry1[0], pos, bLen1, node1isleaf)) #flag = flag1 = node1isleaf in this case.
                    else: #0-branch length: no need to set the flag
                        probVect.append((entry1[0], pos))
                else:
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
                    totLen1 = entry1[2]
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

            flag1, flag2 = getFlag(entry1, isLeaf=node1isleaf), getFlag(entry2, isLeaf=node2isleaf)
            if errorRateSiteSpecific: errorRate = errorRates[pos]
            if entry2[0] == entry1[0] and entry2[0] < 5:  # entry1 and entry2 are two identical nucleotides #Case τ1 = τ2 ϵ {A, C, T, G, R}
                end = min(entry1[1], entry2[1])
                probVect.append((entry2[0], end)) #no flag is needed, it can be deduced that this is a 0  branch length entry and flag=False.
                if returnLK:
                    if entry2[0] == 4: # case entry2 is R
                        cumulPartLk += (totLen1 + totLen2) * (cumulativeRate[end] - cumulativeRate[pos])
                    else:
                        if useRateVariation:
                            cumulPartLk += mutMatrices[pos][entry1[0]][entry1[0]] * (totLen1 + totLen2)
                        else:
                            cumulPartLk += nonMutRates[entry1[0]] * (totLen1 + totLen2)
                    if flag1 or flag2:
                        if errorRateSiteSpecific: cumErrorRate = cumulativeErrorRate[end] - cumulativeErrorRate[pos]
                        else: cumErrorRate = -errorRate * (end - pos)  # Approximation 2
                        cumulPartLk += cumErrorRate * (flag1 + flag2)
                pos = end
            elif (not totLen1) and (not totLen2) and entry1[0] < 5 and entry2[0] < 5 and (not flag2) and (not flag1):  # 0 distance between different nucleotides: merge is not possible
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

                newVec =  getPartialVec(i1, flag1, totLen1, errorRate, upNode=False)

                if entry2[0] == 6:  # entry1 is a nucleotide and entry2 is "O"
                    if totLen2:
                        for j in range4:
                            tot = 0.0
                            for i in range4:
                                tot += mutMatrix[j][i] * entry2[-1][i]
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
                        probVect.append((state, pos)) #no flag, because 0-branch length element and flag=False.
                    if returnLK:
                        cumulPartLk += log(sumV)
                else:  # entry1 and entry2 are nucleotides
                    if entry2[0] == 4:
                        i2 = refIndeces[pos]
                    else:
                        i2 = entry2[0]

                    if totLen2 or (flag2 and errorRate): # In contrast to mergeVectorsUpDown, here, in case totLen2==0, then the newVec will not be multiplied with the partial likelihood vector of entry2
                        partialVec2 = getPartialVec(i2, flag2, totLen2, errorRate,upNode=False)
                        for i in range4:
                            newVec[i] *= partialVec2[i]
                        sumV = sum(newVec)
                        for i in range4: #normalize
                            newVec[i] = newVec[i] / sumV
                        state = simplfy(newVec, refIndeces[pos])
                        pos += 1
                        if state == 6:
                            probVect.append((6, pos, newVec))
                        else:
                            probVect.append((state, pos)) #no flag, because 0-branch length element and flag=False
                        if returnLK:
                            cumulPartLk += log(sumV)
                    else:
                        pos += 1
                        probVect.append((entry2[0], pos))
                        if returnLK:
                            cumulPartLk += log(newVec[i2]) #should equal errorRate/3 if flag1 = true, flag2=False, totLen1=0, totLen2=0

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
                        probVect.append((state, pos)) #flag is False
                    if returnLK:
                        cumulPartLk += log(sumV)
                else:  # entry2 is a nucleotide and entry1 is "O"
                    if entry2[0] == 4:
                        i2 = refIndeces[pos]
                    else:
                        i2 = entry2[0]
                    if totLen2 or (flag2 and errorRate):  # 0-branch length but flag=T will also be covered.
                        partialVec2 = getPartialVec(i2, flag2, totLen2, errorRate, upNode=False)
                        for i in range4:
                            newVec[i] *= partialVec2[i]
                        sumV = sum(newVec)
                        for i in range4:
                            newVec[i] = newVec[i] / sumV
                        state = simplfy(newVec, refIndeces[pos])
                        pos += 1
                        if state == 6:
                            probVect.append((6, pos, newVec))
                        else:
                            probVect.append((state, pos)) #flag is false
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


def reCalculateWithErrors(root,mutMatrix, checkExistingAreCorrect=False,countNodes=False,
                          countPseudoCounts=False,pseudoMutCounts=None,
                          data=None,useRateVariation=False,mutMatrices=None, firstTimeError=False):
    """
    The errorrate version of 'reCalculateAllGenomeLists', Similarly to that, we loop in post-order, starting with calculating the lower-likelihoods of the leaves.
    After calculating all lowers, we will update other genome lists.
    This function is also used to add the error rates for the first time.
    :param firstTimeError: This parameter should be true when the error rates are included for the first time
    - only then should the partial likelihood vectors of the O-entries be updated for the leaf nodes.
    """
    node = root # direction 0 means from parent, direction 1 means from a child
    lastNode = None
    direction = 0
    while node != None:
        if direction == 0:
            if node.children:
                node = node.children[0]
            else: #leaf node:
                if firstTimeError:
                    addErrorTerminalNode(node)
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
                newLower = mergeVectors(node.children[0].probVect, node.children[0].dist, node.children[1].probVect,
                                        node.children[1].dist, mutMatrix, returnLK=False,
                                        useRateVariation=useRateVariation, mutMatrices=mutMatrices,
                                        node1isleaf= (node.children[0].children==[]),node2isleaf=(node.children[1].children==[]) )
                if checkExistingAreCorrect:
                    if areVectorsDifferentDebugging(newLower, node.probVect):
                        print("Inside reCalculateAllGenomeLists(), new lower at node is different from the old one, and it shouldn't be.")
                        print(newLower)
                        print(node.probVect)
                        raise Exception("exit")
                if newLower == None:
                    if not node.children[0].dist:
                        nodeList = []
                        updateBLen(nodeList, node, mutMatrix, useRateVariation=useRateVariation,mutMatrices=mutMatrices)
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
        newVect=rootVector(node.children[1].probVect,node.children[1].dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices, isLeaf=node.children[1].children==[])
        if checkExistingAreCorrect:
            if areVectorsDifferentDebugging(newVect,node.probVectUpRight):
                print("new probVectUpRight at root is different from the old one, and it shouldn't be.")
                print(newVect)
                print(node.probVectUpRight)
                raise Exception("exit")
        node.probVectUpRight=newVect
        newVect=rootVector(node.children[0].probVect,node.children[0].dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices, isLeaf=node.children[0].children==[])
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
                raise Exception("exit")
        node.probVectUpLeft=newVect

        #now traverse the tree downward and update the non-lower genome lists for all other nodes of the tree.
        lastNode=None
        node=node.children[0]
        direction=0
        while node!=None:
            if direction==0:
                if node==node.up.children[0]:
                    vectUp=node.up.probVectUpRight
                    nodeisleaf = not bool(node.up.children[1].children)
                else:
                    vectUp=node.up.probVectUpLeft
                    nodeisleaf = not bool(node.up.children[0].children)
                if node.dist:
                    #newVect=getTot(node,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                    if countPseudoCounts:
                        updatePesudoCounts(vectUp,node.probVect,pseudoMutCounts)

                    newVect=mergeVectorsUpDown(vectUp,node.dist/2,node.probVect,node.dist/2,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,
                                               node1isleaf= False,node2isleaf= node.children==[] ) #the upper likelihood cannot be from a leaf node, except if the branch length element is present
                    if checkExistingAreCorrect:
                        if areVectorsDifferentDebugging(newVect,node.probVectTotUp):
                            print("new probVectTotUp at node is different from the old one, and it shouldn't be.")
                            print(newVect)
                            print(node.probVectTotUp)
                            print(vectUp)
                            print(node.probVect)
                            newVect = mergeVectorsUpDown(vectUp, node.dist / 2, node.probVect, node.dist / 2, mutMatrix,
                                                         useRateVariation=useRateVariation, mutMatrices=mutMatrices,
                                                         node1isleaf=False,
                                                         node2isleaf=node.children == []) # the upper likelihood cannot be from a leaf node, except if the branch length element is present

                            print(newVect)
                    node.probVectTotUp=newVect
                    # if node.dist>=2*minBLenForMidNode:
                    # 	createFurtherMidNodes(node,vectUp,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices)
                else: # Changed this 20/9
                    if hasattr(node, 'probVectTotUp'): # Changed this 20/9
                        node.probVectTotUp = None # Changed this 20/9
                if node.children:
                    newUpRight=mergeVectorsUpDown(vectUp,node.dist,node.children[1].probVect,node.children[1].dist,mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,
                                                  node1isleaf= False, node2isleaf=  node.children[1].children ==[] )
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
                            raise Exception("exit")
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
                                raise Exception("exit")
                        node.probVectUpRight=newUpRight
                    newUpLeft=mergeVectorsUpDown(vectUp,node.dist,node.children[0].probVect,node.children[0].dist,
                                                 mutMatrix,useRateVariation=useRateVariation,mutMatrices=mutMatrices,
                                                 node1isleaf= False, #Vectup cannot be from a child
                                                 node2isleaf= node.children[0].children ==[])
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
                            raise Exception("exit")
                    else:
                        if checkExistingAreCorrect:
                            if areVectorsDifferentDebugging(newUpLeft,node.probVectUpLeft):
                                print("new probVectUpLeft at node is different from the old one, and it shouldn't be.")
                                print(node.children)
                                print(node.dist)
                                print(node.children[0].dist)
                                print(newUpLeft)
                                print()
                                print(node.probVectUpLeft) # original one does not account for errors.
                                print()
                                print(vectUp)
                                print()
                                print(node.children[0].probVect)
                                raise Exception("exit")
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


def errorRateEstimateBranchLengthWithDerivative(probVectP,probVectC,mutMatrix,useRateVariation=False,mutMatrices=None, node2isleaf=False):
    """
    The errorrate version of 'estimateBranchLengthWithDerivative',
    :param node2isleaf: True if probVectC is from a leaf node.
    """
    global errorRate
    c1=0.0
    ais=[]
    indexEntry1, indexEntry2, pos = 0, 0, 0
    entry1=probVectP[indexEntry1]
    entry2=probVectC[indexEntry2]
    end=min(entry1[1],entry2[1])
    while True:
        flag1, flag2 = getFlag(entry1, isLeaf=False), getFlag(entry2, isLeaf=node2isleaf)
        if errorRateSiteSpecific: errorRate = errorRates[pos]
        if entry2[0]==5: # case entry1 is NS
            pos=min(entry1[1],entry2[1]) #no error rate ?
            pass
        elif entry1[0]==5: # case entry2 is N
            #if parent node is type "N", in theory we might have to calculate the contribution of root nucleotides;
            # however, if this node is "N" then every other node in the current tree is "N", so we can ignore this since this contribution cancels out in relative terms.
            pos=min(entry1[1],entry2[1])
            pass #no error rate ?
        else:
            if entry1[0] < 5:
                if len(entry1) == 2:
                    contribLength = False
                elif len(entry1) == 4:  # ErrorRate; instead of ==3, ==4
                    contribLength = entry1[2]
                else:
                    contribLength = entry1[3]
            else:  # ErrorRate; stays the same for O and N entry
                if len(entry1) == 3:
                    contribLength = False
                else:
                    contribLength = entry1[2]
            if entry2[0] < 5:
                if len(entry2) == 4:  # ErrorRate; instead of ==3, ==4 (+ or ==3?)
                    if contribLength:
                        contribLength += entry2[2]
                    else:
                        contribLength = entry2[2]
            else:
                if len(entry2) == 4:  # ErrorRate; stays the same for O and N entries
                    if contribLength:
                        contribLength += entry2[2]
                    else:
                        contribLength = entry2[2]

            if entry1[0]==4: # case entry1 is R
                if entry2[0]==4:
                    end=min(entry1[1],entry2[1])
                    c1+=(cumulativeRate[end]-cumulativeRate[pos])
                    pos=end

                #entry1 is reference and entry2 is of type "O"
                elif entry2[0]==6:
                    if useRateVariation:
                        mutMatrix=mutMatrices[pos]
                    i1=refIndeces[pos]
                    if len(entry1)==5: # ErrorRate: ==4 to ==5
                        if flag1:
                            d = [0] * 4
                            denominator, numerator = 0, 0
                            for i in range4:
                                d[i] = rootFreqs[i] * ((i == i1) * (1 - 4 / 3 * errorRate) + mutMatrix[i][i1] * entry1[2] + errorRate / 3)
                                numerator += d[i] * entry2[-1][i]
                                denominator += mutMatrix[i1][i] * entry2[-1][i]
                            denominator *= rootFreqs[i1]
                            coeff0 = numerator
                            coeff1 = denominator
                            if contribLength:
                                coeff0 += contribLength * coeff1
                        else:
                            coeff0 = rootFreqs[i1] * entry2[-1][i1]
                            coeff1 = 0.0
                            for i in range4:
                                coeff0 += rootFreqs[i] * mutMatrix[i][i1] * entry1[2] * entry2[-1][i]
                                coeff1 += mutMatrix[i1][i] * entry2[-1][i]
                            coeff1 *= rootFreqs[i1]
                            if contribLength:
                                coeff0 += coeff1 * contribLength
                    else: #ErrorRate: nothing changes.
                        coeff0=entry2[-1][i1]
                        coeff1=0.0
                        for j in range4:
                            coeff1+=mutMatrix[i1][j]*entry2[-1][j]
                        if contribLength:
                            coeff0+=coeff1*contribLength
                    if coeff1<0.0:
                        c1+=coeff1/coeff0
                    elif coeff1:
                        coeff0=coeff0/coeff1
                        ais.append(coeff0)
                    pos+=1

                else: #entry1 is R and entry2 is a different but single nucleotide
                    i1 = refIndeces[pos]
                    i2 = entry2[0]
                    if len(entry1)==5: # ErrorRate: ==4 to ==5
                        if useRateVariation: #again f1 can only be true if the error rate is inherited across the root.
                            mutMatrix=mutMatrices[pos]
                        if flag1 or flag2:
                            c1up = entry1[2]
                            pi2pi1 = rootFreqs[i2]/rootFreqs[i1]
                            coeff0 = contribLength + (errorRate / 3 * flag2 + pi2pi1 * (
                                        mutMatrix[i2][i1] * c1up + errorRate / 3 * flag1)) / mutMatrix[i1][i2]
                            coeff1 = 1
                        else:
                            if contribLength:
                                coeff0=rootFreqs[i1]*mutMatrix[i1][i2]*contribLength\
                                       +rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]
                                coeff1=rootFreqs[i1]*mutMatrix[i1][i2]
                            else:
                                coeff0=rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]
                                coeff1=rootFreqs[i1]*mutMatrix[i1][i2]

                        coeff0=coeff0/coeff1
                    else: #ErrorRate : only f2 can be true.
                        if contribLength or flag2:
                            coeff0=contribLength + errorRate*flag2/(3*mutMatrix[i1][i2])
                        else:
                            coeff0=0.0
                    ais.append(coeff0)
                    pos+=1

            # entry1 is of type "O"
            elif entry1[0]==6:
                if useRateVariation:
                    mutMatrix=mutMatrices[pos]
                if entry2[0]==6: # entry1 and entry2 are both of type "O": no Error Rate
                    coeff0=entry1[-1][0]*entry2[-1][0]+entry1[-1][1]*entry2[-1][1]\
                           +entry1[-1][2]*entry2[-1][2]+entry1[-1][3]*entry2[-1][3] #coef0: prob that it stays the same.
                    coeff1=0.0
                    for i in range4:
                        for j in range4:
                            coeff1+=entry1[-1][i]*entry2[-1][j]*mutMatrix[i][j] #coef1: prob that it subsitutes
                    if contribLength: # I think this should be calculated only if coeff1>0
                        coeff0+=coeff1*contribLength
                    if coeff1 < 0.0:
                        c1 += coeff1 / coeff0
                    elif coeff1:
                        coeff0 = coeff0 / coeff1
                        ais.append(coeff0)
                else: #entry1 is "O" and entry2 is a nucleotide
                    if entry2[0]==4:
                        i2=refIndeces[pos]
                    else:
                        i2=entry2[0]
                    coeff1 = 0.0
                    for i in range4:
                        coeff1 += entry1[-1][i] * mutMatrix[i][i2]
                    if flag2: #include error rate if flag2
                        if coeff1 <0:# and errorRate/entry1[-1][i2] <0.01:
                            #second if statment that only does this approximation if the other terms are small
                            coeff0 = entry1[-1][i2]
                            if contribLength:
                                coeff0 += coeff1 * contribLength  #I temperorarly added this: otherwise not be the same as the normal function for cases where error rate =0
                            c1 += coeff1/coeff0
                        elif coeff1:
                            coeff0 = (entry1[-1][i2]+ 1/3*errorRate*(1-4*entry1[-1][i2]))/coeff1 + contribLength
                            ais.append(coeff0)
                    else:
                        coeff0 = entry1[-1][i2]
                        if contribLength:
                            coeff0+=coeff1*contribLength
                        if coeff1 < 0.0:
                            c1 += coeff1 / coeff0
                        elif coeff1:
                            coeff0 = coeff0 / coeff1
                            ais.append(coeff0)

                pos+=1

            else: #entry1 is a non-ref nuc
                if useRateVariation:
                    mutMatrix=mutMatrices[pos]
                if entry2[0]==entry1[0]:
                    c1+=mutMatrix[entry1[0]][entry1[0]] #ErrorRate: nothing changes
                else: #entry1 is a nucleotide and entry2 is not the same as entry1
                    i1=entry1[0]
                    if entry2[0]<5: #entry2 is a nucleotide
                        if entry2[0]==4:
                            i2=refIndeces[pos]
                        else:
                            i2=entry2[0]

                        if len(entry1)==5: # ErrorRate: ==4 to ==5
                            if flag1 or flag2: #although the calculation is the same, I added 'if errorRate' for testing purposes, because the last decimal is a bit different.
                                c1up = entry1[2]
                                pi2pi1 = rootFreqs[i2] / rootFreqs[i1]
                                coeff0 = contribLength + (errorRate/3*flag2 + pi2pi1*(mutMatrix[i2][i1]*c1up + errorRate/3*flag1))/mutMatrix[i1][i2]
                                coeff1 = 1
                            else:
                                if contribLength:
                                    coeff0=rootFreqs[i1]*mutMatrix[i1][i2]*contribLength+rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]
                                    coeff1=rootFreqs[i1]*mutMatrix[i1][i2]
                                else:
                                    coeff0=rootFreqs[i2]*mutMatrix[i2][i1]*entry1[2]
                                    coeff1=rootFreqs[i1]*mutMatrix[i1][i2]
                            coeff0=coeff0/coeff1
                        else:
                            if contribLength or flag2:
                                coeff0=contribLength + errorRate*flag2/(3*mutMatrix[i1][i2])
                            else:
                                coeff0=0.0
                        ais.append(coeff0)

                    else: #entry1 is a nucleotide and entry2 is of type "O"
                        if len(entry1)==5: # ErrorRate: ==4 to ==5
                            if flag1: # ErrorRate:
                                d = [0] * 4
                                denominator, numerator = 0, 0
                                for i in range4:
                                    d[i] = rootFreqs[i] * ((i == i1) * (1 - 4 / 3 * errorRate) + mutMatrix[i][i1] * entry1[2] + errorRate / 3)
                                    numerator += d[i] * entry2[-1][i]
                                    denominator += mutMatrix[i1][i] * entry2[-1][i]
                                denominator *= rootFreqs[i1]
                                coeff0 = numerator
                                coeff1 = denominator
                                if contribLength:
                                    coeff0 += contribLength * coeff1
                            else:
                                coeff0=rootFreqs[i1]*entry2[-1][i1]
                                coeff1=0.0
                                for i in range4:
                                    coeff0+=rootFreqs[i]*mutMatrix[i][i1]*entry1[2]*entry2[-1][i]
                                    coeff1+=mutMatrix[i1][i]*entry2[-1][i]
                                coeff1*=rootFreqs[i1]
                                if contribLength:
                                    coeff0+=coeff1*contribLength
                        else: # errorrate: nothing changes
                            coeff0=entry2[-1][i1]
                            coeff1=0.0
                            for j in range4:
                                coeff1+=mutMatrix[i1][j]*entry2[-1][j]
                            if contribLength:
                                coeff0+=coeff1*contribLength
                        if coeff1<0.0:
                            c1+=coeff1/coeff0
                        elif coeff1:
                            coeff0=coeff0/coeff1
                            ais.append(coeff0)
                pos+=1
        #print(c1, ais, pos)
        if pos==lRef:
            break
        if pos==entry1[1]:
            indexEntry1+=1
            entry1=probVectP[indexEntry1]
        if pos==entry2[1]:
            indexEntry2+=1
            entry2=probVectC[indexEntry2]

    #now optimized branch length based on coefficients
    #minBLenSensitivity=minBLen/100000
    c1=-c1
    n=len(ais)
    if n==0:
        return False
    else:
        tDown = n / c1 - min(ais)
        if tDown <= 0.0:
            return 0.0
        # vDown=calculateDerivative(ais,tDown)
        vDown = 0.0
        for ai in ais:
            vDown += 1.0 / (ai + tDown)
        tUp = n / c1 - max(ais)
        # TODO Myrthe, changed here
        if tUp <= minBLenSensitivity:
            if min(ais):
                tUp = 0.0
            else:
                tUp = minBLenSensitivity
        vUp = 0.0
        for ai in ais:
            vUp += 1.0 / (ai + tUp)
    if vDown>c1+minBLenSensitivity or vUp<c1-minBLenSensitivity:
        if vUp<c1-minBLenSensitivity and (not tUp):
            return 0.0
        print(c1,ais)
        print("Initial values")
        print(tDown,vDown,tUp,vUp)
        print("Initial border parameters don't fit expectations:")
        print(c1,ais)
        print(probVectP)
        print(probVectC)
        print(mutMatrix)
        #raise Exception("exit")

    while tDown-tUp>minBLenSensitivity:
        tMiddle=(tUp+tDown)/2
        vMiddle=calculateDerivative(ais,tMiddle)
        if vMiddle>c1:
            tUp=tMiddle
        else:
            tDown=tMiddle
    return tUp


#for the root, take lower likelihood genome list probVect, and create an overall likelihood (or upper right or upper left) genome list by multiplying likelihoods by root frequencies.
def rootVectorErrorRate(probVect,bLen,mutMatrix,useRateVariation=False,mutMatrices=None, isLeaf=False):
    """
    The errorrate version of 'rootVector',
    :param isLeaf: True if ProbVect is from a leaf node.
    """
    newProbVect=[]
    for entry in probVect:
        if entry[0]==5: #==N
            newProbVect.append(entry)
        elif entry[0]==6:
            if len(entry)==4:
                totBLen=entry[2]
                if bLen:
                    totBLen+=bLen
            else:
                totBLen=bLen
            newVect=[]
            if totBLen:
                if useRateVariation:
                    mutMatrix=mutMatrices[entry[1]-1]
                for i in range4:
                    tot=0.0
                    for j in range4:
                        tot+=mutMatrix[i][j]*entry[-1][j]
                    tot*=totBLen
                    tot+=entry[-1][i]
                    newVect.append(tot*rootFreqs[i])
                totSum=sum(newVect)
                for i in range4:
                    newVect[i]/=totSum
                newProbVect.append((6,entry[1],newVect))
            else:
                for i in range4:
                    newVect.append(entry[-1][i]*rootFreqs[i])
                totSum=sum(newVect)
                for i in range4:
                    newVect[i]/=totSum
                newProbVect.append((6,entry[1],newVect))
        else: #other nucleotide
            if len(entry)==4: #ErrorRate ==4 instead of ==3
                if bLen:
                    newProbVect.append((entry[0],entry[1],entry[2]+bLen,0.0, entry[3])) #ErrorRate: add flag
                else:
                    newProbVect.append((entry[0],entry[1],entry[2],0.0,entry[3])) #ErrorRate: add flag
            else:
                if bLen or isLeaf:
                    newProbVect.append((entry[0],entry[1],bLen,0.0, isLeaf)) #ErrorRate: add flag
                else:
                    newProbVect.append((entry[0],entry[1]))
    return newProbVect

def areVectorsDifferentErrorRate(probVect1,probVect2):
    """
    The errorrate function of 'areVectorsDifferent',
    """
    #should both have flags
    if probVect2==None:
        return True
    indexEntry1, indexEntry2, pos = 0, 0, 0
    entry1=probVect1[indexEntry1]
    entry2=probVect2[indexEntry2]
    while True:
        if entry1[0]!=entry2[0]:
            return True
        if len(entry1)!=len(entry2):
            return True
        if entry1[0]<5:
            if len(entry1)>2:
                assert(len(entry1)!=3)
                if abs(entry1[2] - entry2[2])>thresholdProb: #first bLen element
                     return True
                if len(entry1)==5: #==5 instead of ==4
                    if abs(entry1[3] - entry2[3])>thresholdProb: #2nd blen element
                        return True
                if entry1[-1] != entry2[-1]: #error rate flags
                    return True
        if entry1[0]==6:
            if len(entry1)==4:
                if abs(entry1[2] - entry2[2])>thresholdProb:
                    return True
            for i in range4:
                diffVal=abs(entry1[-1][i] - entry2[-1][i])

                if diffVal:
                    if (not entry1[-1][i]) or (not entry2[-1][i]):
                        return True
                    if diffVal>thresholdDiffForUpdate or (diffVal>thresholdProb and ( (diffVal/entry1[-1][i]>thresholdFoldChangeUpdate)  or  (diffVal/entry2[-1][i]>thresholdFoldChangeUpdate))):
                        return True
        #update pos, end, etc
        pos=min(entry1[1],entry2[1])
        if pos==lRef:
            break
        if pos==entry1[1]:
            indexEntry1+=1
            entry1=probVect1[indexEntry1]
        if pos==entry2[1]:
            indexEntry2+=1
            entry2=probVect2[indexEntry2]
    return False

def traverseTopology(node, function, leafsOnly=False):
    """
    Traverse a topology and perform a function
    :param node: the root node.
    :param function: the function to be performed at every node. It should take only 'node' as argument
    :param leafsOnly: If True, the function will only be performed at the leaf node
    """
    nodesToVisit=[node]
    while nodesToVisit:
        newNode=nodesToVisit.pop()
        for c in newNode.children:
            nodesToVisit.append(c)
        if not newNode.children or not leafsOnly:
            function(newNode)

def countEntries(node):
    numNodes[0]+=1
    for entry in node.probVect: #or should i also check the node.probvectup etc?
        if entry[0] < 4:
            numNodes[1] += 1
        elif entry[0] == 4:
            numNodes[2] += 1
        elif entry[0] == 5:
            numNodes[3] += 1
        else:
            numNodes[4] += 1

def countEntriesAll(root):
    global numNodes
    numNodes = [0]*5
    traverseTopology(root, countEntries)
    print("Number of nodes: "+str(numNodes[0]))
    print("Os per node: "+str(float(numNodes[4])/numNodes[0]))
    print("Nucs per node: "+str(float(numNodes[1])/numNodes[0]))
    print("Ns per node: "+str(float(numNodes[3])/numNodes[0]))

def counttotBLenAll(tree):
    """
    Traverse the tree to count the sum of all the branch lengths.
    """
    root = getRoot(tree)
    global totBLen
    totBLen = 0
    def counttotBLen(node):
        global totBLen
        totBLen += node.dist
    traverseTopology(root, counttotBLen)
    print("Tot branch length: "+str(totBLen))
    return totBLen - root.dist

def countFlagsAll(root):
    """
    Traverse the tree to count the number of actually stored flags, and how many are set to true.
    """
    global totFlags
    totFlags = 0
    global totFlagsTrue
    totFlagsTrue = 0
    global totMinorSeqs
    totMinorSeqs = 0
    global totFlaggedElems
    totFlaggedElems =0
    def countFlags(node):
        global totFlags
        global totFlaggedElems
        global totMinorSeqs
        global totFlagsTrue
        pos=0
        if hasattr(node, 'minorSequences'):
            totMinorSeqs += len(node.minorSequences)
        for entry in node.probVect:  # or should i also check the node.probvectup etc?
            totFlags += 1 if (type(entry[-1]) == bool) else 0
            totFlagsTrue += 1 if (type(entry[-1]) == bool and entry[-1] == True) else 0
            totFlaggedElems += getFlag(entry, node.children==[])*(entry[1]-pos)
            pos = entry[1]
    traverseTopology(root, countFlags)
    print("Tot #flags stored: "+str(totFlags)) # was only 79 in a set of 174seq*1500 nucleotides.
    print("\tTot #flags set to True: "+str(totFlagsTrue)) # was only 79 in a set of 174seq*1500 nucleotides.
    print("\tTot #flags set to False: "+str(totFlags - totFlagsTrue)) # was only 79 in a set of 174seq*1500 nucleotides.
    print("Tot #positions where a flag(=True) will be taken into account: "+str(totFlaggedElems))
    print("#minor sequences: "+str(totMinorSeqs))
    return totFlags

def traverseTwoTopologies(node1, node2):
    """
    To aid testing if error rate functions are correct.
    Traverse two topologies, one with flags, one without flags,
    and test if the results of the error rate functions are the same as the normal MAPLE functions
    :param node1: t1, the tree with the error rates and flags
    :param node2: t0_1, the tree without flags
    :return:
    """
    decimal_places=4
    nodesToVisit=[node1]
    traversalList1 = [node1]
    while nodesToVisit:
        newNode=nodesToVisit.pop()
        for c in newNode.children:
            nodesToVisit.append(c)
            traversalList1.append(c)
    nodesToVisit = [node2]
    traversalList2 = [node2]
    while nodesToVisit:
        newNode=nodesToVisit.pop()
        for c in newNode.children:
            nodesToVisit.append(c)
            traversalList2.append(c)
    for node1, node2 in zip(traversalList1, traversalList2):
        if not (round(findProbRootError(node1.probVect),decimal_places) == round(findProbRootOriginal(node2.probVect),decimal_places)):
            print('findProbRootError incorrect')
            print(findProbRootError(node1.probVect))
            print(findProbRootOriginal(node2.probVect))
        if areVectorsDifferentDebugging(node1.probVect,node2.probVect): # mergeVectors, rootVectors, MergeVectorsupdown and reCalculateAllGenomeLists seem to work as original for errorRate=0
            print("vectors are different")
        if hasattr(node1, 'dist'):
            if round(node1.dist,decimal_places) != round(node2.dist,decimal_places):
                    print("blen is wrong? ")
                    print(node1.dist)
                    print(node2.dist)
        if node1.children:
            if hasattr(node1, 'probVectUpLeft') and hasattr(node1, 'probVectUpRight'):
                if areVectorsDifferentDebugging(node1.probVectUpRight, node2.probVectUpRight):  # mergeVectors, rootVectors, MergeVectorsupdown and reCalculateAllGenomeLists seem to work as original for errorRate=0
                    print("vectors are different")
                if areVectorsDifferentDebugging(node1.probVectUpLeft, node2.probVectUpLeft):  # mergeVectors, rootVectors, MergeVectorsupdown and reCalculateAllGenomeLists seem to work as original for errorRate=0
                    print("vectors are different")
                if not (round(errorRateEstimateBranchLengthWithDerivative(node1.probVectUpRight,node1.children[0].probVect, mutMatrix),decimal_places) ==
                        round(estimateBranchLengthWithDerivativeOriginal(node2.probVectUpRight,node2.children[0].probVect, mutMatrix),decimal_places)):
                    print('errorRateEstimateBranchLengthWithDerivative with probVectUpRicht incorrect')
                    print(errorRateEstimateBranchLengthWithDerivative(node1.probVectUpRight,node1.children[0].probVect, mutMatrix))
                    print(estimateBranchLengthWithDerivativeOriginal(node2.probVectUpRight,node2.children[0].probVect, mutMatrix))
                if not (round(errorRateEstimateBranchLengthWithDerivative(node1.probVectUpLeft,node1.children[1].probVect, mutMatrix),decimal_places) ==
                        round(estimateBranchLengthWithDerivativeOriginal(node2.probVectUpLeft,node2.children[1].probVect, mutMatrix),decimal_places)):
                    print('errorRateEstimateBranchLengthWithDerivative with probVectUpLeft incorrect')
                    print(errorRateEstimateBranchLengthWithDerivative(node1.probVectUpLeft,node1.children[1].probVect, mutMatrix))
                    print(estimateBranchLengthWithDerivativeOriginal(node2.probVectUpLeft,node2.children[1].probVect, mutMatrix))
            if not (round(appendProbNodeErrorRate(node1.children[0].probVect, node1.children[1].probVect, node1.children[1].dist, mutMatrix),decimal_places)) ==\
                   (round(appendProbNodeOriginal(node2.children[0].probVect, node2.children[1].probVect, node2.children[1].dist, mutMatrix),decimal_places)):#, useRateVariation=True, mutMatrices=mutMatrices))
                print('appendProbNodeErrorRate incorrect')
                print(appendProbNodeErrorRate(node1.children[0].probVect, node1.children[1].probVect, node1.children[1].dist, mutMatrix))
                print(appendProbNodeOriginal(node2.children[0].probVect, node2.children[1].probVect, node2.children[1].dist, mutMatrix))
            if node1.up and node2.up:
                upVect1 = node1.up.probVectUpRight if node1==node1.up.children[0] else node1.up.probVectUpLeft
                upVect2 = node2.up.probVectUpRight if node2==node2.up.children[0] else node2.up.probVectUpLeft
                if not (round(appendProbNodeErrorRate(upVect1, node1.children[1].probVect, node1.children[1].dist, mutMatrix),decimal_places)) == \
                       (round(appendProbNodeOriginal(upVect2, node2.children[1].probVect, node2.children[1].dist,mutMatrix),decimal_places)):  # , useRateVariation=True, mutMatrices=mutMatrices))
                        print('appendProbNodeErrorRate incorrect')
                        print(round(appendProbNodeErrorRate(upVect1, node1.children[1].probVect, node1.children[1].dist, mutMatrix),decimal_places))
                        print(round(appendProbNodeOriginal(upVect2, node2.children[1].probVect, node2.children[1].dist,mutMatrix),decimal_places))  # , useRateVariation=True, mutMatrices=mutMatrices)



totalEntries = 0
def countNumEntries(node):
    global totalEntries
    totalEntries += len(node.probVect)

def overallLeaf(node):
    """
    This function should aid error detection and correction.
    It calculates the overall likelihood of each leaf node its upper likelihood vector and the lower vector, in which the error rates are taken into account.
    The function mergeVectorsUpDownErrorDetection The vectup will be transfered down over the branch length leading up the node. This is     #Similar to what is normally calculated for ProbVectTotUp of leaf nodes.

    After calculating the Overall likelihood, the function loops over the overall and lower likelihood simulteanously, to see what was observed at the position.
    It subsequently  prints out the posterior error probabilities for each position or fragment.
    """
    # if a flagged entry is at an inner node, one can calculate error probabilities, but most then transfer it to the leaf node from which this arrived.
    global mutMatrices, useRateVariation
    if node.up:
        print(node)
        if node == node.up.children[0]:
            vectUp = node.up.probVectUpRight
        else:
            vectUp = node.up.probVectUpLeft
        node.overall = mergeVectorsUpDownErrorDetection(vectUp, node.dist, node.probVect, 0, mutMatrix, node2isleaf=True, mutMatrices=None, useRateVariation=False)
        print(node.overall)
        print(" ".join([["A", "C", "G", "T", "R", "N", "O"][x[0]] + str(x[1]) for x in node.probVect]))
        print(" ".join([["A", "C", "G", "T", "R", "N", "O"][x[0]] + str(x[1]) for x in vectUp]))
        print(" ".join([["A", "C", "G", "T", "R", "N", "O"][x[0]] + str(x[1]) for x in node.overall]))
        if not node.children: #loop over overall and lower likelihood, to see what was observed at the position.
            indexEntry1, indexEntry2, pos = 0, 0, 0
            errorProbs = []
            entry1 = node.overall[indexEntry1]
            entry2 = node.probVect[indexEntry2]
            while True:
                if entry2[0] < 5 and (entry1[0] == 6 or (entry1[0] <5 and entry1[0]!=entry2[0])):
                    if entry2[0] == 4:
                        i2 = refIndeces[pos]
                    else:
                        i2 = entry2[0]
                    if entry1[0]==6:
                        errorProbs.append((pos,1 - entry1[-1][i2]))
                    else:
                        errorProbs.append((pos, 1))
                else:
                    errorProbs.append((pos, 0))
                pos = min(entry1[1], entry2[1])
                if pos == lRef:
                    break
                if pos == entry1[1]:
                    indexEntry1 += 1
                    entry1 = node.overall[indexEntry1]
                if pos == entry2[1]:
                    indexEntry2 += 1
                    entry2 = node.probVect[indexEntry2]
            node.errorProbs = errorProbs
            print(errorProbs)


def activateErrorFunctions(activate=True, activatefirsttime=False):
    "#replace normal functions with pointers to their error rate variants. "
    global areVectorsDifferentOriginal, estimateBranchLengthWithDerivativeOriginal, mergeVectorsUpDownOriginal, mergeVectorsOriginal, rootVectorOriginal, appendProbNodeOriginal, findProbRootOriginal, reCalculateAllGenomeListsOriginal
    global estimateBranchLengthWithDerivative, areVectorsDifferent, mergeVectorsUpDown, mergeVectors, rootVector, appendProbNode, findProbRoot, reCalculateAllGenomeLists
    import copy #todo check if this works without the copy library, by simply stating estimateBranchLengthWithDerivativeOriginal =estimateBranchLengthWithDerivative
    if activate:
        if activatefirsttime:
            estimateBranchLengthWithDerivativeOriginal = copy.deepcopy(estimateBranchLengthWithDerivative)
            mergeVectorsOriginal = copy.deepcopy(mergeVectors)
            mergeVectorsUpDownOriginal = copy.deepcopy(mergeVectorsUpDown)
            rootVectorOriginal = copy.deepcopy(rootVector)
            appendProbNodeOriginal = copy.deepcopy(appendProbNode)
            findProbRootOriginal = copy.deepcopy(findProbRoot)
            reCalculateAllGenomeListsOriginal = copy.deepcopy(reCalculateAllGenomeLists)
            areVectorsDifferentOriginal = copy.deepcopy(areVectorsDifferent)
            #
        areVectorsDifferent = areVectorsDifferentErrorRate
        estimateBranchLengthWithDerivative = errorRateEstimateBranchLengthWithDerivative
        mergeVectors = mergeVectorsError
        mergeVectorsUpDown = mergeVectorsUpDownError
        rootVector = rootVectorErrorRate
        appendProbNode = appendProbNodeErrorRate
        findProbRoot = findProbRootError
        reCalculateAllGenomeLists =reCalculateWithErrors
    else:
        estimateBranchLengthWithDerivative = estimateBranchLengthWithDerivativeOriginal
        mergeVectors = mergeVectorsOriginal
        mergeVectorsUpDown = mergeVectorsUpDownOriginal
        rootVector = rootVectorOriginal
        appendProbNode = appendProbNodeOriginal
        findProbRoot = findProbRootOriginal
        reCalculateAllGenomeLists =reCalculateAllGenomeListsOriginal
        areVectorsDifferent = areVectorsDifferentOriginal


if errorRate or errorRateSiteSpecific:
    activateErrorFunctions(activate=True, activatefirsttime=True) #replace normal functions with pointers to their error rate variants.
    if debugging: print("start recalculating with ErrorRate")
    reCalculateWithErrors(t1, mutMatrix, useRateVariation=rateVariation,mutMatrices=mutMatrices,firstTimeError=True)#, errorRate= errorRate)
    if debugging: print("ErrorRate included")
    if debugging: reCalculateWithErrors(t1, mutMatrix, checkExistingAreCorrect=True, useRateVariation=rateVariation,mutMatrices=mutMatrices,firstTimeError=False)

#--------------------------------------------

#if asked, first run a round of short range topology search
timeTopology=0.0
if fastTopologyInitialSearch and (inputTree=="" or largeUpdate):
    start=time()
    #if inputTree=="" or largeUpdate:
    setAllDirty(t1)
    newRoot,improvement=startTopologyUpdates(t1,mutMatrix,checkEachSPR=debugging,strictTopologyStopRules=strictTopologyStopRulesInitial,allowedFailsTopology=allowedFailsTopologyInitial,thresholdLogLKtopology=thresholdLogLKtopologyInitial,thresholdTopologyPlacement=thresholdTopologyPlacementInitial,useRateVariation=rateVariation,mutMatrices=mutMatrices)
    if newRoot!=None:
        t1=newRoot
    timeForUpdatingTopology=(time()-start)
    print("Time for initial traversal of the tree for only short range topology changes or branch length: "+str(timeForUpdatingTopology))
    timeTopology+=timeForUpdatingTopology
    print("LK improvement apparently brought: "+str(improvement))

    #run improvements only on the nodes that have been affected by some changes in the last round, and so on
    start=time()
    subRound=0
    while subRound<20:
        print("Topological subround "+str(subRound+1))
        newRoot,improvement=startTopologyUpdates(t1,mutMatrix,checkEachSPR=debugging,strictTopologyStopRules=strictTopologyStopRulesInitial,allowedFailsTopology=allowedFailsTopologyInitial,thresholdLogLKtopology=thresholdLogLKtopologyInitial,thresholdTopologyPlacement=thresholdTopologyPlacementInitial,useRateVariation=rateVariation,mutMatrices=mutMatrices)
        if newRoot!=None:
            t1=newRoot
        print("LK improvement apparently brought: "+str(improvement))
        if improvement<thresholdLogLKwholeTopologyImprovement:
            break
        subRound+=1
    timeForUpdatingTopology=(time()-start)
    print("Time for the subrounds of this traversal of the tree: "+str(timeForUpdatingTopology))
    timeTopology+=timeForUpdatingTopology

#now run topological improvements
debugging=False
for i in range(numTopologyImprovements):
    print("Starting topological impromevement attempt traversing number "+str(i+1))
    start=time()
    if inputTree=="" or largeUpdate:
        setAllDirty(t1)
    newRoot,improvement=startTopologyUpdates(t1,mutMatrix,checkEachSPR=debugging,useRateVariation=rateVariation,mutMatrices=mutMatrices)
    if newRoot!=None:
        t1=newRoot
    timeForUpdatingTopology=(time()-start)
    print("Time for this traversal of the tree to update the topology or branch length: "+str(timeForUpdatingTopology))
    timeTopology+=timeForUpdatingTopology
    print("LK improvement apparently brought: "+str(improvement))

    if runOnlyExample:
        print("Tree after improvement "+str(i))
        newickString=createNewick(t1)
        print(newickString)

    if improvement<thresholdLogLKwholeTopologyImprovement:
        print("Small improvement, stopping topological search.")
        break

    #run improvements only on the nodes that have been affected by some changes in the last round, and so on
    start=time()
    subRound=0
    while subRound<20:
        print("Topological subround "+str(subRound+1))
        newRoot,improvement=startTopologyUpdates(t1,mutMatrix,checkEachSPR=debugging,useRateVariation=rateVariation,mutMatrices=mutMatrices)
        if newRoot!=None:
            t1=newRoot
        print("LK improvement apparently brought: "+str(improvement))
        if improvement<thresholdLogLKwholeTopologyImprovement:
            break
        subRound+=1
    timeForUpdatingTopology=(time()-start)
    print("Time for the subrounds of this traversal of the tree: "+str(timeForUpdatingTopology))
    timeTopology+=timeForUpdatingTopology
    if not (inputTree=="" or largeUpdate):
        break

testingBLen=False
if testingBLen:
    #calculate total likelihood
    if calculateLKfinalTree:
        totalLK=calculateTreeLikelihood(t1,mutMatrix,useRateVariation=rateVariation,mutMatrices=mutMatrices)
        print("totalLK before branch length optimization: "+str(totalLK))

#optimize branch lengths of the final tree
if optimizeBranchLengths:
    start=time()
    setAllDirty(t1)
    print("Branch length optimization")
    improvement=traverseTreeToOptimizeBranchLengths(t1,mutMatrix,testing=testingBLen,useRateVariation=rateVariation,mutMatrices=mutMatrices)
    if testingBLen:
        print("Branch length optimization round 1 approximate expected improvement: "+str(improvement))
    subRound=0
    while subRound<20:
        if improvement<thresholdLogLKwholeTopologyImprovement:
            break
        subRound+=1
        improvement=traverseTreeToOptimizeBranchLengths(t1,mutMatrix,useRateVariation=rateVariation,mutMatrices=mutMatrices)
        if testingBLen:
            print("branch length finalization subround "+str(subRound+1)+" improvement "+str(improvement))
        else:
            print("branch length finalization subround "+str(subRound+1))
    timeForBranchOptimization=(time()-start)
    print("Time for updating branch lengths: "+str(timeForBranchOptimization))


#calculate total likelihood
if calculateLKfinalTree:
    totalLK=calculateTreeLikelihood(t1,mutMatrix,useRateVariation=rateVariation,mutMatrices=mutMatrices)
    print("totalLK: "+str(totalLK))

#Free space by deleting the genome lists in the tree.
nextLeaves=[t1]
while nextLeaves:
    node=nextLeaves.pop()
    node.probVect=None
    if node.dist:
        node.probVectTot=None
        if node.up!=None:
            node.probVectTotUp=None
            node.furtherMidNodes=None
    if node.children:
        node.probVectUpRight=None
        node.probVectUpLeft=None
        for c in node.children:
            nextLeaves.append(c)
print("Deleted genome lists.")

if not debugging and inputTree=="":
    #Read input file to collect all sample names
    fileI=open(inputFile)
    line=fileI.readline()
    if refFile=="":
        line=fileI.readline()
        while line[0]!=">":
            line=fileI.readline()
    sampleNames=[]
    while line!="" and line!="\n":
        name=line.replace(">","").replace("\n","")
        line=fileI.readline()
        while line!="" and line!="\n" and line[0]!=">":
            line=fileI.readline()
        sampleNames.append(name)
    fileI.close()
    print("Sample names read.")

    #put sample names in the tree
    nextLeaves=[t1]
    name=""
    while nextLeaves:
        node=nextLeaves.pop()
        if not node.children:
            if debugging:
                name=sampleNames[int(node.name[1:])]  #strip the 'S' off here, which is appended sopmewhere earlier when debugging is True.
                sampleNames[int(node.name[1:])]=None
            else:
                name=sampleNames[node.name]
                sampleNames[node.name]=None
            node.name=name
            for m in range(len(node.minorSequences)):
                if debugging:
                    old_name = int(node.minorSequences[m][1:])
                else: old_name = int(node.minorSequences[m])
                name=sampleNames[old_name]
                sampleNames[old_name]=None
                node.minorSequences[m]=name
        else:
            for c in node.children:
                nextLeaves.append(c)
    print("Sample names assigned.")


if binaryTree:
    newickString=createBinaryNewick(t1)
else:
    newickString=createNewick(t1)

if runOnlyExample:
    print("Final tree:")
    print(newickString)
    raise Exception("exit")

file=open(outputFile+"_tree.tree","w")
file.write(newickString)
file.close()
file=open(outputFile+"_subs.txt","w")
for i in range4:
    for j in range4:
        file.write(str(mutMatrix[i][j])+"\t")
    file.write("\n")
if rateVariation:
    file.write("\n\n"+"Site rates:")
    for i in range(lRef):
        file.write(str(i+1)+"\n"+str(siteRates[i])+"\n")
file.close()
print("Missed minor samples: "+str(totalMissedMinors[0]))
print("Final Substitution matrix:")
print(mutMatrix)
print("Time spent finding placement nodes: "+str(timeFinding))
print("Time spent placing samples on the tree: "+str(timePlacing))
print("Time spent in total updating the topology and branch lengths: "+str(timeTopology))
print("successfully (?) finished code")
runTime = time() - startTime
############################
if benchmarkingFile:
    # Required arguments: - benchmarkingFile to append to!
    trueTree = readNewick(trueTree, divideBranchLengthsBy=lRef)[0]
    leafNameDict, nodeTable, leafCount, numBranches, leafDistDict, branchLengthDict, sumBranchLengths = prepareTreeComparison(trueTree, rooted=True, addRootRFL = False)
    estimatedTree = readNewick(outputFile+"_tree.tree")[0] #simply using t1 did not seem to work (the RF would for instance not be 0, when comparing with the identical tree for instance), so I load the newly estimated tree from the newick format.
    numDiffs, normalisedRF, leafCount, foundBranches, missedBranches, notFoundBranches, RFL =\
        RobinsonFouldsWithDay1985(estimatedTree, leafNameDict, nodeTable,  leafCount,numBranches, leafDistDict, branchLengthDict, sumBranchLengths, rooted=True, addRootRFL = False)
    totBlenTrue = counttotBLenAll(trueTree)
    totBlenEstimated = counttotBLenAll(estimatedTree)

    if not os.path.exists(benchmarkingFile): # write a header if the file does not yet exist.
        file = open(benchmarkingFile, "w")
        file.write( "timeOfJob\t"  +"inputFile\t" + "repeat\t" + "errorRateInInference\t" + "errorRateInSimulation\t" + "siteSpecificInference\t" + "siteSpecificSimulation\t"  +"lRef\t"+ "leaves\t" + "||\t" +
                    "runtime\t" + "LK\t" + "RF\t" + "normalisedRF\t" + "foundBranches\t" + "missedBranches\t" + "notFoundBranches\t" + "RFL\t" + "totalBranchLength\t" + "totalBranchLengthTrue\n")
        file.close()
    repeat = "None"
    errorRateSimulated = "None"
    errorRateSiteSpecificSimulated = False
    try:
        splitFileName = inputFile[:-4].split("_")
        for item in splitFileName:
            if 'repeat' in item:
                repeat = item[6:]
            elif 'errors' in item:
                errorRateSimulated = item[6:]
            elif 'sitespecific' in item:
                errorRateSiteSpecificSimulated = True
    except:
        None
    file = open(benchmarkingFile, "a")
    file.write(str(time()) + "\t" + str(inputFile) + "\t" + repeat + "\t" + str(errorRate) + "\t"  + str(errorRateSimulated) + "\t"   + str(bool(errorRateSiteSpecific)) + "\t" +  str(bool(errorRateSiteSpecificSimulated)) +  "\t"  + str(lRef) + "\t" + str(leafCount) + "\t||\t" +
               str(runTime) + "\t" + str(totalLK) + "\t" + str(numDiffs) + "\t" + str(normalisedRF)  + "\t" + str(foundBranches) + "\t" + str(missedBranches) + "\t" + str(notFoundBranches) + "\t" + str(RFL) + "\t" + str(totBlenEstimated)+ "\t"  + str(totBlenTrue)+ "\n")
    print(str(time()) + "\t" + str(inputFile) + "\t" + repeat + "\t" + str(errorRate) + "\t"  + str(errorRateSimulated) + "\t"   + str(bool(errorRateSiteSpecific)) + "\t" +  str(bool(errorRateSiteSpecificSimulated)) +  "\t"  + str(lRef) + "\t" + str(leafCount) + "\t||\t" +
               str(runTime) + "\t" + str(totalLK) + "\t" + str(numDiffs) + "\t" + str(normalisedRF)  + "\t" + str(foundBranches) + "\t" + str(missedBranches) + "\t" + str(notFoundBranches) + "\t" + str(RFL) + "\t" + str(totBlenEstimated)+ "\t"  + str(totBlenTrue)+ "\n")
    file.close()

############################
exit()




















