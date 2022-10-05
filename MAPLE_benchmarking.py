import sys
import os.path
from time import time
import argparse
import random
import time

# todo: call my simulate errors from terminal.
# call my MAPLE from terminal.


#©EMBL-European Bioinformatics Institute, 2021

#Prepare files to run MAPLE simulations/benchmarking
#assumes tree for simulations has been downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/
#and that a modified version of phastSim is used that allows to read such tree.
#This has been run beforehand using the command lines
#cd /Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/simulations/phastSim
#python3 bin/phastSim --outpath /Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/simulations/ --treeFile /Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/simulations/public-latest.all.nwk --scale 0.00003344 --reference /Users/demaio/Documents/GitHub/phastSim/phastSim/example/MN908947.3.fasta --categoryProbs 0.25 0.25 0.25 0.25 --categoryRates 0.1 0.5 1.0 2.0 --outputFile phastSim_genomes_4categories_UNREST --createNewick 
#python3 bin/phastSim --outpath /Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/simulations/ --treeFile /Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/simulations/public-latest.all.nwk --scale 0.00003344 --reference /Users/demaio/Documents/GitHub/phastSim/phastSim/example/MN908947.3.fasta --outputFile phastSim_genomes_UNREST --createNewick 
#python3 bin/phastSim --outpath /Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/simulations/ --treeFile /Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/simulations/public-latest.all.nwk --scale 0.00003344 --reference /Users/demaio/Documents/GitHub/phastSim/phastSim/example/MN908947.3.fasta --alpha 0.1 --outputFile phastSim_genomes_alpha_UNREST --createNewick 
#and then moving the resulting files on the cluster in the appropriate default folders below.


#To use this script, first run 
#bsub -M 20000 /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 MAPLE_benchmarking.py --createTotalData
#/hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/MAPLE_benchmarking.py --createBashScript
#sh /nfs/research/goldman/demaio/fastLK/simulationsNew/createSubsampleInputFiles.sh

#then to run MAPLE
#sh /nfs/research/goldman/demaio/fastLK/simulationsNew/submitMAPLE0.1.8.sh
#etc
#sh /nfs/research/goldman/demaio/fastLK/simulationsNew/submitRF_MAPLE0.1.8.sh
#etc
#sh /nfs/research/goldman/demaio/fastLK/simulationsNew/submitIQtreeLK_MAPLE0.1.4.sh
#etc
#sh /nfs/research/goldman/demaio/fastLK/simulationsNew/submitLK_MAPLE0.1.5.sh
#etc
#/hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/MAPLE_benchmarking.py --collectResults

parser = argparse.ArgumentParser(description='Create files for the benchmarking of MAPLE, both input files and cluster shell scripts.')
parser.add_argument('--inputRealData',default="/nfs/research/goldman/demaio/fastLK/realData/2021-03-31_unmasked_differences_reduced.txt_consensus-based.txt", help='input real data file; should contain the difference of all samples with respet to the reference.')
parser.add_argument('--inputRealDataReference',default="/nfs/research/goldman/demaio/fastLK/realData/2021-03-31_unmasked_differences_reduced.txt_consensus.fa", help='input real data file with the reference genome used.')
parser.add_argument('--pathToSimulationFolder',default="/nfs/research/goldman/demaio/fastLK/simulationsNew/", help='path to the phastSim python script.')
parser.add_argument('--inputRealTree',default="/nfs/research/goldman/demaio/fastLK/simulationsNew/public-latest.all.nwk", help='path to the (real) input tree for phastSim simulations.')
parser.add_argument('--inputSimulationReference',default="/nfs/research/goldman/demaio/fastLK/simulationsNew/MN908947.3.fasta", help='path to the reference genome to be used for phastSim simulations.')

parser.add_argument("--createTotalData", help="Extract realized trees from global ones and mask simulated alignments similarly to the real one.", action="store_true")

parser.add_argument("--createBashScript", help="Create bash script to run the file creation on the cluster in parallel.", action="store_true")

parser.add_argument("--createInputData", help="Subsample real and simulated data to create input datasets for benchmarking.", action="store_true")
parser.add_argument("--subSampleNum",help="Number of subsamples to extract.",  type=int, default=10)
parser.add_argument("--repeat",help="Which repeat to simulate, typically between 1-10.",  type=int, default=1)
#subsamples of 0) real data; 1) simulations basic scenario; 2) simulated with rate variation; 3) simulated with high (alpha) rate variation; 4) simulations with N's.
parser.add_argument("--scenario",help="Which scenario to subsample, typically between 0-4.",  type=int, default=0)

parser.add_argument("--collectResults", help="Collect the results from the benchmarking.", action="store_true")
args = parser.parse_args()


#subsample datasets
subSampleNums=[1000,2000,5000,10000,20000,50000,100000,200000,500000]
subSampleFasta=[1000,2000,5000,10000,20000]

#define binary tree structure for quickly finding samples 
class BinarySearchTree:
	def __init__(self):
		self.root = None
	def __iter__(self):
		return self.root.__iter__()
	def put(self,key):
		if self.root:
			self._put(key,self.root)
		else:
			self.root = TreeNodeSearch(key)
	def _put(self,key,currentNode):
		if key < currentNode.key:
			if currentNode.hasLeftChild():
				self._put(key,currentNode.leftChild)
			else:
				currentNode.leftChild = TreeNodeSearch(key)
		else:
			if currentNode.hasRightChild():
				self._put(key,currentNode.rightChild)
			else:
				currentNode.rightChild = TreeNodeSearch(key)
	def get(self,key):
		if self.root:
			return self._get(key,self.root)
		else:
			return False

	def _get(self,key,currentNode):
		if not currentNode:
			return False
		elif currentNode.key == key:
			return True
		elif key < currentNode.key:
			return self._get(key,currentNode.leftChild)
		else:
			return self._get(key,currentNode.rightChild)

	def __getitem__(self,key):
		return self.get(key)

class TreeNodeSearch:
	def __init__(self,key,left=None,right=None):
		self.key = key
		self.leftChild = left
		self.rightChild = right
	def hasLeftChild(self):
		return self.leftChild
	def hasRightChild(self):
		return self.rightChild

#Class defining nodes of the tree
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
def readNewick(nwFile,multipleTrees=False,dirtiness=True):
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
					node.dist=float(distStr)
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
					node.dist=float(distStr)
					distStr=""
				node=node.up
				index+=1
			else:
				name+=nwString[index]
				index+=1
		if not finished:
			print("Error, final character ; not found in newick string in file "+nwFile+".")
			exit()

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

#create newick string for a tree 
def createNewick(node):
	nextNode=node
	stringList=[]
	direction=0
	lastNode=None
	while nextNode!=None:
		if nextNode.children:
			if direction==0:
				#print("Will go into first child")
				stringList.append("(")
				lastNode=nextNode
				nextNode=nextNode.children[0]
			elif direction==1 and lastNode!=nextNode.children[-1]:
				#print("Will go into another child")
				stringList.append(",")
				childNum=0
				#print("Num children: "+str(len(nextNode.children)))
				while nextNode.children[childNum]!=lastNode:
					#print("Not child: "+str(childNum))
					childNum+=1
				lastNode=nextNode
				nextNode=nextNode.children[childNum+1]
				direction=0
			else:
				#print("Coming from last child")
				if nextNode.dist:
					stringList.append("):"+str(nextNode.dist))
				else:
					stringList.append("):"+str(0.0))
				#if nextNode.up!=None:
				direction=1
				lastNode=nextNode
					# if nextNode.up.children[0]==nextNode:
					# 	direction=1
					# else:
					# 	direction=2
				nextNode=nextNode.up
		else:
			#print("Terminal node "+nextNode.name)
			if nextNode.dist:
				stringList.append(nextNode.name+":"+str(nextNode.dist))
			else:
				stringList.append(nextNode.name+":"+str(0.0))
			#if nextNode.up!=None:
			direction=1
			lastNode=nextNode
				# if nextNode.up.children[0]==nextNode:
				# 	direction=1
				# else:
				# 	direction=2
			nextNode=nextNode.up
	stringList.append(";")
	return "".join(stringList)

#create postorder traversal list of nodes for an input tree
def postorderList(phylo):
	numNodes=0
	if not phylo.children:
		return [phylo]
	nodeList=[]
	lastNode=phylo
	node=phylo.children[0]
	while node!=phylo or lastNode!=phylo.children[-1]:
		if lastNode==node.up:
			numNodes+=1
			if node.children:
				lastNode=node
				node=node.children[0]
			else:
				nodeList.append(node)
				lastNode=node
				node=node.up
		else:
			childIndex=0
			while lastNode!=node.children[childIndex]:
				childIndex+=1
			if childIndex==len(node.children)-1:
				lastNode=node
				nodeList.append(node)
				node=node.up
			else:
				lastNode=node
				node=node.children[childIndex+1]
	nodeList.append(phylo)
	print("Numbers of nodes:")
	print(numNodes+1)
	print(len(nodeList))
	return nodeList

#Robinson-Foulds distance (1981) using a simplification of the algorithm from Day 1985.
#this function prepare the data to compare trees to a reference one t1.
#I split in two functions so that I don't have to repeat these steps for the reference tree if I compare to the same reference tree multiple times.
def prepareTreeComparison(t1,rooted=False,minimumBLen=0.000006):
	#dictionary of values given to sequence names
	leafNameDict={}
	#list of sequence names sorted according to value
	leafNameDictReverse=[]
	#table containing clusters in the tree
	nodeTable=[]
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
				if node==t1:
					nodeTable[lastR][0]=lastL
					nodeTable[lastR][1]=lastR
				else:
					if node.dist>minimumBLen:
						numBranches+=1
						if rooted or lastL>0:
							if node==node.up.children[-1]:
								nodeTable[lastL][0]=lastL
								nodeTable[lastL][1]=lastR
							else:
								nodeTable[lastR][0]=lastL
								nodeTable[lastR][1]=lastR
						else: # re-root at leaf 0, so flip the values for the current branch if it contains leaf 0.
								flippedL=lastR+1
								flippedR=nLeaves-1
								nodeTable[flippedL][0]=flippedL
								nodeTable[flippedL][1]=flippedR
			else:
				nextNode=node.children[node.exploredChildren]
				movingFrom=0
		node=nextNode
	return leafNameDict, nodeTable, leafCount, numBranches,

#Robinson-Foulds distance (1981) using a simplification of the algorithm from Day 1985.
#this function compares the current tree t2 to a previous one for which prepareTreeComparison() was run.
def RobinsonFouldsWithDay1985(t2,leafNameDict,nodeTable,leafCount,numBranches,rooted=False,minimumBLen=0.000006):
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
								if nodeTable[lastL][0]==lastL and nodeTable[lastL][1]==lastR:
									foundBranches+=1
								elif nodeTable[lastR][0]==lastL and nodeTable[lastR][1]==lastR:
									foundBranches+=1
								else:
									missedBranches+=1
							else: # re-root at leaf 0, so flip the values for the current branch if it contains leaf 0.
								flippedL=lastR+1
								flippedR=leafCount-1
								if nodeTable[flippedL][0]==flippedL and nodeTable[flippedL][1]==flippedR:
									foundBranches+=1
								elif nodeTable[flippedR][0]==flippedL and nodeTable[flippedR][1]==flippedR:
									foundBranches+=1
								else:
									missedBranches+=1
						else:
							missedBranches+=1
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
	return numDiffs, float(numDiffs)/(2*(leafCount-3)), leafCount, foundBranches, missedBranches, (numBranches-foundBranches)



def readConciseAlignment(fileName,numbersFirst=True,shift01pos=False):
	start = time.time()
	fileI=open(fileName)
	line=fileI.readline()
	nSeqs=0
	data={}
	addendum=0
	if shift01pos:
		addendum=1
	while line!="" and line!="\n":
		nSeqs+=1
		seqList=[]
		name=line.replace(">","").replace("\n","")
		line=fileI.readline()
		while line!="" and line!="\n" and line[0]!=">":
			linelist=line.split()
			if len(linelist)>2:
				print("Format problem")
				exit()
				#entry=(linelist[0],int(linelist[1]),int(linelist[2]))
			else:
				if numbersFirst:
					entry=(linelist[1],int(linelist[0])+addendum)
				else:
					entry=(linelist[0],int(linelist[1])+addendum)
			seqList.append(entry)
			line=fileI.readline()
		data[name]=seqList
	fileI.close()
	time2 = time.time() - start
	print("Time to read DNA reduced data file: "+str(time2))
	print(str(nSeqs)+" sequences in file.")
	return data

def readConciseAlignmentReal(fileName):
	start = time.time()
	fileI=open(fileName)
	line=fileI.readline()
	nSeqs=0
	data={}
	while line!="" and line!="\n":
		nSeqs+=1
		seqList=[]
		name=line.replace(">","").replace("\n","")
		line=fileI.readline()
		while line!="" and line!="\n" and line[0]!=">":
			linelist=line.split()
			if len(linelist)==3 and int(linelist[2])>1:
				entry=("N",int(linelist[1]),int(linelist[2]))
			else:
				if linelist[0]=="-":
					entry=("N",int(linelist[1]))
				else:
					entry=(linelist[0].upper(),int(linelist[1]))
			seqList.append(entry)
			line=fileI.readline()
		data[name]=seqList
	fileI.close()
	time2 = time.time() - start
	print("Time to read DNA reduced data file: "+str(time2))
	print(str(nSeqs)+" sequences in file.")
	return data

#assign ambiguities from leafAmb to tip "leaf"
def applyAmbiguities(leaf,leafAmb):
	#count number of isolated ambiguities and create subvector of leafAmb with just long stretches of "N"s.
	onlyAmb=[]
	numAmb=0
	lastPos=1
	for entry in leafAmb:
		if (len(entry)==2 or entry[2]==1) and (entry[0]=="N" or (not entry[0] in allelesList)):
			numAmb+=1
		elif entry[0]=="N" and len(entry)>2 and entry[2]>1:
			if entry[1]>lastPos:
				onlyAmb.append(["R",lastPos,entry[1]-lastPos])
			onlyAmb.append(entry)
			lastPos=entry[1]+entry[2]
	if lastPos<=len(ref)-12:
		onlyAmb.append(["R",lastPos,len(ref)-12+1-lastPos])
	#count number of non-reference entries in "leaf", and sample those to me masked.
	numDiff=0
	lastPos=1
	leafR=[]
	for entry in leaf:
		if entry[1]<=len(ref)-12:
			if (len(entry)<3 or entry[2]==1) and entry[0]!="N":
				numDiff+=1
			if entry[1]>lastPos:
				leafR.append(["R",lastPos,entry[1]-lastPos])
			leafR.append(entry)
			if len(entry)==3:
				lastPos=entry[1]+entry[2]
			else:
				lastPos=entry[1]+1
	if lastPos<=len(ref)-12:
		leafR.append(["R",lastPos,len(ref)-12+1-lastPos])
	if numAmb>=numDiff:
		masked=range(numAmb)
	else:
		masked=random.sample(range(numDiff), numAmb)

	#now create "newLeaf" by masking "leaf"
	indexEntry1=0
	indexEntry2=0
	pos=1
	entry1=leafR[indexEntry1]
	pos1=entry1[1]
	if entry1[0]!="N" and entry1[0]!="R":
		end1=pos1
	else:
		end1=pos1+entry1[2]-1

	entry2=onlyAmb[indexEntry2]
	pos2=entry2[1]
	end2=pos2+entry2[2]-1
	
	end=min(end1,end2)
	length=end+1-pos
	newLeaf=[]
	diffIndex=-1
	while True:
		if entry1[0]!="N" and entry1[0]!="O" and entry1[0]!="R":
			diffIndex+=1
			if diffIndex in masked or entry2[0]=="N":
				newLeaf.append(["N",pos,length])
			else:
				newLeaf.append([entry1[0],pos,length])
		elif entry2[0]=="N":
			newLeaf.append(["N",pos,length])
		else:
			newLeaf.append([entry1[0],pos,length])

		pos+=length
		if pos>lRef-12:
			break
		if pos>end1:
			indexEntry1+=1
			entry1=leafR[indexEntry1]
			pos1=entry1[1]
			if entry1[0]!="N" and entry1[0]!="R":
				end1=pos1
			else:
				end1=pos1+entry1[2]-1
		if pos>end2:
			indexEntry2+=1
			entry2=onlyAmb[indexEntry2]
			pos2=entry2[1]
			if entry2[0]!="N" and entry2[0]!="R":
				end2=pos2
			else:
				end2=pos2+entry2[2]-1
		end=min(end1,end2)
		length=end+1-pos

	return newLeaf

alleles={"A":0,"C":1,"G":2,"T":3}
allelesList=["A","C","G","T"]
allelesLow={"a":0,"c":1,"g":2,"t":3}
allelesListLow=["a","c","g","t"]
ambiguities={"y":[0.0,1.0,0.0,1.0],"r":[1.0,0.0,1.0,0.0],"w":[1.0,0.0,0.0,1.0],"s":[0.0,1.0,1.0,0.0],"k":[0.0,0.0,1.0,1.0],"m":[1.0,1.0,0.0,0.0],"d":[1.0,0.0,1.0,1.0],"v":[1.0,1.0,1.0,0.0],"h":[1.0,1.0,0.0,1.0],"b":[0.0,1.0,1.0,1.0]}

#collect reference
def collectReference(fileName):
	file=open(fileName)
	line=file.readline()
	ref=""
	while line!="":
		line=file.readline()
		ref+=line.replace("\n","")
	#lRef=len(ref)
	#print("Ref genome length: "+str(lRef))
	file.close()
	return ref

ref = collectReference(args.inputRealDataReference)
lRef=len(ref)
simuRef = collectReference(args.inputSimulationReference)
lSimuRef=len(simuRef)

#create the input data for the benchmarking analysis
if args.createTotalData:
	random.seed(a=1)
	alignments=[args.pathToSimulationFolder+"phastSim_genomes_UNREST.txt",args.pathToSimulationFolder+"phastSim_genomes_4categories_UNREST.txt",args.pathToSimulationFolder+"phastSim_genomes_alpha_UNREST.txt"]
	#create realized trees, collapsing branches without mutations
	for alignment in alignments:
		print("Preparing mutation-informed simulated global tree from "+alignment)
		file=open(alignment.replace("txt","tree"))
		treeLine=file.readline()
		treeLine=treeLine.replace("&mutations={}","")
		charList=[]
		character=treeLine[0]
		index=0
		while character!="\n":
			if character=="[":
				index+=1
				if treeLine[index]=="]":
					while treeLine[index]!="," and treeLine[index]!=")" and treeLine[index]!=";":
						index+=1
					character=treeLine[index]
					charList.append(":0.0")
				else:
					while treeLine[index]!="]":
						index+=1
					index+=1
					character=treeLine[index]
			else:
				charList.append(character)
				index+=1
				character=treeLine[index]
		charList.append("\n")
		newTreeString="".join(charList)
		file=open(alignment.replace(".txt","_realizedTree.tree"),"w")
		file.write(newTreeString)
		file.close()
		print("Extracted realized tree")

	#apply N's to the simulated data 
	diffFile=(alignments[0].replace(".txt",""))+"_Ns.txt"
	data=readConciseAlignment(alignments[0],shift01pos=True)
	dataReal=readConciseAlignmentReal(args.inputRealData)
	samples=data.keys()
	print(str(len(samples))+" sequences in the concise DNA data file; creating simulated alignment with N's.")
	#Now create version of the data with ambiguous characters
	keys2=list(dataReal.keys())
	count=0
	fileO=open(diffFile,"w")
	for name1 in samples:
		fileO.write(">"+name1+"\n")
		i2=random.randint(0,len(keys2)-1)
		name2=keys2[i2]
		newLeaf=applyAmbiguities(data[name1],dataReal[name2])
		for m in newLeaf:
			if m[0]!="R":
				if len(m)==2:
					fileO.write(m[0]+"\t"+str(m[1])+"\n")
				else:
					fileO.write(m[0]+"\t"+str(m[1])+"\t"+str(m[2])+"\n")
		count+=1
		if count==1:
			print(data[name1])
			print(dataReal[name2])
			print(newLeaf)
		data[name1]=None
		if (count%100000)==0:
			print(count)
	fileO.close()
	del dataReal
	del data



if args.createBashScript:
	#creating bash script to run this python script in parallel on the cluster
	file=open(args.pathToSimulationFolder+"createSubsampleInputFiles.sh","w")
	folders=["realDataSubsamples","simulationsSubsamples","simulations4catSubsamples","simulationsAlphaSubsamples","simulationsNsSubsamples"]
	for j in [1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000]:
		file.write("for i in $(seq 1 10)\n"+"do \n\t")
		for scenario in range(len(folders)):
			folderNameSimu=args.pathToSimulationFolder+folders[scenario]+'/'+str(j)+"subsamples/"
			file.write("bsub -M "+str(int(15000+j/20))+" -o "+folderNameSimu+"repl\"$i\"_fileCreation_console_output.txt -e "+folderNameSimu+"repl\"$i\"_fileCreation_console_error.txt"
			+" /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/MAPLE_benchmarking.py --createInputData --scenario "+str(scenario)+" --subSampleNum "+str(j)+" --repeat \"$i\" \n\n\t")
		file.write("done\n\n")
	file.close()
	print("Created bash file "+args.pathToSimulationFolder+"createSubsampleInputFiles.sh")

	#creating folders for output files
	for scenario in range(len(folders)):
		folder=folders[scenario]
		for j in [1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000]:
			for i in range(10):
				folderNameSimu=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/output_repl"+str(i+1)+"/"
				if not os.path.isdir(folderNameSimu):
					os.mkdir(folderNameSimu)

	




#create subsample files 
if args.createInputData:
	scenario=args.scenario
	subsampleTreeInference=args.subSampleNum
	#subsamples of 1) real data; 2) simulations basic scenario; 3) simulated with rate variation; 4) simulated with high (alpha) rate variation; 5) simulations with N's.
	folders=["realDataSubsamples","simulationsSubsamples","simulations4catSubsamples","simulationsAlphaSubsamples","simulationsNsSubsamples"]
	alignments=[args.inputRealData, args.pathToSimulationFolder+"phastSim_genomes_UNREST.txt", args.pathToSimulationFolder+"phastSim_genomes_4categories_UNREST.txt", args.pathToSimulationFolder+"phastSim_genomes_alpha_UNREST.txt", args.pathToSimulationFolder+"phastSim_genomes_UNREST_Ns.txt"]
	phylos=["",args.pathToSimulationFolder+"phastSim_genomes_UNREST_realizedTree.tree", args.pathToSimulationFolder+"phastSim_genomes_4categories_UNREST_realizedTree.tree", args.pathToSimulationFolder+"phastSim_genomes_alpha_UNREST_realizedTree.tree", args.pathToSimulationFolder+"phastSim_genomes_UNREST_realizedTree.tree"]
	seed=args.repeat
	repeat=seed

	if scenario==0:
		ref = collectReference(args.inputRealDataReference)
		lRef=len(ref)
	else:
		ref = collectReference(args.inputSimulationReference)
		lRef=len(simuRef)
	#for scenario in range(len(folders)):
	folder=folders[scenario]
	pathToRepeat=args.pathToSimulationFolder+folder+'/'+str(subsampleTreeInference)+"subsamples/repeat"+str(repeat)+"_"+str(subsampleTreeInference)+"samples_"+folders[scenario]
	print("Considering scenario "+folder)
	if not os.path.isdir(args.pathToSimulationFolder+folder):
		os.mkdir(args.pathToSimulationFolder+folder)
	if scenario==0 or scenario==4:
		data=readConciseAlignmentReal(alignments[scenario])
	else:
		data=readConciseAlignment(alignments[scenario],shift01pos=True)
	samples=data.keys()
	print(str(len(samples))+" sequences in the total concise DNA data file")

	#for subsampleTreeInference in subSampleNums:
	print("Creating subsample files for "+str(subsampleTreeInference)+" subsamples")
	if not os.path.isdir(args.pathToSimulationFolder+folder+'/'+str(subsampleTreeInference)+"subsamples/"):
		os.mkdir(args.pathToSimulationFolder+folder+'/'+str(subsampleTreeInference)+"subsamples/")
	subsample=subsampleTreeInference
	random.seed(seed)
	newSamples=random.sample(samples,subsampleTreeInference)

	#create subsampled MAPLE file
	diffFile=pathToRepeat+".txt"
	fileO=open(diffFile,"w")
	numSeq=0
	for s in newSamples:
		numSeq+=1
		if numSeq==1:
			name1=s
			print(name1)
		if numSeq==2:
			name2=s
			print(name2)
		fileO.write(">"+s+"\n")
		for m in data[s]:
			if m[0]!="R":
				if len(m)==2:
					fileO.write(m[0]+"\t"+str(m[1])+"\n")
				else:
					fileO.write(m[0]+"\t"+str(m[1])+"\t"+str(m[2])+"\n")
		data[s]=None
	fileO.close()
	print("MAPLE file created")
	del data
	del samples

	#create fasta and phylip files for subsample
	#if subsampleTreeInference in subSampleFasta:
	data=readConciseAlignmentReal(diffFile)
	newSamples=data.keys()
	phylipFile=pathToRepeat+".phy"
	fastaFile=pathToRepeat+".fa"
	fileFa=open(fastaFile,"w")
	fileO=open(phylipFile,"w")
	lRef=len(ref)
	fileO.write(str(subsample)+"\t"+str(lRef)+"\n")
	for s in newSamples:
		fileO.write(s+" ")
		fileFa.write(">"+s+"\n")
		refList=list(ref)
		for m in data[s]:
			if m[0]!="R":
				if len(m)==2:
					refList[m[1]-1]=m[0]
				else:
					for i in range(m[2]):
						refList[m[1]+i-1]=m[0]
		data[s]=None
		seq="".join(refList)
		fileO.write(seq+"\n")
		fileFa.write(seq+"\n")
	fileO.close()
	fileFa.close()
	print("Fasta file created")
	#del data

	#create initial (dummy) tree for USheER
	file=open(pathToRepeat+"_initialTree.tree","w")
	file.write("("+name1+":10,"+name2+":10):1;\n")
	file.close()
	#create VCN file for USheER
	newFastaFile=pathToRepeat+"_withRef.fa"
	vcfFile=pathToRepeat+".vcf"
	if scenario==0:
		os.system("cat "+args.inputRealDataReference+" "+fastaFile+" > "+newFastaFile+"\n")
	else:
		os.system("cat "+args.inputSimulationReference+" "+fastaFile+" > "+newFastaFile+"\n")
	os.system("/hps/software/users/goldman/nicola/miniconda3/envs/usher-env/bin/faToVcf "+newFastaFile+" "+vcfFile+"\n")
	print("VCF file created")

	#extract subtree of the larger simulated tree
	if scenario>0:
		#create search tree for current sample names:
		print("Total number of leaves to be found: "+str(len(newSamples)))
		binTree=BinarySearchTree()
		for s in newSamples:
			binTree.put(s)

		phylo = readNewick(phylos[scenario])[0]
		nodeList=postorderList(phylo)

		#extract subtree for subsample
		numNode=0
		numDescList=[0,0,0,0,0]
		for node in nodeList:
			numNode+=1
			#print("numNode "+str(numNode))
			if len(node.children)==0:
				numDescList[0]+=1
				#print(node.name)
				if binTree[node.name]:
					node.subtree=Tree(name=node.name,dist=node.dist)
					#if len(newSamples)<15:
					#	print("found "+node.name)
				else:
					node.subtree=None
			else:
				if len(node.children)<4:
					numDescList[len(node.children)]+=1
				else:
					numDescList[-1]+=1
				#print(len(node.children))
				numDesc=0
				for c in node.children:
					if c.subtree!=None:
						numDesc+=1
						child=c
				if numDesc==0:
					node.subtree=None
				elif numDesc==1:
					node.subtree=child.subtree
					node.subtree.dist+=node.dist
				else:
					node.subtree=Tree(dist=node.dist)
					for c in node.children:
						if c.subtree:
							node.subtree.add_child(c.subtree)
							c.subtree.up=node.subtree
							#node.subtree.children[-1].dist=c.subtree.dist

		#print("Distribution of number of descendants:")
		#print(numDescList)
		#print("num nodes explored: "+str(numNode))
		#print("root "+str(numDesc))
		#print("Replicate "+str(repeat))
		#print("num samples "+str(subsampleTreeInference))
		# if len(newSamples)<15:
		# 	print(newSamples)
		# 	print(len(phylo.subtree.children))
		# 	nodeList=postorderList(phylo.subtree)
		# 	print(nodeList[0].name)
		# 	print(nodeList[1].name)
		# 	print(len(nodeList))
		# 	print(pathToRepeat+"_realized.nw")
		newickString=createNewick(phylo.subtree)
		#if len(newSamples)<15:
		#	print(newickString)
		phylo=None
		#del newSamples
		del data
		phyloFile=open(pathToRepeat+"_realized.nw","w")
		phyloFile.write(newickString+"\n")
		phyloFile.close()
		print("Subtree extracted and written to file "+pathToRepeat+"_realized.nw")













if args.createBashScript:
	numSamples=[1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000]
	numSamples=[10000]
	MAPLEversions=["0.1.4","0.1.5","0.1.6","0.1.7","0.1.9"]
	MAPLEversions=["0.1.9"]
	testOptions=True

	#creating bash script to run this python script in parallel on the cluster to create subsampled data
	file=open(args.pathToSimulationFolder+"createSubsampleInputFiles.sh","w")
	folders=["realDataSubsamples","simulationsSubsamples","simulations4catSubsamples","simulationsAlphaSubsamples","simulationsNsSubsamples"]
	for j in numSamples:
		file.write("for i in $(seq 1 10)\n"+"do \n\t")
		for scenario in range(len(folders)):
			folderNameSimu=args.pathToSimulationFolder+folders[scenario]+'/'+str(j)+"subsamples/"
			file.write("bsub -M "+str(int(15000+j/20))+" -o "+folderNameSimu+"repl\"$i\"_fileCreation_console_output.txt -e "+folderNameSimu+"repl\"$i\"_fileCreation_console_error.txt"
			+" /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/MAPLE_benchmarking.py --createInputData --scenario "+str(scenario)+" --subSampleNum "+str(j)+" --repeat \"$i\" \n\n\t")
		file.write("done\n\n")
	file.close()

	#creating bash scripts for MAPLE
	MAPLEoptions=[""," --fast"," --rateVariation"," --fast --rateVariation"," --model UNREST"," --fast --model UNREST"," --model UNREST --rateVariation"," --fast --model UNREST --rateVariation"]
	MAPLEoptionsNames=["","_fast","_rateVar","_fast_rateVar","_unrest","_fast_unrest","_unrest_rateVar","_fast_unrest_rateVar"]
	foldersForTreeFile=["","simulationsSubsamples","simulations4catSubsamples","simulationsAlphaSubsamples","simulationsSubsamples"]
	for versNum in range(len(MAPLEversions)):
		version=MAPLEversions[versNum]
		file=open(args.pathToSimulationFolder+"submitMAPLE"+version+".sh","w")
		for j in numSamples:
			file.write("for i in $(seq 1 10)\n do\n\n")
			for scenario in range(len(folders)):
				if scenario==0:
					refFile=args.inputRealDataReference
				else:
					refFile=args.inputSimulationReference
				folder=folders[scenario]
				folderNameSimu=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/output_repl\"$i\"/"
				pathToFile=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/repeat\"$i\"_"+str(j)+"samples_"+folder

				for option in range(len(MAPLEoptions)):
					if option<2 or testOptions:
						#remove existing result if already exists - this helps in case of re-running an analysis, and in case the new version fails, to wrongly use the old tree as estimate of the new version.
						file.write("rm -f "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_tree.tree || true\n")
						#run MAPLE
						file.write("bsub -M "+str(int(400+j/20))+" -o "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_console_output.txt -e "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_console_error.txt "
						+"/hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/MAPLEv"+version+".py --reference "+refFile+MAPLEoptions[option]+" --input "
						+pathToFile+".txt --overwrite --output "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"\n")

			file.write("done\n\n")
		file.close()
		print("Created MAPLE bash script "+args.pathToSimulationFolder+"submitMAPLE"+version+".sh")

		#create bash script for running robinson-foulds distance estimations
		versionForRF="0.1.8"
		file=open(args.pathToSimulationFolder+"submitRF_MAPLE"+version+".sh","w")
		for j in numSamples:
			file.write("for i in $(seq 1 10)\n do\n\n")
			for scenario in range(len(folders)):
				if scenario>0:
					refFile=args.inputSimulationReference
					folder=folders[scenario]
					folderNameSimu=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/output_repl\"$i\"/"
					pathToFile=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/repeat\"$i\"_"+str(j)+"samples_"+folder
					pathToTreeFile=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/repeat\"$i\"_"+str(j)+"samples_"+folder
					#pathToTreeFile=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/repeat\"$i\"_"+str(j)+"samples_"+foldersForTreeFile[scenario]
					#/nfs/research/goldman/demaio/fastLK/simulationsNew/simulationsNsSubsamples/1000subsamples/repeat6_1000samples_simulationsSubsamples_realized.nw
					for option in range(len(MAPLEoptions)):
						if option<2 or testOptions:
							#remove existing result if already exists - this helps in case of re-running an analysis, and in case the new version fails, to wrongly use the old tree as estimate of the new version.
							file.write("rm -f "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_RFdistances.txt || true\n")
							#run MAPLE RF alculation
							file.write("bsub -M "+str(int(400+j/20))+" -o "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_RF_console_output.txt -e "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_RF_console_error.txt "
							+"/hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/MAPLEv"+versionForRF+".py --inputRFtrees "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_tree.tree --inputTree "+pathToTreeFile+"_realized.nw "
							+" --overwrite --output "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+" --genomeLength "+ lRef +" --errorTesting True "+"\n") # todo; where to get the genomeLength from?

			file.write("done\n\n")
		file.close()
		print("Created RF bash script "+args.pathToSimulationFolder+"submitRF_MAPLE"+version+".sh")

		#create bash script for running likelihood evaluations (in MAPLE)
		versionForLK="0.1.9"
		LKoptions=["",""," --rateVariation"," --rateVariation",""]
		furtherLKoptions=[""," --model UNREST"]
		furtherLKoptionsNames=["","_unrest"]
		file=open(args.pathToSimulationFolder+"submitLK_MAPLE"+version+".sh","w")
		for j in numSamples:
			file.write("for i in $(seq 1 10)\n do\n\n")
			for scenario in range(len(folders)):
				if scenario==0:
					refFile=args.inputRealDataReference
				else:
					refFile=args.inputSimulationReference
				folder=folders[scenario]
				folderNameSimu=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/output_repl\"$i\"/"
				pathToFile=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/repeat\"$i\"_"+str(j)+"samples_"+folder

				for option in range(len(MAPLEoptions)):
					if option<2 or testOptions:
						for furtherOptionNum in range(len(furtherLKoptions)):
							#remove existing result if already exists - this helps in case of re-running an analysis, and in case the new version fails, to wrongly use the old tree as estimate of the new version.
							file.write("rm -f "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_LKestimation"+furtherLKoptionsNames[furtherOptionNum]+"_LK.txt || true\n")
							#run MAPLE LK calculation
							file.write("bsub -M "+str(int(400+j/20))+" -o "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_LK_console_output.txt -e "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_LK_console_error.txt "
							+"/hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/MAPLEv"+versionForLK+".py --reference "+refFile+LKoptions[scenario]+furtherLKoptions[furtherOptionNum]+" --inputTree "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_tree.tree --input "
							+pathToFile+".txt --calculateLKfinalTree --numTopologyImprovements 0 --noFastTopologyInitialSearch  --overwrite --output "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_LKestimation"+furtherLKoptionsNames[furtherOptionNum]+"\n")

							#for real data evaluate also under ratevariation
							if scenario==0:
								#remove existing result if already exists - this helps in case of re-running an analysis, and in case the new version fails, to wrongly use the old tree as estimate of the new version.
								file.write("rm -f "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_LKestimation"+furtherLKoptionsNames[furtherOptionNum]+"_rateVar_LK.txt || true\n")
								#run MAPLE LK calculation
								file.write("bsub -M "+str(int(400+j/20))+" -o "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_LK_console_output.txt -e "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_LK_console_error.txt "
								+"/hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/MAPLEv"+versionForLK+".py --reference "+refFile+LKoptions[scenario]+furtherLKoptions[furtherOptionNum]+" --rateVariation --inputTree "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_tree.tree --input "
								+pathToFile+".txt --calculateLKfinalTree --numTopologyImprovements 0 --noFastTopologyInitialSearch  --overwrite --output "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_LKestimation"+furtherLKoptionsNames[furtherOptionNum]+"_rateVar\n")

			file.write("done\n\n")
		file.close()
		print("Created LK bash script "+args.pathToSimulationFolder+"submitLK_MAPLE"+version+".sh")




		#create bash script for running likelihood evaluations (in IQtree)
		folders=["realDataSubsamples","simulationsSubsamples","simulations4catSubsamples","simulationsAlphaSubsamples","simulationsNsSubsamples"]
		LKoptions=["","","+G","+G",""]
		file=open(args.pathToSimulationFolder+"submitIQtreeLK_MAPLE"+version+".sh","w")
		for j in [1000, 2000, 5000, 10000, 20000]:
			file.write("for i in $(seq 1 10)\n do\n\n")
			for scenario in range(len(folders)):
				folder=folders[scenario]
				folderNameSimu=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/output_repl\"$i\"/"
				pathToFile=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/repeat\"$i\"_"+str(j)+"samples_"+folder

				for option in range(len(MAPLEoptions)):
					if option<2 or testOptions:
						if scenario==2 or scenario==3:
							memory=int(500+j*3.0)
						else:
							memory=int(500+j*1.7)
						file.write("cd "+folderNameSimu+"\n\t"+"bsub -M "+str(memory)+" -o "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_IQtreeLK_console_output.txt -e "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_IQtreeLK_console_error.txt "
						+"/hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+pathToFile+".phy -st DNA -te "+folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_tree.tree -pre MAPLE"+version+MAPLEoptionsNames[option]+"_IQtreeLKestimation -m GTR"+LKoptions[scenario]+" -quiet -nt 1 -keep-ident -redo -blmin 0.000000005 \n")

			file.write("done\n\n")
		file.close()
		print("Created IQtreeLK bash script "+args.pathToSimulationFolder+"submitIQtreeLK_MAPLE"+version+".sh")




if args.collectResults:
	printOptions=True
	#MAPLEversions=["0.1.4","0.1.5","0.1.6","0.1.7","0.1.8"]
	MAPLEversions=["0.1.7"]
	#numSamples=[1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000]
	numSamples=[100000]

	#collecting the results from the benchmarking
	folders=["realDataSubsamples","simulationsSubsamples","simulations4catSubsamples","simulationsAlphaSubsamples","simulationsNsSubsamples"]
	MAPLEoptionsNames=["","_fast","_rateVar","_fast_rateVar","_unrest","_fast_unrest","_unrest_rateVar","_fast_unrest_rateVar"]
	times={}
	memories={}
	maxMemories={}
	lks={}
	lksUnrest={}
	lksSiteVar={}
	lksSiteVarUnrest={}
	IQlks={}
	RFs={}
	for versNum in range(len(MAPLEversions)):
		version=MAPLEversions[versNum]
		times[versNum]={}
		memories[versNum]={}
		maxMemories[versNum]={}
		lks[versNum]={}
		lksUnrest[versNum]={}
		lksSiteVar[versNum]={}
		lksSiteVarUnrest[versNum]={}
		IQlks[versNum]={}
		RFs[versNum]={}
		for j in numSamples:
			times[versNum][j]={}
			memories[versNum][j]={}
			maxMemories[versNum][j]={}
			lks[versNum][j]={}
			lksUnrest[versNum][j]={}
			lksSiteVar[versNum][j]={}
			lksSiteVarUnrest[versNum][j]={}
			IQlks[versNum][j]={}
			RFs[versNum][j]={}
			for scenario in range(len(folders)):
				times[versNum][j][scenario]={}
				memories[versNum][j][scenario]={}
				maxMemories[versNum][j][scenario]={}
				lks[versNum][j][scenario]={}
				lksUnrest[versNum][j][scenario]={}
				lksSiteVar[versNum][j][scenario]={}
				lksSiteVarUnrest[versNum][j][scenario]={}
				IQlks[versNum][j][scenario]={}
				RFs[versNum][j][scenario]={}
				folder=folders[scenario]
				for option in range(len(MAPLEoptionsNames)):
					times[versNum][j][scenario][option]=[]
					memories[versNum][j][scenario][option]=[]
					maxMemories[versNum][j][scenario][option]=[]
					lks[versNum][j][scenario][option]=[]
					lksUnrest[versNum][j][scenario][option]=[]
					lksSiteVar[versNum][j][scenario][option]=[]
					lksSiteVarUnrest[versNum][j][scenario][option]=[]
					IQlks[versNum][j][scenario][option]=[]
					RFs[versNum][j][scenario][option]=[]
					if option<2 or printOptions:
						notFound=""
						notFoundCount=0
						notFoundTime=""
						notFoundTimeCount=0
						notFoundRF=""
						notFoundRFCount=0
						notFoundLK=""
						notFoundLKCount=0
						notFoundLKunrest=""
						notFoundLKunrestCount=0
						notFoundLKsiteVarUnrest=""
						notFoundLKsiteVarUnrestCount=0
						notFoundLKsiteVar=""
						notFoundLKsiteVarCount=0
						notFoundIQLK=""
						notFoundIQLKCount=0
						timeFile=open(args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/MAPLE"+version+MAPLEoptionsNames[option]+"_times.txt","w")
						memoryFile=open(args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/MAPLE"+version+MAPLEoptionsNames[option]+"_memory.txt","w")
						maxMemoryFile=open(args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/MAPLE"+version+MAPLEoptionsNames[option]+"_maxMemory.txt","w")
						lkFile=open(args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/MAPLE"+version+MAPLEoptionsNames[option]+"_lk.txt","w")
						lkUnrestFile=open(args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/MAPLE"+version+MAPLEoptionsNames[option]+"_lkUnrest.txt","w")
						if scenario==0:
							lkSiteVarUnrestFile=open(args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/MAPLE"+version+MAPLEoptionsNames[option]+"_lkSiteVarUnrest.txt","w")
							lkSiteVarFile=open(args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/MAPLE"+version+MAPLEoptionsNames[option]+"_lkSiteVar.txt","w")
						iqtreeLKFile=open(args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/MAPLE"+version+MAPLEoptionsNames[option]+"_IQtreeLK.txt","w")
						if scenario>0:
							rfFile=open(args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/MAPLE"+version+MAPLEoptionsNames[option]+"_RF.txt","w")
						for i in range(10):
							folderNameSimu=args.pathToSimulationFolder+folder+'/'+str(j)+"subsamples/output_repl"+str(i+1)+"/"
							if not os.path.isfile(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_tree.tree"):
								notFound+="-"+str(i+1)
								notFoundCount+=1
							else:

								#collect runtime and memory demands
								file=open(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_console_output.txt")
								line=file.readline()
								timeRun=-1.0
								while line!="":
									if "Final Substitution matrix:" in line:
										#bestLK=float(line.split()[1].replace(" ","").replace("\n",""))
										while line!="Resource usage summary:\n":
											line=file.readline()
										line=file.readline()
										line=file.readline()
										timeRun=float(line.split()[3])
										line=file.readline()
										#print("FastLK tree file with "+str(j)+" samples, replicate "+str(i+1)+" speed "+s+" ")
										#print(line)
										if line.split()[3]=="-":
											maxMemoryUsed=float("NaN")
											aveMemoryUsed=float("NaN")
											continue
										maxMemoryUsed=float(line.split()[3])
										line=file.readline()
										if line.split()[3]=="-":
											aveMemoryUsed=float("NaN")
											continue
										aveMemoryUsed=float(line.split()[3])
									line=file.readline()
								if timeRun<0.0:
									notFoundTime+="-"+str(i+1)
									notFoundTimeCount+=1
									#print("MAPLE"+versions[s]+" with "+str(j)+" samples "+simulationsMAPLE[simu]+", replicate "+str(i+1)+" - running time not found.")
								else:
									times[versNum][j][scenario][option].append(timeRun)
									memories[versNum][j][scenario][option].append(aveMemoryUsed)
									maxMemories[versNum][j][scenario][option].append(maxMemoryUsed)
									#lks[versNum][j][scenario][option]=[]
									#IQlks[versNum][j][scenario][option]=[]
									#lks[j][s].append(bestLK)
								file.close()

								#collect LKs
								if not os.path.isfile(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_LKestimation_LK.txt"):
									notFoundLK+="-"+str(i+1)
									notFoundLKCount+=1
								else:
									file=open(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_LKestimation_LK.txt")
									line=file.readline()
									lk=float(line.split()[0])
									lks[versNum][j][scenario][option].append(lk)
									file.close()

								#collect LKs estimated under UNREST model
								if not os.path.isfile(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_LKestimation_unrest_LK.txt"):
									notFoundLKunrest+="-"+str(i+1)
									notFoundLKunrestCount+=1
								else:
									file=open(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_LKestimation_unrest_LK.txt")
									line=file.readline()
									lk=float(line.split()[0])
									lksUnrest[versNum][j][scenario][option].append(lk)
									file.close()

								if scenario==0:
									#collect also likelihoods estimated under rate variation for the real data
									if not os.path.isfile(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_LKestimation_unrest_rateVar_LK.txt"):
										notFoundLKsiteVarUnrest+="-"+str(i+1)
										notFoundLKsiteVarUnrestCount+=1
									else:
										file=open(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_LKestimation_unrest_rateVar_LK.txt")
										line=file.readline()
										lk=float(line.split()[0])
										lksSiteVarUnrest[versNum][j][scenario][option].append(lk)
										file.close()

									if not os.path.isfile(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_LKestimation_rateVar_LK.txt"):
										notFoundLKsiteVar+="-"+str(i+1)
										notFoundLKsiteVarCount+=1
									else:
										file=open(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_LKestimation_rateVar_LK.txt")
										line=file.readline()
										lk=float(line.split()[0])
										lksSiteVar[versNum][j][scenario][option].append(lk)
										file.close()

								#collect IQtree LKs
								if not os.path.isfile(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_IQtreeLKestimation.log"):
									notFoundIQLK+="-"+str(i+1)
									notFoundIQLKCount+=1
								else:
									file=open(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_IQtreeLKestimation.log")
									line=file.readline()
									while not ("BEST SCORE FOUND :" in line) and line!="":
										line=file.readline()
									if line=="":
										notFoundIQLK+="-"+str(i+1)
										notFoundIQLKCount+=1
									else:
										IQlk=float(line.split()[4])
										IQlks[versNum][j][scenario][option].append(IQlk)
									file.close()

								#collect RF distances
								if scenario>0:
									if not os.path.isfile(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_RFdistances.txt"):
										notFoundRF+="-"+str(i+1)
										notFoundRFCount+=1
									else:
										file=open(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_RFdistances.txt")
										line=file.readline()
										line=file.readline()
										if line.split()[1]=='None':
											print(folderNameSimu+"MAPLE"+version+MAPLEoptionsNames[option]+"_RFdistances.txt")
											print(line)
											notFoundRF+="-"+str(i+1)
											notFoundRFCount+=1
										else:
											RFdistance=float(line.split()[1])
											RFs[versNum][j][scenario][option].append(RFdistance)
										file.close()

						if printOptions or option<2:
							if notFoundCount>0:
								print("Scenario "+folder+" MAPLE"+version+", "+str(j)+" samples, options "+MAPLEoptionsNames[option]+", tree file not found for replicates "+notFound)
							if notFoundTimeCount>0:
								print("Scenario "+folder+" MAPLE"+version+", "+str(j)+" samples, options "+MAPLEoptionsNames[option]+", runtime not found for replicates "+notFoundTime)
							if notFoundLKCount>0:
								print("Scenario "+folder+" MAPLE"+version+", "+str(j)+" samples, options "+MAPLEoptionsNames[option]+", LK file not found for replicates "+notFoundLK)
							if notFoundLKunrestCount>0:
								print("Scenario "+folder+" MAPLE"+version+", "+str(j)+" samples, options "+MAPLEoptionsNames[option]+", LK under unrest file not found for replicates "+notFoundLKunrest)
							if notFoundLKsiteVarUnrestCount>0:
								print("Scenario "+folder+" MAPLE"+version+", "+str(j)+" samples, options "+MAPLEoptionsNames[option]+", LK under siteVar-unrest file not found for replicates "+notFoundLKsiteVarUnrest)
							if notFoundLKsiteVarCount>0:
								print("Scenario "+folder+" MAPLE"+version+", "+str(j)+" samples, options "+MAPLEoptionsNames[option]+", LK under siteVar file not found for replicates "+notFoundLKsiteVar)
							if notFoundIQLKCount>0:
								print("Scenario "+folder+" MAPLE"+version+", "+str(j)+" samples, options "+MAPLEoptionsNames[option]+", IQtree LK file not found for replicates "+notFoundIQLK)
							if scenario>0 and notFoundRFCount>0:
								print("Scenario "+folder+" MAPLE"+version+", "+str(j)+" samples, options "+MAPLEoptionsNames[option]+", RF file not found for replicates "+notFoundRF)
							if len(times[versNum][j][scenario][option])>0:
								print("Scenario "+folder+" MAPLE"+version+", "+str(j)+" samples, options "+MAPLEoptionsNames[option]+" : times, memory, max memory, LKs, unrest  LKs, IQtree LKs, RFdistances")
								print(times[versNum][j][scenario][option])
								print(str(min(times[versNum][j][scenario][option]))+" \t"+str(sum(times[versNum][j][scenario][option])/len(times[versNum][j][scenario][option]))+" \t"+str(max(times[versNum][j][scenario][option])))
								print(memories[versNum][j][scenario][option])
								print(str(min(memories[versNum][j][scenario][option]))+" \t"+str(sum(memories[versNum][j][scenario][option])/len(memories[versNum][j][scenario][option]))+" \t"+str(max(memories[versNum][j][scenario][option])))
								print(maxMemories[versNum][j][scenario][option])
								print(str(min(maxMemories[versNum][j][scenario][option]))+" \t"+str(sum(maxMemories[versNum][j][scenario][option])/len(maxMemories[versNum][j][scenario][option]))+" \t"+str(max(maxMemories[versNum][j][scenario][option])))
							if len(lks[versNum][j][scenario][option])>0:
								print(lks[versNum][j][scenario][option])
								print(str(min(lks[versNum][j][scenario][option]))+" \t"+str(sum(lks[versNum][j][scenario][option])/len(lks[versNum][j][scenario][option]))+" \t"+str(max(lks[versNum][j][scenario][option])))
							if len(lksUnrest[versNum][j][scenario][option])>0:
								print(lksUnrest[versNum][j][scenario][option])
								print(str(min(lksUnrest[versNum][j][scenario][option]))+" \t"+str(sum(lksUnrest[versNum][j][scenario][option])/len(lksUnrest[versNum][j][scenario][option]))+" \t"+str(max(lksUnrest[versNum][j][scenario][option])))
							if scenario==0:
								if len(lksSiteVarUnrest[versNum][j][scenario][option])>0:
									print(lksSiteVarUnrest[versNum][j][scenario][option])
									print(str(min(lksSiteVarUnrest[versNum][j][scenario][option]))+" \t"+str(sum(lksSiteVarUnrest[versNum][j][scenario][option])/len(lksSiteVarUnrest[versNum][j][scenario][option]))+" \t"+str(max(lksSiteVarUnrest[versNum][j][scenario][option])))
								if len(lksSiteVar[versNum][j][scenario][option])>0:
									print(lksSiteVar[versNum][j][scenario][option])
									print(str(min(lksSiteVar[versNum][j][scenario][option]))+" \t"+str(sum(lksSiteVar[versNum][j][scenario][option])/len(lksSiteVar[versNum][j][scenario][option]))+" \t"+str(max(lksSiteVar[versNum][j][scenario][option])))
							if len(IQlks[versNum][j][scenario][option])>0:
								print(IQlks[versNum][j][scenario][option])
								print(str(min(IQlks[versNum][j][scenario][option]))+" \t"+str(sum(IQlks[versNum][j][scenario][option])/len(IQlks[versNum][j][scenario][option]))+" \t"+str(max(IQlks[versNum][j][scenario][option])))
							if scenario>0 and len(RFs[versNum][j][scenario][option])>0:
								print(RFs[versNum][j][scenario][option])
								print(str(min(RFs[versNum][j][scenario][option]))+" \t"+str(sum(RFs[versNum][j][scenario][option])/len(RFs[versNum][j][scenario][option]))+" \t"+str(max(RFs[versNum][j][scenario][option])))
							print("\n")
						for time in times[versNum][j][scenario][option]:
							timeFile.write(str(time)+"\t")
						#timeFile.write(times[versNum][j][scenario][option])
						timeFile.close()
						for time in memories[versNum][j][scenario][option]:
							memoryFile.write(str(time)+"\t")
						#memoryFile.write(memories[versNum][j][scenario][option])
						memoryFile.close()
						for time in maxMemories[versNum][j][scenario][option]:
							maxMemoryFile.write(str(time)+"\t")
						#maxMemoryFile.write(maxMemories[versNum][j][scenario][option])
						maxMemoryFile.close()
						for time in lks[versNum][j][scenario][option]:
							lkFile.write(str(time)+"\t")
						#lkFile.write(lks[versNum][j][scenario][option])
						lkFile.close()
						for time in lksUnrest[versNum][j][scenario][option]:
							lkUnrestFile.write(str(time)+"\t")
						lkUnrestFile.close()
						if scenario==0:
							for time in lksSiteVarUnrest[versNum][j][scenario][option]:
								lkSiteVarUnrestFile.write(str(time)+"\t")
							lkSiteVarUnrestFile.close()
							for time in lksSiteVar[versNum][j][scenario][option]:
								lkSiteVarFile.write(str(time)+"\t")
							lkSiteVarFile.close()
						for time in IQlks[versNum][j][scenario][option]:
							iqtreeLKFile.write(str(time)+"\t")
						#iqtreeLKFile.write(IQlks[versNum][j][scenario][option])
						iqtreeLKFile.close()
						if scenario>0:
							for time in RFs[versNum][j][scenario][option]:
								rfFile.write(str(time)+"\t")
							#rfFile.write(RFs[versNum][j][scenario][option])
							rfFile.close()
						
						





















exit()

nRepeats=10
nRepeatsSimu=10
folderName="/nfs/research/goldman/demaio/fastLK/realData/subsamples/"
folderNameSimu="/nfs/research/goldman/demaio/fastLK/simulations/subsamples/" # folderNameSimu = args.pathToSimulationFolder + folder + '/' + str(j) + "subsamples/output_repl\"$i\"/"
speeds=["slowest","slow","medium","fast","fastest"]
allowedFails=[5,5,5,4,3]
thresholdLogLKs=[120.0, 100.0, 80.0, 60.0, 40.0]
bLenAdjustments=[10,10,10,1,1]
numTopologyImprovements=[5,3,2,1,0]
allowedFailsTopology=[6,4,3,2,1]
thresholdLogLKtopology=[150.0,100.0,80.0,60.0,40.0]
thresholdTopologyPlacement=[-0.1,-0.2,-0.5,-1.0,-2.0]
bLenFactors=[4.0,4.0,3.0,2.0,1.0]
#REDUCE THESE TO TRY AND REDUCE MEMORY COST, for example using
#bLenFactors=[4.0,3.0,2.0,1.0]
#bLenFactors=[8.0,8.0,3.0,2.0]
sampleSizes = [100, 1000, 2000, 5000]#, 10000, 20000, 50000, 100000, 200000, 500000]
errorRates = [0, 0.0001, 0.0005]# args.errorRates #
readFileUser = "demaio"
writeFileUser = "myrthe"
folderNameCode ="/nfs/research/goldman/" +writeFileUser+  "/fastLK/code/" # folderNameSimu = args.pathToSimulationFolder + folder + '/' + str(j) + "subsamples/output_repl\"$i\"/"

# file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/createFoldersSimu.sh","w")
# for j in [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000]:
# 	for s in range(len(speeds)):
# 		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"mkdir "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/\n"+"done\n\n")
# 		#file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"bsub -M "+str(int(200+j/8))+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"_console_error.txt /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/estimatePhylogenyIterativeFastLK.py --reference /nfs/research/goldman/demaio/fastLK/realData/2021-03-31_unmasked_differences_reduced.txt_consensus.fa --input "+folderName+str(j)+"subsamples/diffFile_seed\"$i\"_"+str(j)+"samples.txt --binaryTree --overwrite --allowedFails "+str(allowedFails[s])+" --thresholdLogLK "+str(thresholdLogLKs[s])+" --bLenAdjustment "+str(bLenAdjustments[s])+" --bLenFactor "+str(bLenFactors[s])+" --numTopologyImprovements "+str(numTopologyImprovements[s])+" --thresholdTopologyPlacement "+str(thresholdTopologyPlacement[s])+" --allowedFailsTopology "+str(allowedFailsTopology[s])+" --thresholdLogLKtopology "+str(thresholdLogLKtopology[s])+" --output "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"\n"+"done\n\n")
# 	file.write("\n\n")
# file.close()


"""SIMULATE ERRORS""" #Myrthe
console_output = folderNameCode+"results/console_output.txt"
console_errors = folderNameCode+"results/console_errors.txt"
file= open("simulateErrors.sh","w") #open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/createSmallMAPLEfiles.sh","w")
# file.write("mkdir " + folderNameCode + "\n")
file.write("rm -r " + console_output + "\n")# + "; mkdir " + console_output + "\n")
file.write("rm -r " + console_errors + "\n")#+"; mkdir " + console_errors + "\n")
file.write(" \n#SIMULATE ERRORS \n")
for errorRate in errorRates:
	#if errorRate:
		for j in sampleSizes:
			pathRead = "/nfs/research/goldman/"+ readFileUser + "/fastLK/simulations/subsamples/" + str(j) + 'subsamples/' #+"subsamples/"
			pathWrite = "/nfs/research/goldman/"+ writeFileUser + "/fastLK/simulations/subsamples/" + str(j) + 'subsamples/' #+"subsamples/"
			fileNameIn = "fastaFile_repeat\"$i\"_"+str(j)+"samples_Ns.fa"
			fileNameOut = "fastaFile_repeat\"$i\"_"+str(j)+"samples_Ns_errors"+ str(errorRate) + ".fa"
			file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"
				          + "mkdir -p " + pathWrite + "\n\t"
			           # perhaps make sure we only do so if it doesn't exist already.
						+"bsub -M " +str(int(400+j/2))+" -o "+ console_output + " -e " + console_errors # job creation ; output file and error file. ,
			            + " /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 "
			            +  folderNameCode + "MAPLE_simulate_errors.py " #when using python instead?
					    + " --input " + pathRead + fileNameIn
			           + " --output " + pathWrite + fileNameOut
			           + " --errorRate " + str(errorRate)+"\n done\n\n")

"""SIMULATE SITE-SPECIFIC ERRORS""" #Myrthe
console_output = folderNameCode+"results/console_output.txt"
console_errors = folderNameCode+"results/console_errors.txt"
file= open("simulateErrors.sh","w") #open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/createSmallMAPLEfiles.sh","w")
# file.write("mkdir " + folderNameCode + "\n")
file.write("rm -r " + console_output + "\n")# + "; mkdir " + console_output + "\n")
file.write("rm -r " + console_errors + "\n")#+"; mkdir " + console_errors + "\n")
file.write(" \n#SIMULATE ERRORS \n")
for errorRate in errorRates:
	if errorRate: #errorRate should not be 0.
		for j in sampleSizes:
			pathRead = "/nfs/research/goldman/"+ readFileUser + "/fastLK/simulations/subsamples/" + str(j) + 'subsamples/' #+"subsamples/"
			pathWrite = "/nfs/research/goldman/"+ writeFileUser + "/fastLK/simulations/subsamples/" + str(j) + 'subsamples/' #+"subsamples/"
			fileNameIn = "fastaFile_repeat\"$i\"_"+str(j)+"samples_Ns.fa"
			fileNameOut = "fastaFile_repeat\"$i\"_"+str(j)+"samples_Ns_errors"+ str(errorRate) + ".fa"
			file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"
				          + "mkdir -p " + pathWrite + "\n\t"
			           # perhaps make sure we only do so if it doesn't exist already.
						+"bsub -M " +str(int(400+j/2))+" -o "+ console_output + " -e " + console_errors # job creation ; output file and error file. ,
			            + " /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 "
			            +  folderNameCode + "MAPLE_simulate_errors.py " #when using python instead?
					    + " --input " + pathRead + fileNameIn
			           + " --output " + pathWrite + fileNameOut
			           + " --siteSpecific "
			           + " --errorRate " + str(errorRate)+"\n done\n\n")
#afterwarda an error rate file is stored. file[:-3] + "_siteSpecificErrors.txt","w")


"""CALL MAPLE FILE""" #Myrthe
#fixing the accidental overwriting of some small MAPLE files
file.write(" \n#CREATE MAPLE FILES \n")
for errorRate in errorRates:
	#if errorRate:
		for j in sampleSizes:
			pathWrite = "/nfs/research/goldman/"+ writeFileUser + "/fastLK/simulations/subsamples/" + str(j) + 'subsamples/' #+"subsamples/"
			fileNameIn = "fastaFile_repeat\"$i\"_"+str(j)+"samples_Ns_errors"+ str(errorRate) + ".fa"
			fileNameOut = "diffFile_repeat\"$i\"_"+str(j)+"samples_Ns_errors"+ str(errorRate) + ".txt"
			file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"
			           #+ "mkdir -p " + pathWrite + "\n\t"
			           + "bsub -M "+str(int(400+j/2))+" -o "+ console_output
				+" -e " + console_errors
				+" /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/createMapleFile.py "
				+" --path " + pathWrite
				+" --fasta " + fileNameIn
				+" --output " + fileNameOut
				+" --overwrite\n" #Myrthe: added '_errors' in input file
				+"done\n\n")#--reference MN908947.3.fasta
file.close()

j=5000
errorRate = 0.0001
filename = "diffFile_repeat\"$i\"_"+str(j)+"samples_Ns_errors"+ str(errorRate) + ".txt"
pathWrite = "/nfs/research/goldman/" + writeFileUser + "/fastLK/simulations/subsamples/" + str(
	j) + 'subsamples/'  # +"subsamples/"

# Check if files exist
print("for i in $(seq 1 10)\n" + "do \n\t",
      "FILE= ", pathWrite + filename,
      """
	  
		  if [ -f "$FILE" ]; then
			  echo "$FILE exists."
			  FILESIZE=$(stat -c%s "$FILE")
			  MINSIZE = 100
			  if (( FILESIZE < MINSIZE)); then 
				  echo "$FILE too small"
				  head -3 "$FILE"  
			  fi
		  else
			  echo "$FILE does not exist."
		  fi
	  done
	  """
      )

"""RUN MAPLE """

version = "0.1.9_error"
file = open("submitMAPLEv" + version + ".sh", "w")
file.write("rm -r " + console_output + "\n")# + "; mkdir " + console_output + "\n")
file.write("rm -r " + console_errors + "\n")#+"; mkdir " + console_errors + "\n")
file.write(" \n # creating bash scripts for MAPLE \n")
for simulationErrorRate in [0] + errorRates:
	for inferenceErrorRate in [0] + errorRates:
		for j in sampleSizes:
			pathRead = "/nfs/research/goldman/" + readFileUser + "/fastLK/simulations/subsamples/" + str(j) + 'subsamples/'  # +"subsamples/"
			pathWrite = "/nfs/research/goldman/" + writeFileUser + "/fastLK/simulations/subsamples/" + str(j) + 'subsamples/'
			fileName = "diffFile_repeat\"$i\"_"+str(j)+"samples_Ns" + "_errors"+ str(simulationErrorRate) + ".txt" #if errorRate else "diffFile_repeat\"$i\"_"+str(j)+"samples_Ns"  + ".txt"
			treeFile = "treeFile_repeat\"$i\"_" + str(j) +"samples.nw"
			# possible to do remove existing result if already exists - this helps in case of re-running an analysis, and in case the new version fails, to wrongly use the old tree as estimate of the new version.
			# run MAPLE
			file.write("for i in $(seq 1 10)\n do\n\t")
			file.write("bsub -M " + str(int(400 + j / 20))
			            + " -o " + console_output
			            + " -e " + console_errors
			            + " /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 " + folderNameCode +"MAPLEv" + version + ".py "
			            + " --input " + pathWrite + fileName
						+ " --errorRate " + str(inferenceErrorRate)
						+ " --benchmarkingFile " + folderNameCode+"results/benchmarkingFile.tsv"
						+ " --trueTree " + pathRead + treeFile
			            + " --calculateLKfinalTree --overwrite"
			            + " --output " + pathWrite + "outputFile_infError_" + str(inferenceErrorRate) + fileName[8:]
			            + "\n")
			file.write("done\n\n")
file.close()
print("Created MAPLE bash script submitMAPLEv" + version + ".sh")

"""RUN MAPLE ORIGINAL"""
versions = ["0.1.9_original", "0.1.9_error_site_specific"]  #originial
resultFile= "benchmarkingFile2.tsv"# andere benhmark=results folder.2
file = open("submitMAPLEv" + ".sh", "w")
file.write("rm -r " + console_output + "\n")# + "; mkdir " + console_output + "\n")
file.write("rm -r " + console_errors + "\n")#+"; mkdir " + console_errors + "\n")
for version in versions:
	file.write(" \n # creating bash scripts for MAPLE \n")
	for inferenceErrorRate in [0]:
		for simulationErrorRate in [0]:# + errorRates:
			for j in [5000]:##sampleSizes:
				pathRead = "/nfs/research/goldman/" + readFileUser + "/fastLK/simulations/subsamples/" + str(j) + 'subsamples/'  # +"subsamples/"
				pathWrite = "/nfs/research/goldman/" + writeFileUser + "/fastLK/simulations/subsamples/" + str(j) + 'subsamples/'
				fileName = "diffFile_repeat\"$i\"_"+str(j)+"samples_Ns" + "_errors"+ str(simulationErrorRate) + ".txt" #if errorRate else "diffFile_repeat\"$i\"_"+str(j)+"samples_Ns"  + ".txt"
				treeFile = "treeFile_repeat\"$i\"_" + str(j) +"samples.nw"
				# possible to do remove existing result if already exists - this helps in case of re-running an analysis, and in case the new version fails, to wrongly use the old tree as estimate of the new version.
				# run MAPLE
				file.write("for i in $(seq 1 10)\n do\n\t")
				file.write("bsub -M " + str(int(400 + j / 20))
				            + " -o " + console_output
				            + " -e " + console_errors
				            + " /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 " + folderNameCode +"MAPLEv" + version + ".py "
				            + " --input " + pathWrite + fileName
							+ " --errorRate " + str(inferenceErrorRate)
							+ " --benchmarkingFile " + folderNameCode+"results/" + resultFile
							+ " --trueTree " + pathRead + treeFile
				            + " --calculateLKfinalTree --overwrite"
				            + " --output " + pathWrite + "outputFile_infError_" + str(inferenceErrorRate) + fileName[8:]
				            + "\n")
				file.write("done\n\n")
file.close()
print("Created MAPLE bash script submitMAPLEv" + version + ".sh")




# for file siez 5000, rate 0.0001, only 4 results exist.

# MAPLE_options_errortesting = {"--benchmarkingFile ": 'benchmarkingFile.txt' --calculateLKfinalTree True }
#
# --input fastaFile_repeat1_100samples_withRef_Ns_errors.txt
# --overwrite
# --output output/output + input #dependend on input
# --errorRate 0.0005 #dependend on error rate.
# --calculateLKfinalTree True
# --trueTree data/treeFile_repeat1_100samples.nw #dependend on input
# --benchmarkingFile benchmarkingFile.tsv


#######################

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitMAPLEsmall.sh","w")
#, 5000, 10000, 20000, 50000, 100000, 200000, 500000
for j in sampleSizes:
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+
		"rm -f "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_Ns_tree.tree || true\n")

		file.write("bsub -M "+str(int(400+j/20))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_Ns_console_output.txt -e "
		+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_Ns_console_error.txt /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/MAPLEv0.1.1.py"
		+ " --reference /nfs/research/goldman/demaio/fastLK/simulations/MN908947.3.fasta --input "
		+folderNameSimu+str(j)+"subsamples/diffFile_repeat\"$i\"_"+str(j)+"samples_Ns.txt --overwrite --output "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_Ns\n")

		file.write("rm -f "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_4cat_tree.tree || true\n")

		file.write("bsub -M "+str(int(400+j/20))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_4cat_console_output.txt -e "
		+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_4cat_console_error.txt /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/MAPLEv0.1.1.py --reference /nfs/research/goldman/demaio/fastLK/simulations/MN908947.3.fasta --input "
		+folderNameSimu+str(j)+"subsamples/diffFile_4cat_repeat\"$i\"_"+str(j)+"samples.txt --overwrite --output "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_4cat\n"
		+"done\n\n")

file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitMAPLEsmall_evaluations.sh","w")
#[1000, 2000, 5000, 10000, 20000]
for j in [1000,2000]:
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+
		"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*3.0))+" -o "+folderNameSimu+str(j)+
		"subsamples/output_repl\"$i\"/MAPLE004_4cat_evaluation_console_output.txt -e "+folderNameSimu+str(j)+
		"subsamples/output_repl\"$i\"/MAPLE004_4cat_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+
		"subsamples/phylipFile_4cat_repeat\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_4cat_tree.tree -pre "+str(j)+
		"subsamples_repl\"$i\"_MAPLE004Evaluation_4cat -m GTR+G -quiet -nt 1 -keep-ident -redo -blmin 0.000000005 \n")
		file.write("cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*1.7))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_Ns_evaluation_console_output.txt -e "+
		folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_Ns_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+
		"subsamples/phylipFile_repeat\"$i\"_"+str(j)+"samples_Ns.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_Ns_tree.tree -pre "+str(j)+
		"subsamples_repl\"$i\"_MAPLE004Evaluation_Ns -m GTR -quiet -nt 1 -keep-ident -redo -blmin 0.000000005 \n"
		+"done\n\n")
		file.write("\n\n")
file.close()


#write .sh file for new stimations with MAPLE >= v0.04
file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitMAPLEall.sh","w")
#, 5000, 10000, 20000, 50000, 100000, 200000, 500000
for j in [1000,2000]:
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"rm -f "
		+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_GTR_tree.tree || true\n")

		file.write("bsub -M "+str(int(400+j/20))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_GTR_console_output.txt -e "
		+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_GTR_console_error.txt /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/MAPLEv0.1.1.py --reference /nfs/research/goldman/demaio/fastLK/simulations/MN908947.3.fasta --input "
		+folderNameSimu+str(j)+"subsamples/diffFile_repeat\"$i\"_"+str(j)+"samples.txt --overwrite --output "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_GTR\n")

		file.write("rm -f "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_Ns_tree.tree || true\n")

		file.write("bsub -M "+str(int(400+j/20))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_Ns_console_output.txt -e "
		+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_Ns_console_error.txt /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/MAPLEv0.1.1.py --reference /nfs/research/goldman/demaio/fastLK/simulations/MN908947.3.fasta --input "
		+folderNameSimu+str(j)+"subsamples/diffFile_repeat\"$i\"_"+str(j)+"samples_Ns.txt --overwrite --output "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_Ns\n")

		file.write("rm -f "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_4cat_tree.tree || true\n")

		file.write("bsub -M "+str(int(400+j/20))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_4cat_console_output.txt -e "
		+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_4cat_console_error.txt /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/MAPLEv0.1.1.py --reference /nfs/research/goldman/demaio/fastLK/simulations/MN908947.3.fasta --input "
		+folderNameSimu+str(j)+"subsamples/diffFile_4cat_repeat\"$i\"_"+str(j)+"samples.txt --overwrite --output "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_4cat\n") 
		
		file.write("rm -f "+folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_tree.tree || true\n")

		file.write("bsub -M "+str(int(400+j/20))+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_console_output.txt -e "
		+folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_console_error.txt /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/MAPLEv0.1.1.py --reference /nfs/research/goldman/demaio/fastLK/realData/2021-03-31_unmasked_differences_reduced.txt_consensus.fa --input "
		+folderName+str(j)+"subsamples/diffFile_seed\"$i\"_"+str(j)+"samples.txt --overwrite --output "+folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE004\n"+"done\n\n")

file.close()




file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitMAPLEfast_all.sh","w")
#, 5000, 10000, 20000, 50000, 100000, 200000, 500000
for j in [1000,2000]:
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"rm -f "
		+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004fast_GTR_tree.tree || true\n")

		file.write("bsub -M "+str(int(400+j/20))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004fast_GTR_console_output.txt -e "
		+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004fast_GTR_console_error.txt /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/MAPLEv0.0.6.py --reference /nfs/research/goldman/demaio/fastLK/simulations/MN908947.3.fasta --input "
		+folderNameSimu+str(j)+"subsamples/diffFile_repeat\"$i\"_"+str(j)+"samples.txt --binaryTree --fast --overwrite --output "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004fast_GTR\n")

		file.write("rm -f "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004fast_Ns_tree.tree || true\n")

		file.write("bsub -M "+str(int(400+j/20))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004fast_Ns_console_output.txt -e "
		+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004fast_Ns_console_error.txt /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/MAPLEv0.0.6.py --reference /nfs/research/goldman/demaio/fastLK/simulations/MN908947.3.fasta --input "
		+folderNameSimu+str(j)+"subsamples/diffFile_repeat\"$i\"_"+str(j)+"samples_Ns.txt --binaryTree --fast --overwrite --output "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004fast_Ns\n")

		file.write("rm -f "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004fast_4cat_tree.tree || true\n")

		file.write("bsub -M "+str(int(400+j/20))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004fast_4cat_console_output.txt -e "
		+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004fast_4cat_console_error.txt /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/MAPLEv0.0.6.py --reference /nfs/research/goldman/demaio/fastLK/simulations/MN908947.3.fasta --input "
		+folderNameSimu+str(j)+"subsamples/diffFile_4cat_repeat\"$i\"_"+str(j)+"samples.txt --binaryTree --fast --overwrite --output "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004fast_4cat\n") 
		
		file.write("rm -f "+folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE004fast_tree.tree || true\n")

		file.write("bsub -M "+str(int(400+j/20))+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE004fast_console_output.txt -e "
		+folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE004fast_console_error.txt /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/MAPLEv0.0.6.py --reference /nfs/research/goldman/demaio/fastLK/realData/2021-03-31_unmasked_differences_reduced.txt_consensus.fa --input "
		+folderName+str(j)+"subsamples/diffFile_seed\"$i\"_"+str(j)+"samples.txt --binaryTree --fast --overwrite --output "+folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE004fast\n"+"done\n\n")

file.close()

#run new MAPLE evaluations
file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitMAPLEall_evaluations.sh","w")
#[1000, 2000, 5000, 10000, 20000]
for j in [1000,2000]:
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+
		"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*1.7))+" -o "+folderNameSimu+str(j)+
		"subsamples/output_repl\"$i\"/MAPLE004_GTR_evaluation_console_output.txt -e "+folderNameSimu+str(j)+
		"subsamples/output_repl\"$i\"/MAPLE004_GTR_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_repeat\"$i\"_"+str(j)+
		"samples.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_GTR_tree.tree -pre "+str(j)+"subsamples_repl\"$i\"_MAPLE004Evaluation_GTR -m GTR -quiet -nt 1 -keep-ident -redo -blmin 0.000000005 \n")
		file.write("cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*3.0))+" -o "+folderNameSimu+str(j)+
		"subsamples/output_repl\"$i\"/MAPLE004_4cat_evaluation_console_output.txt -e "+folderNameSimu+str(j)+
		"subsamples/output_repl\"$i\"/MAPLE004_4cat_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+
		"subsamples/phylipFile_4cat_repeat\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_4cat_tree.tree -pre "+str(j)+
		"subsamples_repl\"$i\"_MAPLE004Evaluation_4cat -m GTR+G -quiet -nt 1 -keep-ident -redo -blmin 0.000000005 \n")
		file.write("cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*1.7))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_Ns_evaluation_console_output.txt -e "+
		folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_Ns_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+
		"subsamples/phylipFile_repeat\"$i\"_"+str(j)+"samples_Ns.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_Ns_tree.tree -pre "+str(j)+
		"subsamples_repl\"$i\"_MAPLE004Evaluation_Ns -m GTR -quiet -nt 1 -keep-ident -redo -blmin 0.000000005 \n")
		file.write("cd "+folderName+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(j)+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_evaluation_console_output.txt -e "+
		folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderName+str(j)+
		"subsamples/phylipFile_seed\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE004_tree.tree -pre "+str(j)+
		"subsamples_repl\"$i\"_MAPLE004Evaluation -m GTR -quiet -nt 1 -keep-ident -redo -blmin 0.000000005 \n"+"done\n\n")
		file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitMAPLEfast_all_evaluations.sh","w")
#[1000, 2000, 5000, 10000, 20000]
for j in [1000,2000]:
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+
		"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*1.7))+" -o "+folderNameSimu+str(j)+
		"subsamples/output_repl\"$i\"/MAPLE004fast_GTR_evaluation_console_output.txt -e "+folderNameSimu+str(j)+
		"subsamples/output_repl\"$i\"/MAPLE004fast_GTR_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_repeat\"$i\"_"+str(j)+
		"samples.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004fast_GTR_tree.tree -pre "+str(j)+"subsamples_repl\"$i\"_MAPLE004fastEvaluation_GTR -m GTR -quiet -nt 1 -keep-ident -redo -blmin 0.000000005 \n")
		file.write("cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*3.0))+" -o "+folderNameSimu+str(j)+
		"subsamples/output_repl\"$i\"/MAPLE004fast_4cat_evaluation_console_output.txt -e "+folderNameSimu+str(j)+
		"subsamples/output_repl\"$i\"/MAPLE004fast_4cat_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+
		"subsamples/phylipFile_4cat_repeat\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004fast_4cat_tree.tree -pre "+str(j)+
		"subsamples_repl\"$i\"_MAPLE004fastEvaluation_4cat -m GTR+G -quiet -nt 1 -keep-ident -redo -blmin 0.000000005 \n")
		file.write("cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*1.7))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004fast_Ns_evaluation_console_output.txt -e "+
		folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004fast_Ns_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+
		"subsamples/phylipFile_repeat\"$i\"_"+str(j)+"samples_Ns.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE004fast_Ns_tree.tree -pre "+str(j)+
		"subsamples_repl\"$i\"_MAPLE004fastEvaluation_Ns -m GTR -quiet -nt 1 -keep-ident -redo -blmin 0.000000005 \n")
		file.write("cd "+folderName+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(j)+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE004fast_evaluation_console_output.txt -e "+
		folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE004fast_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderName+str(j)+
		"subsamples/phylipFile_seed\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE004fast_tree.tree -pre "+str(j)+
		"subsamples_repl\"$i\"_MAPLE004fastEvaluation -m GTR -quiet -nt 1 -keep-ident -redo -blmin 0.000000005 \n"+"done\n\n")
		file.write("\n\n")
file.close()

#new UShER runs
speedsUsher=["medium"]
file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitUsherNew.sh","w")
notOnlyNs=False
#1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000
for j in [1000,2000]:
	for s in range(len(speedsUsher)):
			
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/200))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usherNew_"+speedsUsher[s]+"_console_output.txt -e "+
		folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usherNew_"+speedsUsher[s]+"_console_error.txt usher -v "+folderNameSimu+str(j)+"subsamples/vcfFile_repeat\"$i\"_"+str(j)+
		"samples.vcf  -d "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_output/ -T 1 -t "+folderNameSimu+str(j)+"subsamples/initialTree_repeat\"$i\"_"+str(j)+"samples.tree \n"+"done\n\n")

		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/200))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usherNew_4cat_"+speedsUsher[s]+"_console_output.txt -e "+
		folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usherNew_4cat_"+speedsUsher[s]+"_console_error.txt usher -v "+folderNameSimu+str(j)+"subsamples/vcfFile_4cat_repeat\"$i\"_"+str(j)+
		"samples.vcf  -d "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_output_4cat/ -T 1 -t "+folderNameSimu+str(j)+"subsamples/initialTree_4cat_repeat\"$i\"_"+str(j)+"samples.tree \n"+"done\n\n")
	
		if notOnlyNs:

			file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/50))+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/usherNew_"+speedsUsher[s]+"_console_output.txt -e "+
			folderName+str(j)+"subsamples/output_repl\"$i\"/usherNew_"+speedsUsher[s]+"_console_error.txt usher -v "+folderName+str(j)+"subsamples/vcfFile_seed\"$i\"_"+str(j)+
			"samples.vcf  -d "+folderName+str(j)+"subsamples/output_repl\"$i\"/usher_output/ -T 1 -t "+folderName+str(j)+"subsamples/initialTree_seed\"$i\"_"+str(j)+"samples.tree  \n"+"done\n\n")

			file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/50))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usherNew_"+speedsUsher[s]+"_Ns_console_output.txt -e "+
			folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usherNew_"+speedsUsher[s]+"_Ns_console_error.txt usher -v "+folderNameSimu+str(j)+"subsamples/vcfFile_repeat\"$i\"_"+str(j)+
			"samples_Ns.vcf  -d "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_output_Ns/ -T 1 -t "+folderNameSimu+str(j)+"subsamples/initialTree_repeat\"$i\"_"+str(j)+"samples_Ns.tree \n"+"done\n\n")
		
	file.write("\n\n")
file.close()

speedsMAT=["medium"]
file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitMatOptimizeNew.sh","w")
for j in [1000,2000]:
	for s in range(len(speedsMAT)):
		if notOnlyNs:
			file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/5))+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_console_output.txt -e "+
			folderName+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_console_error.txt matOptimize -n -v "+folderName+str(j)+"subsamples/vcfFile_seed\"$i\"_"+str(j)+"samples.vcf  -o "+
			folderName+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_output.pb -T 1 -t "+folderName+str(j)+"subsamples/output_repl\"$i\"/usher_output/final-tree.nh  \n"+"done\n\n")
			
			file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/5))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_Ns_console_output.txt -e "+
			folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_Ns_console_error.txt matOptimize -n -v "+folderNameSimu+str(j)+"subsamples/vcfFile_repeat\"$i\"_"+str(j)+"samples_Ns.vcf  -o "+
			folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_output_Ns.pb -T 1 -t "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_output_Ns/final-tree.nh \n"+"done\n\n")
		
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/5))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_console_output.txt -e "+
		folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_console_error.txt matOptimize -n -v "+folderNameSimu+str(j)+"subsamples/vcfFile_repeat\"$i\"_"+str(j)+"samples.vcf  -o "+
		folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_output.pb -T 1 -t "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_output/final-tree.nh \n"+"done\n\n")

		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/5))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_4cat_"+speedsMAT[s]+"_console_output.txt -e "+
		folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_4cat_"+speedsMAT[s]+"_console_error.txt matOptimize -n -v "+folderNameSimu+str(j)+"subsamples/vcfFile_4cat_repeat\"$i\"_"+str(j)+"samples.vcf  -o "+
		folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_output_4cat.pb -T 1 -t "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_output_4cat/final-tree.nh \n"+"done\n\n")
	
		
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitMatOptimize_convertOutputTreeNew.sh","w")
for j in [1000,2000]:
	for s in range(len(speedsMAT)):
		if notOnlyNs:
			file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/5))+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/matOptimizeConvert_"+speedsMAT[s]+"_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/matOptimizeConvert_"+speedsMAT[s]+"_console_error.txt matUtils extract -i "+folderName+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_output.pb -T 1 -d "+folderName+str(j)+"subsamples/output_repl\"$i\"/ -t matOptimize_"+speedsMAT[s]+"_output.nh  \n"+"done\n\n")

			file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/5))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimizeConvert_"+speedsMAT[s]+"_Ns_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimizeConvert_"+speedsMAT[s]+"_Ns_console_error.txt matUtils extract -i "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_output_Ns.pb -T 1 -d "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/ -t matOptimize_"+speedsMAT[s]+"_output_Ns.nh \n"+"done\n\n")
	
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/5))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimizeConvert_"+speedsMAT[s]+"_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimizeConvert_"+speedsMAT[s]+"_console_error.txt matUtils extract -i "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_output.pb -T 1 -d "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/ -t matOptimize_"+speedsMAT[s]+"_output.nh \n"+"done\n\n")

		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/5))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimizeConvert_4cat_"+speedsMAT[s]+"_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimizeConvert_4cat_"+speedsMAT[s]+"_console_error.txt matUtils extract -i "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_output_4cat.pb -T 1 -d "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/ -t matOptimize_"+speedsMAT[s]+"_output_4cat.nh \n"+"done\n\n")
	
		
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitUsherEvaluationNew.sh","w")
for j in [1000, 2000]:
	for s in range(len(speedsUsher)):
		if notOnlyNs:
			file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"cd "+folderName+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(j)+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/usher_"+speedsUsher[s]+"_evaluation_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/usher_"+speedsUsher[s]+"_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderName+str(j)+"subsamples/phylipFile_seed\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderName+str(j)+"subsamples/output_repl\"$i\"/usher_output/final-tree.nh -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsUsher[s]+"_usherEvaluation -m GTR -quiet -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
			
			file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*1.7))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_"+speedsUsher[s]+"_Ns_evaluation_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_"+speedsUsher[s]+"_Ns_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_repeat\"$i\"_"+str(j)+"samples_Ns.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_output_Ns/final-tree.nh -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsUsher[s]+"_Ns_usherEvaluation -m GTR -quiet -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")

		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*3.0))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_4cat_"+speedsUsher[s]+"_evaluation_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_4cat_"+speedsUsher[s]+"_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_4cat_repeat\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_output_4cat/final-tree.nh -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsUsher[s]+"_usherEvaluation_4cat -m GTR+G -redo -quiet -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
		
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*1.7))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_"+speedsUsher[s]+"_evaluation_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_"+speedsUsher[s]+"_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_repeat\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_output/final-tree.nh -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsUsher[s]+"_usherEvaluation -m GTR -quiet -nt 1 -keep-ident -redo -blmin 0.000000005 \n"+"done\n\n")
		
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitMatOptimizeEvaluationNew.sh","w")
for j in [1000, 2000]:
	for s in range(len(speedsMAT)):
		if notOnlyNs:
			file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"cd "+folderName+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*1.7))+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_evaluation_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderName+str(j)+"subsamples/phylipFile_seed\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderName+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_output.nh -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsMAT[s]+"_matOptimizeEvaluation -m GTR -quiet -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")

			file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*1.7))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_Ns_evaluation_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_Ns_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_repeat\"$i\"_"+str(j)+"samples_Ns.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_output_Ns.nh -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsMAT[s]+"_Ns_matOptimizeEvaluation -m GTR -quiet -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")

		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*3.0))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_4cat_"+speedsMAT[s]+"_evaluation_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_4cat_"+speedsMAT[s]+"_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_4cat_repeat\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_output_4cat.nh -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsMAT[s]+"_matOptimizeEvaluation_4cat -m GTR+G -quiet -redo -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
		
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*1.7))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_evaluation_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_repeat\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_output.nh -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsMAT[s]+"_matOptimizeEvaluation -m GTR -quiet -nt 1 -keep-ident -redo -blmin 0.000000005 \n"+"done\n\n")
		
	file.write("\n\n")
file.close()





# file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitFastLKnew.sh","w")
# #for j in [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000]:
# for j in [20000, 100000]:
# 	for s in range(len(speeds)):
# 		file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"bsub -M "+str(int(200+j/5))+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastLK_"+speeds[s]+"_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastLK_"+speeds[s]+"_console_error.txt /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/estimatePhylogenyIterativeFastLK.py --reference /nfs/research/goldman/demaio/fastLK/realData/2021-03-31_unmasked_differences_reduced.txt_consensus.fa --input "+folderName+str(j)+"subsamples/diffFile_seed\"$i\"_"+str(j)+"samples.txt --allowedFails "+str(allowedFails[s])+" --thresholdLogLK "+str(thresholdLogLKs[s])+" --bLenAdjustment "+str(bLenAdjustments[s])+" --bLenFactor "+str(bLenFactors[s])+"  --output "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastLK_"+speeds[s]+"\n"+"done\n\n")
# 	file.write("\n\n")
# file.close()

# file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitFastLKevaluation.sh","w")
# for j in [100, 200, 500, 1000, 2000, 5000, 10000, 20000]:
# 	for s in range(len(speeds)):
# 		file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"cd "+folderName+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(j)+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastLK_"+speeds[s]+"_evaluation_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastLK_"+speeds[s]+"_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderName+str(j)+"subsamples/phylipFile_seed\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastLK_"+speeds[s]+"_tree.tree -pre "+str(j)+"subsamples_repl\"$i\"_"+speeds[s]+"_fastLKevaluation -m GTR -quiet -redo -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
# 	file.write("\n\n")
# file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitFastLKtopology.sh","w")
for j in [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000]:
	for s in range(len(speeds)):
		file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/3))+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"_console_error.txt /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/estimatePhylogenyIterativeFastLK.py --reference /nfs/research/goldman/demaio/fastLK/realData/2021-03-31_unmasked_differences_reduced.txt_consensus.fa --input "+folderName+str(j)+"subsamples/diffFile_seed\"$i\"_"+str(j)+"samples.txt --binaryTree --allowedFails "+str(allowedFails[s])+" --thresholdLogLK "+str(thresholdLogLKs[s])+" --bLenAdjustment "+str(bLenAdjustments[s])+" --bLenFactor "+str(bLenFactors[s])+" --numTopologyImprovements "+str(numTopologyImprovements[s])+" --thresholdTopologyPlacement "+str(thresholdTopologyPlacement[s])+" --allowedFailsTopology "+str(allowedFailsTopology[s])+" --thresholdLogLKtopology "+str(thresholdLogLKtopology[s])+" --output "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"\n"+"done\n\n")
		#file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"bsub -M "+str(int(200+j/8))+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"_console_error.txt /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/estimatePhylogenyIterativeFastLK.py --reference /nfs/research/goldman/demaio/fastLK/realData/2021-03-31_unmasked_differences_reduced.txt_consensus.fa --input "+folderName+str(j)+"subsamples/diffFile_seed\"$i\"_"+str(j)+"samples.txt --binaryTree --overwrite --allowedFails "+str(allowedFails[s])+" --thresholdLogLK "+str(thresholdLogLKs[s])+" --bLenAdjustment "+str(bLenAdjustments[s])+" --bLenFactor "+str(bLenFactors[s])+" --numTopologyImprovements "+str(numTopologyImprovements[s])+" --thresholdTopologyPlacement "+str(thresholdTopologyPlacement[s])+" --allowedFailsTopology "+str(allowedFailsTopology[s])+" --thresholdLogLKtopology "+str(thresholdLogLKtopology[s])+" --output "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"\n"+"done\n\n")
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitFastLKtopologyEvaluation.sh","w")
for j in [100, 200, 500, 1000, 2000, 5000, 10000, 20000]:
	for s in range(len(speeds)):
		file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"cd "+folderName+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(j)+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"_evaluation_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderName+str(j)+"subsamples/phylipFile_seed\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"_tree.tree -pre "+str(j)+"subsamples_repl\"$i\"_"+speeds[s]+"_fastLKtopologyEvaluation -m GTR -quiet -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
		#file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"cd "+folderName+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(j)+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"_evaluation_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderName+str(j)+"subsamples/phylipFile_seed\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"_tree.tree -pre "+str(j)+"subsamples_repl\"$i\"_"+speeds[s]+"_fastLKtopologyEvaluation -m GTR -quiet -redo -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
	file.write("\n\n")
file.close()


#run fastLK on simulated data
file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitFastLKtopology_simulations.sh","w")
for j in [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000]:
	for s in range(len(speeds)):
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/3))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"_console_error.txt /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/estimatePhylogenyIterativeFastLK.py --reference  /nfs/research/goldman/demaio/fastLK/simulations/MN908947.3.fasta --input "+folderNameSimu+str(j)+"subsamples/diffFile_repeat\"$i\"_"+str(j)+"samples.txt --binaryTree --allowedFails "+str(allowedFails[s])+" --thresholdLogLK "+str(thresholdLogLKs[s])+" --model UNREST --bLenAdjustment "+str(bLenAdjustments[s])+" --bLenFactor "+str(bLenFactors[s])+" --numTopologyImprovements "+str(numTopologyImprovements[s])+" --thresholdTopologyPlacement "+str(thresholdTopologyPlacement[s])+" --allowedFailsTopology "+str(allowedFailsTopology[s])+" --thresholdLogLKtopology "+str(thresholdLogLKtopology[s])+" --overwrite  --output "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"\n"+"done\n\n")

		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/3))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_4cat_"+speeds[s]+"_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_4cat_"+speeds[s]+"_console_error.txt /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/estimatePhylogenyIterativeFastLK.py --reference  /nfs/research/goldman/demaio/fastLK/simulations/MN908947.3.fasta --input "+folderNameSimu+str(j)+"subsamples/diffFile_4cat_repeat\"$i\"_"+str(j)+"samples.txt --binaryTree --allowedFails "+str(allowedFails[s])+" --thresholdLogLK "+str(thresholdLogLKs[s])+" --model UNREST --bLenAdjustment "+str(bLenAdjustments[s])+" --bLenFactor "+str(bLenFactors[s])+" --numTopologyImprovements "+str(numTopologyImprovements[s])+" --thresholdTopologyPlacement "+str(thresholdTopologyPlacement[s])+" --allowedFailsTopology "+str(allowedFailsTopology[s])+" --thresholdLogLKtopology "+str(thresholdLogLKtopology[s])+" --overwrite --output "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_4cat_"+speeds[s]+"\n"+"done\n\n")
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitFastLKtopology_simulations_GTR.sh","w")
for j in [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000]:
	for s in range(len(speeds)):
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/3))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_GTR_"+speeds[s]+"_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_GTR_"+speeds[s]+"_console_error.txt /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/estimatePhylogenyIterativeFastLK.py --reference  /nfs/research/goldman/demaio/fastLK/simulations/MN908947.3.fasta --input "+folderNameSimu+str(j)+"subsamples/diffFile_repeat\"$i\"_"+str(j)+"samples.txt --binaryTree --allowedFails "+str(allowedFails[s])+" --thresholdLogLK "+str(thresholdLogLKs[s])+" --bLenAdjustment "+str(bLenAdjustments[s])+" --bLenFactor "+str(bLenFactors[s])+" --numTopologyImprovements "+str(numTopologyImprovements[s])+" --thresholdTopologyPlacement "+str(thresholdTopologyPlacement[s])+" --allowedFailsTopology "+str(allowedFailsTopology[s])+" --thresholdLogLKtopology "+str(thresholdLogLKtopology[s])+" --overwrite --output "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_GTR_"+speeds[s]+"\n"+"done\n\n")
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitFastLKtopology_simulations_Ns.sh","w")
for j in [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000]:
	for s in range(len(speeds)):
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/3))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"_Ns_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"_Ns_console_error.txt /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/estimatePhylogenyIterativeFastLK.py --reference  /nfs/research/goldman/demaio/fastLK/simulations/MN908947.3.fasta --input "+folderNameSimu+str(j)+"subsamples/diffFile_repeat\"$i\"_"+str(j)+"samples_Ns.txt --binaryTree --allowedFails "+str(allowedFails[s])+" --thresholdLogLK "+str(thresholdLogLKs[s])+" --bLenAdjustment "+str(bLenAdjustments[s])+" --bLenFactor "+str(bLenFactors[s])+" --numTopologyImprovements "+str(numTopologyImprovements[s])+" --thresholdTopologyPlacement "+str(thresholdTopologyPlacement[s])+" --allowedFailsTopology "+str(allowedFailsTopology[s])+" --thresholdLogLKtopology "+str(thresholdLogLKtopology[s])+" --overwrite  --output "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"_Ns\n"+"done\n\n")
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitFastLKtopology_simulationsEvaluation.sh","w")
for j in [100, 200, 500, 1000, 2000, 5000, 10000, 20000]:
	for s in range(len(speeds)):
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*1.7))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_GTR_"+speeds[s]+"_evaluation_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_GTR_"+speeds[s]+"_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_repeat\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_GTR_"+speeds[s]+"_tree.tree -pre "+str(j)+"subsamples_repl\"$i\"_"+speeds[s]+"_fastLKtopologyEvaluation_GTR -m GTR -quiet -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*3.0))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_4cat_"+speeds[s]+"_evaluation_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_4cat_"+speeds[s]+"_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_4cat_repeat\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_4cat_"+speeds[s]+"_tree.tree -pre "+str(j)+"subsamples_repl\"$i\"_"+speeds[s]+"_fastLKtopologyEvaluation_4cat -m GTR+G -quiet -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*1.7))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"_evaluation_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_repeat\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"_tree.tree -pre "+str(j)+"subsamples_repl\"$i\"_"+speeds[s]+"_fastLKtopologyEvaluation -m GTR -quiet -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
		#file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"cd "+folderName+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(j)+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"_evaluation_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderName+str(j)+"subsamples/phylipFile_seed\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"_tree.tree -pre "+str(j)+"subsamples_repl\"$i\"_"+speeds[s]+"_fastLKtopologyEvaluation -m GTR -quiet -redo -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitFastLKtopology_simulationsEvaluation_Ns.sh","w")
for j in [100, 200, 500, 1000, 2000, 5000, 10000, 20000]:
	for s in range(len(speeds)):
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*1.7))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"_Ns_evaluation_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"_Ns_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_repeat\"$i\"_"+str(j)+"samples_Ns.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastLKtopology_"+speeds[s]+"_Ns_tree.tree -pre "+str(j)+"subsamples_repl\"$i\"_"+speeds[s]+"_Ns_fastLKtopologyEvaluation -m GTR -quiet -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
	file.write("\n\n")
file.close()




#test new MAPLE version
speedsMAPLE=["4-160-4-120","5-160-4-120","4-180-4-120","4-200-4-120","nonStrictInitialStopRules","4-160-5-120","4-160-4-100","4-160-4-140","4-160-4-160","StrictTopologyStopRules","double"]
allowedFailsMAPLE=[4,5,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4]
thresholdLogLKsMAPLE=[160.0,160.0,180.0,200.0,160.0,160.0,160.0,160.0,160.0,160.0,160.0,160.0,160.0,160.0,160.0,160.0,160.0,160.0,160.0,160.0,160.0,160.0,160.0,160.0,160.0,160.0,160.0,160.0]
strictRule=["","","",""," --nonStrictInitialStopRules","","","","","","","","","","","","","","","","","","","","","","","",""]
#strictRule=[" --nonStrictInitialStopRules"," --nonStrictInitialStopRules"," --nonStrictInitialStopRules"," --nonStrictInitialStopRules"," --nonStrictInitialStopRules"," --nonStrictInitialStopRules"," --nonStrictInitialStopRules"," --nonStrictInitialStopRules"," --nonStrictInitialStopRules"," --nonStrictInitialStopRules"," --nonStrictInitialStopRules"," --nonStrictInitialStopRules"," --nonStrictInitialStopRules"," --nonStrictInitialStopRules"]
bLenFactorsMAPLE=[4.1,4.1,4.1,4.1,4.1,4.1,4.1,4.1,4.1,4.1,4.1,4.1,4.1,4.1,4.1,4.1,4.1,4.1,4.1,4.1,4.1,4.1,4.1,4.1,4.1,4.1,4.1,4.1]
thresholdProbs=[0.00000001,0.00000001,0.00000001,0.00000001,0.00000001,0.00000001,0.00000001,0.00000001,0.00000001,0.00000001,0.00000001,0.00000001,0.00000001,0.00000001,0.00000001,0.00000001]
updateSubstMatrixEveryThisSamples=[25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25]
thresholdDiffForUpdate=[0.0000001,0.0000001,0.0000001,0.0000001,0.0000001,0.0000001,0.0000001,0.0000001,0.0000001,0.0000001,0.0000001,0.0000001,0.0000001,0.0000001,0.0000001,0.0000001]
thresholdFoldChangeUpdate=[1.001,1.001,1.001,1.001,1.001,1.001,1.001,1.001,1.001,1.001,1.001,1.001,1.001,1.001,1.001,1.001]
minBLen=[0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2]
maxBLen=[40.0,40.0,40.0,40.0,40.0,40.0,40.0,40.0,40.0,40.0,40.0,40.0,40.0,40.0,40.0,40.0]
thresholdLogLKconsecutivePlacement=[0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01]
bLenAdjustmentsMAPLE=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
numTopologyImprovementsMAPLE=[1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1]
allowedFailsTopologyMAPLE=[4,4,4,4,4,5,4,4,4,4,4,4,4,4,4,4]
thresholdLogLKtopologyMAPLE=[120.0,120.0,120.0,120.0,120.0,120.0,100.0,140.0,160.0,120.0,120.0,120.0,120.0,120.0]
thresholdTopologyPlacementMAPLE=[-0.00001,-0.00001,-0.00001,-0.00001,-0.00001,-0.00001,-0.00001,-0.00001,-0.00001,-0.00001,-0.00001,-0.00001,-0.00001,-0.00001,-0.00001,-0.00001]
#nonStrictTopologyStopRules=["","","","","",""," --nonStrictTopologyStopRules"," --nonStrictTopologyStopRules"," --nonStrictTopologyStopRules"," --nonStrictTopologyStopRules"," --nonStrictTopologyStopRules"," --nonStrictTopologyStopRules","",""]
nonStrictTopologyStopRules=[" --nonStrictTopologyStopRules"," --nonStrictTopologyStopRules"," --nonStrictTopologyStopRules"," --nonStrictTopologyStopRules"," --nonStrictTopologyStopRules"," --nonStrictTopologyStopRules"," --nonStrictTopologyStopRules"," --nonStrictTopologyStopRules"," --nonStrictTopologyStopRules",""," --nonStrictTopologyStopRules","","","","","","","","","","","","","",""]

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/MAPLE_testingParameters.sh","w")
file2=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/MAPLE_testingParameters_evaluation.sh","w")
#for j in [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000]:
for j in [20000,50000,100000,200000]:
	for s in range(len(speedsMAPLE)):
		
		testSimu=False
		if testSimu==True:
			file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"rm -f "
			+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE_GTR_"+speedsMAPLE[s]+"_tree.tree || true\n")
			#+"done\n\n")

			#file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/6))+" -o "
			file.write("bsub -M "+str(int(400+j/15))+" -o "
			+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE_GTR_"+speedsMAPLE[s]+"_console_output.txt -e "
			+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE_GTR_"+speedsMAPLE[s]
			+"_console_error.txt /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/MAPLEv0.0.4.py --reference /nfs/research/goldman/demaio/fastLK/simulations/MN908947.3.fasta --input "
			+folderNameSimu+str(j)+"subsamples/diffFile_repeat\"$i\"_"+str(j)+"samples.txt --binaryTree --allowedFails "+str(allowedFailsMAPLE[s])+" --thresholdLogLK "
			+str(thresholdLogLKsMAPLE[s])+" --bLenAdjustment "+str(bLenAdjustmentsMAPLE[s])+" --minBLenForMidNode "+str(bLenFactorsMAPLE[s])+" --numTopologyImprovements "
			+str(numTopologyImprovementsMAPLE[s])+" --thresholdTopologyPlacement "+f'{thresholdTopologyPlacementMAPLE[s]:.10f}'+" --allowedFailsTopology "
			+str(allowedFailsTopologyMAPLE[s])+" --thresholdLogLKtopology "+str(thresholdLogLKtopologyMAPLE[s])+strictRule[s]+" --thresholdProb "+str(thresholdProbs[s])+" --updateSubstMatrixEveryThisSamples "+str(updateSubstMatrixEveryThisSamples[s])
			+" --thresholdDiffForUpdate "+str(thresholdDiffForUpdate[s])+" --thresholdFoldChangeUpdate "+str(thresholdFoldChangeUpdate[s])+" --minBLen "+str(minBLen[s])+" --maxBLen "+str(maxBLen[s])+nonStrictTopologyStopRules[s]+" --overwrite --thresholdLogLKconsecutivePlacement "+str(thresholdLogLKconsecutivePlacement[s])+" --output "
			+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE_GTR_"+speedsMAPLE[s]+"\n") # 
			#+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE_GTR_"+speedsMAPLE[s]+"\n"+"done\n\n") # 

			#file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"rm -f "
			#+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE_Ns_"+speedsMAPLE[s]+"_tree.tree || true\n"+"done\n\n")
			file.write("rm -f "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE_Ns_"+speedsMAPLE[s]+"_tree.tree || true\n")

			#file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/6))+" -o "
			file.write("bsub -M "+str(int(400+j/15))+" -o "
			+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE_Ns_"+speedsMAPLE[s]+"_console_output.txt -e "
			+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE_Ns_"+speedsMAPLE[s]
			+"_console_error.txt /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/MAPLEv0.0.4.py --reference /nfs/research/goldman/demaio/fastLK/simulations/MN908947.3.fasta --input "
			+folderNameSimu+str(j)+"subsamples/diffFile_repeat\"$i\"_"+str(j)+"samples_Ns.txt --binaryTree --allowedFails "+str(allowedFailsMAPLE[s])+" --thresholdLogLK "
			+str(thresholdLogLKsMAPLE[s])+" --bLenAdjustment "+str(bLenAdjustmentsMAPLE[s])+" --minBLenForMidNode "+str(bLenFactorsMAPLE[s])+" --numTopologyImprovements "
			+str(numTopologyImprovementsMAPLE[s])+" --thresholdTopologyPlacement "+f'{thresholdTopologyPlacementMAPLE[s]:.10f}'+" --allowedFailsTopology "
			+str(allowedFailsTopologyMAPLE[s])+" --thresholdLogLKtopology "+str(thresholdLogLKtopologyMAPLE[s])+strictRule[s]+" --thresholdProb "+str(thresholdProbs[s])+" --updateSubstMatrixEveryThisSamples "+str(updateSubstMatrixEveryThisSamples[s])
			+" --thresholdDiffForUpdate "+str(thresholdDiffForUpdate[s])+" --thresholdFoldChangeUpdate "+str(thresholdFoldChangeUpdate[s])+" --minBLen "+str(minBLen[s])+" --maxBLen "+str(maxBLen[s])+nonStrictTopologyStopRules[s]+" --overwrite --thresholdLogLKconsecutivePlacement "+str(thresholdLogLKconsecutivePlacement[s])+" --output "
			+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE_Ns_"+speedsMAPLE[s]+"\n") # 
			#+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE_Ns_"+speedsMAPLE[s]+"\n"+"done\n\n") # 

			#file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"rm -f "
			#+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE_4cat_"+speedsMAPLE[s]+"_tree.tree || true\n"+"done\n\n")
			file.write("rm -f "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE_4cat_"+speedsMAPLE[s]+"_tree.tree || true\n")

			#file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/6))+" -o "
			file.write("bsub -M "+str(int(400+j/15))+" -o "
			+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE_4cat_"+speedsMAPLE[s]+"_console_output.txt -e "
			+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE_4cat_"+speedsMAPLE[s]
			+"_console_error.txt /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/MAPLEv0.0.4.py --reference /nfs/research/goldman/demaio/fastLK/simulations/MN908947.3.fasta --input "
			+folderNameSimu+str(j)+"subsamples/diffFile_4cat_repeat\"$i\"_"+str(j)+"samples.txt --binaryTree --allowedFails "+str(allowedFailsMAPLE[s])+" --thresholdLogLK "
			+str(thresholdLogLKsMAPLE[s])+" --bLenAdjustment "+str(bLenAdjustmentsMAPLE[s])+" --minBLenForMidNode "+str(bLenFactorsMAPLE[s])+" --numTopologyImprovements "
			+str(numTopologyImprovementsMAPLE[s])+" --thresholdTopologyPlacement "+f'{thresholdTopologyPlacementMAPLE[s]:.10f}'+" --allowedFailsTopology "
			+str(allowedFailsTopologyMAPLE[s])+" --thresholdLogLKtopology "+str(thresholdLogLKtopologyMAPLE[s])+strictRule[s]+" --thresholdProb "+str(thresholdProbs[s])+" --updateSubstMatrixEveryThisSamples "+str(updateSubstMatrixEveryThisSamples[s])
			+" --thresholdDiffForUpdate "+str(thresholdDiffForUpdate[s])+" --thresholdFoldChangeUpdate "+str(thresholdFoldChangeUpdate[s])+" --minBLen "+str(minBLen[s])+" --maxBLen "+str(maxBLen[s])+nonStrictTopologyStopRules[s]+" --overwrite --thresholdLogLKconsecutivePlacement "+str(thresholdLogLKconsecutivePlacement[s])+" --output "
			+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE_4cat_"+speedsMAPLE[s]+"\n") #
			#+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE_4cat_"+speedsMAPLE[s]+"\n"+"done\n\n") #

			#file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"rm -f "
			#+folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE_"+speedsMAPLE[s]+"_tree.tree || true\n"+"done\n\n")
			file.write("rm -f "+folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE_"+speedsMAPLE[s]+"_tree.tree || true\n")

		else:
			file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"rm -f "
			+folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE_"+speedsMAPLE[s]+"_tree.tree || true\n")

		#if j<20000:
		#file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/6))+" -o "
		file.write("bsub -M "+str(int(400+j/15))+" -o "
		+folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE_"+speedsMAPLE[s]+"_console_output.txt -e "
		+folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE_"+speedsMAPLE[s]
		+"_console_error.txt /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/MAPLEv0.0.4.py --reference /nfs/research/goldman/demaio/fastLK/realData/2021-03-31_unmasked_differences_reduced.txt_consensus.fa --input "
		+folderName+str(j)+"subsamples/diffFile_seed\"$i\"_"+str(j)+"samples.txt --binaryTree --allowedFails "+str(allowedFailsMAPLE[s])+" --thresholdLogLK "
		+str(thresholdLogLKsMAPLE[s])+" --bLenAdjustment "+str(bLenAdjustmentsMAPLE[s])+" --minBLenForMidNode "+str(bLenFactorsMAPLE[s])+" --numTopologyImprovements "
		+str(numTopologyImprovementsMAPLE[s])+" --thresholdTopologyPlacement "+f'{thresholdTopologyPlacementMAPLE[s]:.10f}'+" --allowedFailsTopology "
		+str(allowedFailsTopologyMAPLE[s])+" --thresholdLogLKtopology "+str(thresholdLogLKtopologyMAPLE[s])+strictRule[s]+" --thresholdProb "+str(thresholdProbs[s])+" --updateSubstMatrixEveryThisSamples "+str(updateSubstMatrixEveryThisSamples[s])
		+" --thresholdDiffForUpdate "+str(thresholdDiffForUpdate[s])+" --thresholdFoldChangeUpdate "+str(thresholdFoldChangeUpdate[s])+" --minBLen "+str(minBLen[s])+" --maxBLen "+str(maxBLen[s])+nonStrictTopologyStopRules[s]+" --calculateLKfinalTree --overwrite --thresholdLogLKconsecutivePlacement "+str(thresholdLogLKconsecutivePlacement[s])+" --output "
		+folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE_"+speedsMAPLE[s]+"\n"+"done\n\n") #
		#+folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE_"+speedsMAPLE[s]+"\n"+"done\n\n") #

		if j<20000:
			file2.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"cd "+folderName+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(j)+" -o "
			+folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE_"+speedsMAPLE[s]+"_evaluation_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE_"
			+speedsMAPLE[s]+"_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderName+str(j)+"subsamples/phylipFile_seed\"$i\"_"
			+str(j)+"samples.phy -st DNA -te "+folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE_"+speedsMAPLE[s]+"_tree.tree -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsMAPLE[s]
			+"_MAPLEEvaluation -m GTR -quiet -nt 1 -keep-ident -redo -blmin 0.00000000501 \n"+"done\n\n")

	file.write("\n\n")
	file2.write("\n\n")
file.close()
file2.close()

mapleVersions=["0.1.1","0.1.3"]

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/MAPLE_testingVersion.sh","w")
#for j in [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000]:
for j in [100000]:
	for s in range(len(mapleVersions)):
		
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"rm -f "
		+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE"+mapleVersions[s]+"_GTR_tree.tree || true\n")
		#+"done\n\n")

		#file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/6))+" -o "
		file.write("bsub -M "+str(int(400+j/15))+" -o "
		+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE"+mapleVersions[s]+"_GTR_console_output.txt -e "
		+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE"+mapleVersions[s]+"_GTR_console_error.txt /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/MAPLEv"+mapleVersions[s]+".py --reference /nfs/research/goldman/demaio/fastLK/simulations/MN908947.3.fasta --input "
		+folderNameSimu+str(j)+"subsamples/diffFile_repeat\"$i\"_"+str(j)+"samples.txt --calculateLKfinalTree --overwrite --output "
		+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE"+mapleVersions[s]+"_GTR\n") # 
		#+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE_GTR_"+speedsMAPLE[s]+"\n"+"done\n\n") # 

		#file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"rm -f "
		#+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE_Ns_"+speedsMAPLE[s]+"_tree.tree || true\n"+"done\n\n")
		file.write("rm -f "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE"+mapleVersions[s]+"_Ns_tree.tree || true\n")

		#file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/6))+" -o "
		file.write("bsub -M "+str(int(400+j/15))+" -o "
		+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE"+mapleVersions[s]+"_Ns_console_output.txt -e "
		+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE"+mapleVersions[s]+"_Ns_console_error.txt /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/MAPLEv"+mapleVersions[s]+".py --reference /nfs/research/goldman/demaio/fastLK/simulations/MN908947.3.fasta --input "
		+folderNameSimu+str(j)+"subsamples/diffFile_repeat\"$i\"_"+str(j)+"samples_Ns.txt --calculateLKfinalTree --overwrite --output "
		+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE"+mapleVersions[s]+"_Ns\n") # 
		#+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE_Ns_"+speedsMAPLE[s]+"\n"+"done\n\n") # 

		#file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"rm -f "
		#+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE_4cat_"+speedsMAPLE[s]+"_tree.tree || true\n"+"done\n\n")
		file.write("rm -f "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE"+mapleVersions[s]+"_4cat_tree.tree || true\n")

		#file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/6))+" -o "
		file.write("bsub -M "+str(int(400+j/15))+" -o "
		+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE"+mapleVersions[s]+"_4cat_console_output.txt -e "
		+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE"+mapleVersions[s]+"_4cat_console_error.txt /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/MAPLEv"+mapleVersions[s]+".py --reference /nfs/research/goldman/demaio/fastLK/simulations/MN908947.3.fasta --input "
		+folderNameSimu+str(j)+"subsamples/diffFile_4cat_repeat\"$i\"_"+str(j)+"samples.txt --calculateLKfinalTree --overwrite --output "
		+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE"+mapleVersions[s]+"_4cat\n") #
		#+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/MAPLE_4cat_"+speedsMAPLE[s]+"\n"+"done\n\n") #

		#file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"rm -f "
		#+folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE_"+speedsMAPLE[s]+"_tree.tree || true\n"+"done\n\n")
		file.write("rm -f "+folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE"+mapleVersions[s]+"_tree.tree || true\n")


		#if j<20000:
		#file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/6))+" -o "
		file.write("bsub -M "+str(int(400+j/15))+" -o "
		+folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE"+mapleVersions[s]+"_console_output.txt -e "
		+folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE"+mapleVersions[s]+"_console_error.txt /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 /nfs/research/goldman/demaio/fastLK/code/MAPLEv"+mapleVersions[s]+".py --reference /nfs/research/goldman/demaio/fastLK/realData/2021-03-31_unmasked_differences_reduced.txt_consensus.fa --input "
		+folderName+str(j)+"subsamples/diffFile_seed\"$i\"_"+str(j)+"samples.txt --calculateLKfinalTree --overwrite --output "
		+folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE"+mapleVersions[s]+"\n"+"done\n\n") #
		#+folderName+str(j)+"subsamples/output_repl\"$i\"/MAPLE_"+speedsMAPLE[s]+"\n"+"done\n\n") #

	file.write("\n\n")
file.close()









#UShER
speedsUsher=["medium"]
file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitUsher.sh","w")
for j in [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000]:
	for s in range(len(speedsUsher)):
		file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/5))+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/usher_"+speedsUsher[s]+"_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/usher_"+speedsUsher[s]+"_console_error.txt usher -v "+folderName+str(j)+"subsamples/vcfFile_seed\"$i\"_"+str(j)+"samples.vcf  -d "+folderName+str(j)+"subsamples/output_repl\"$i\"/usher_output/ -T 1 -t "+folderName+str(j)+"subsamples/initialTree_seed\"$i\"_"+str(j)+"samples.tree  \n"+"done\n\n")
		
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/5))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_"+speedsUsher[s]+"_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_"+speedsUsher[s]+"_console_error.txt usher -v "+folderNameSimu+str(j)+"subsamples/vcfFile_repeat\"$i\"_"+str(j)+"samples.vcf  -d "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_output/ -T 1 -t "+folderNameSimu+str(j)+"subsamples/initialTree_repeat\"$i\"_"+str(j)+"samples.tree \n"+"done\n\n")

		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/5))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_4cat_"+speedsUsher[s]+"_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_4cat_"+speedsUsher[s]+"_console_error.txt usher -v "+folderNameSimu+str(j)+"subsamples/vcfFile_4cat_repeat\"$i\"_"+str(j)+"samples.vcf  -d "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_output_4cat/ -T 1 -t "+folderNameSimu+str(j)+"subsamples/initialTree_4cat_repeat\"$i\"_"+str(j)+"samples.tree \n"+"done\n\n")
	
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/5))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_"+speedsUsher[s]+"_Ns_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_"+speedsUsher[s]+"_Ns_console_error.txt usher -v "+folderNameSimu+str(j)+"subsamples/vcfFile_repeat\"$i\"_"+str(j)+"samples_Ns.vcf  -d "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_output_Ns/ -T 1 -t "+folderNameSimu+str(j)+"subsamples/initialTree_repeat\"$i\"_"+str(j)+"samples_Ns.tree \n"+"done\n\n")
	
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitUsherEvaluation.sh","w")
for j in [100, 200, 500, 1000, 2000, 5000, 10000, 20000]:
	for s in range(len(speedsUsher)):
		file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"cd "+folderName+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(j)+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/usher_"+speedsUsher[s]+"_evaluation_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/usher_"+speedsUsher[s]+"_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderName+str(j)+"subsamples/phylipFile_seed\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderName+str(j)+"subsamples/output_repl\"$i\"/usher_output/final-tree.nh -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsUsher[s]+"_usherEvaluation -m GTR -quiet -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
		
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*3.0))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_4cat_"+speedsUsher[s]+"_evaluation_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_4cat_"+speedsUsher[s]+"_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_4cat_repeat\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_output_4cat/final-tree.nh -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsUsher[s]+"_usherEvaluation_4cat -m GTR+G -quiet -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
		
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*1.7))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_"+speedsUsher[s]+"_evaluation_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_"+speedsUsher[s]+"_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_repeat\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_output/final-tree.nh -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsUsher[s]+"_usherEvaluation -m GTR -quiet -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
		
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*1.7))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_"+speedsUsher[s]+"_Ns_evaluation_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_"+speedsUsher[s]+"_Ns_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_repeat\"$i\"_"+str(j)+"samples_Ns.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_output_Ns/final-tree.nh -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsUsher[s]+"_Ns_usherEvaluation -m GTR -quiet -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
	file.write("\n\n")
file.close()

speedsMAT=["medium"]
file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitMatOptimize.sh","w")
for j in [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000]:
	for s in range(len(speedsMAT)):
		file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/5))+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_console_error.txt matOptimize -n -v "+folderName+str(j)+"subsamples/vcfFile_seed\"$i\"_"+str(j)+"samples.vcf  -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_output.pb -T 1 -t "+folderName+str(j)+"subsamples/output_repl\"$i\"/usher_output/final-tree.nh  \n"+"done\n\n")
		
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/5))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_console_error.txt matOptimize -n -v "+folderNameSimu+str(j)+"subsamples/vcfFile_repeat\"$i\"_"+str(j)+"samples.vcf  -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_output.pb -T 1 -t "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_output/final-tree.nh \n"+"done\n\n")

		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/5))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_4cat_"+speedsMAT[s]+"_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_4cat_"+speedsMAT[s]+"_console_error.txt matOptimize -n -v "+folderNameSimu+str(j)+"subsamples/vcfFile_4cat_repeat\"$i\"_"+str(j)+"samples.vcf  -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_output_4cat.pb -T 1 -t "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_output_4cat/final-tree.nh \n"+"done\n\n")
	
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/5))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_Ns_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_Ns_console_error.txt matOptimize -n -v "+folderNameSimu+str(j)+"subsamples/vcfFile_repeat\"$i\"_"+str(j)+"samples_Ns.vcf  -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_output_Ns.pb -T 1 -t "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/usher_output_Ns/final-tree.nh \n"+"done\n\n")
	
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitMatOptimize_convertOutputTree.sh","w")
for j in [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000]:
	for s in range(len(speedsMAT)):
		file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/5))+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/matOptimizeConvert_"+speedsMAT[s]+"_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/matOptimizeConvert_"+speedsMAT[s]+"_console_error.txt matUtils extract -i "+folderName+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_output.pb -T 1 -d "+folderName+str(j)+"subsamples/output_repl\"$i\"/ -t matOptimize_"+speedsMAT[s]+"_output.nh  \n"+"done\n\n")
		
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/5))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimizeConvert_"+speedsMAT[s]+"_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimizeConvert_"+speedsMAT[s]+"_console_error.txt matUtils extract -i "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_output.pb -T 1 -d "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/ -t matOptimize_"+speedsMAT[s]+"_output.nh \n"+"done\n\n")

		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/5))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimizeConvert_4cat_"+speedsMAT[s]+"_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimizeConvert_4cat_"+speedsMAT[s]+"_console_error.txt matUtils extract -i "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_output_4cat.pb -T 1 -d "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/ -t matOptimize_"+speedsMAT[s]+"_output_4cat.nh \n"+"done\n\n")
	
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(400+j/5))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimizeConvert_"+speedsMAT[s]+"_Ns_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimizeConvert_"+speedsMAT[s]+"_Ns_console_error.txt matUtils extract -i "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_output_Ns.pb -T 1 -d "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/ -t matOptimize_"+speedsMAT[s]+"_output_Ns.nh \n"+"done\n\n")
	
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitMatOptimizeEvaluation.sh","w")
for j in [100, 200, 500, 1000, 2000, 5000, 10000, 20000]:
	for s in range(len(speedsMAT)):
		file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"cd "+folderName+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*1.7))+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_evaluation_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderName+str(j)+"subsamples/phylipFile_seed\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderName+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_output.nh -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsMAT[s]+"_matOptimizeEvaluation -m GTR -quiet -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
		
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*3.0))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_4cat_"+speedsMAT[s]+"_evaluation_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_4cat_"+speedsMAT[s]+"_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_4cat_repeat\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_output_4cat.nh -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsMAT[s]+"_matOptimizeEvaluation_4cat -m GTR+G -quiet -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
		
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*1.7))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_evaluation_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_repeat\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_output.nh -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsMAT[s]+"_matOptimizeEvaluation -m GTR -quiet -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
		
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*1.7))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_Ns_evaluation_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_Ns_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_repeat\"$i\"_"+str(j)+"samples_Ns.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/matOptimize_"+speedsMAT[s]+"_output_Ns.nh -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsMAT[s]+"_Ns_matOptimizeEvaluation -m GTR -quiet -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
	file.write("\n\n")
file.close()






#IQtree
speedsIQtree=["low","medium","fast"]
speedString=["-pers 0.1 -nstop 500 -blmin 0.000000005","-blmin 0.000000005","-fast"]
file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitIQtree.sh","w")
#for j in [100, 200, 500, 1000, 2000, 5000, 10000, 20000]:
for j in [100, 200, 500, 1000, 2000, 5000, 10000]:
	for s in range(len(speedsIQtree)):
		file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"cd "+folderName+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(j*1.7))+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/IQtree_"+speedsIQtree[s]+"_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/IQtree_"+speedsIQtree[s]+"_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderName+str(j)+"subsamples/phylipFile_seed\"$i\"_"+str(j)+"samples.phy -st DNA -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsIQtree[s]+"_iqtree -m GTR -quiet -redo -nt 1 "+speedString[s]+" \n"+"done\n\n")
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitIQtreeEvaluation.sh","w")
for j in [100, 200, 500, 1000, 2000, 5000, 10000]:
	for s in range(len(speedsIQtree)):
		file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"cd "+folderName+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(j)+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/IQtree_"+speedsIQtree[s]+"_evaluation_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/IQtree_"+speedsIQtree[s]+"_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderName+str(j)+"subsamples/phylipFile_seed\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderName+str(j)+"subsamples/output_repl\"$i\"/"+str(j)+"subsamples_repl\"$i\"_"+speedsIQtree[s]+"_iqtree.treefile -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsIQtree[s]+"_IQtreeEvaluation -m GTR -quiet -redo -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
	file.write("\n\n")
file.close()



#IQtree starting from fastLK tree
speedsIQtree=["slow","medium","fast"]
speedString=["-pers 0.1 -nstop 500 -blmin 0.000000005","-blmin 0.000000005","-fast"]
file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitIQtree_startFromFastLK.sh","w")
for j in [100, 200, 500, 1000, 2000, 5000, 10000]:
#for j in [100, 200, 500, 1000, 2000]:
	for s in range(len(speedsIQtree)):
		file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"cd "+folderName+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(j*1.7))+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/IQtree_"+speedsIQtree[s]+"_startFromFastLK_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/IQtree_"+speedsIQtree[s]+"_startFromFastLK_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderName+str(j)+"subsamples/phylipFile_seed\"$i\"_"+str(j)+"samples.phy -st DNA -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsIQtree[s]+"_startFromFastLK_iqtree -m GTR -quiet -redo -nt 1 "+speedString[s]+" -t "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastLK_slow2_tree.tree \n"+"done\n\n")
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitIQtree_startFromFastLKEvaluation.sh","w")
for j in [100, 200, 500, 1000, 2000, 5000, 10000]:
	for s in range(len(speedsIQtree)):
		file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"cd "+folderName+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(j)+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/IQtree_"+speedsIQtree[s]+"_startFromFastLK_evaluation_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/IQtree_"+speedsIQtree[s]+"_startFromFastLK_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderName+str(j)+"subsamples/phylipFile_seed\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderName+str(j)+"subsamples/output_repl\"$i\"/"+str(j)+"subsamples_repl\"$i\"_"+speedsIQtree[s]+"_startFromFastLK_iqtree.treefile -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsIQtree[s]+"_startFromFastLK_IQtreeEvaluation -m GTR -quiet -redo -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
	file.write("\n\n")
file.close()


#IQtree on simulated data
file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitIQtree_simulations.sh","w")
#for j in [100, 200, 500, 1000, 2000, 5000, 10000, 20000]:
for j in [1000, 2000]:
	for s in range(len(speedsIQtree)):
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*1.7))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/IQtree_"+speedsIQtree[s]+"_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/IQtree_"+speedsIQtree[s]+"_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_repeat\"$i\"_"+str(j)+"samples.phy -st DNA -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsIQtree[s]+"_iqtree -m GTR -quiet -redo -nt 1 "+speedString[s]+" \n"+"done\n\n")# -redo
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*3.0))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/IQtree_"+speedsIQtree[s]+"_4cat_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/IQtree_"+speedsIQtree[s]+"_4cat_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_4cat_repeat\"$i\"_"+str(j)+"samples.phy -st DNA -pre "+str(j)+"subsamples_repl\"$i\"_4cat_"+speedsIQtree[s]+"_iqtree -m GTR+G -quiet -redo -nt 1 "+speedString[s]+" \n"+"done\n\n")
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitIQtree_simulations_Ns.sh","w")
for j in [100, 200, 500, 1000, 2000, 5000, 10000, 20000]:
	for s in range(len(speedsIQtree)):
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*1.7))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/IQtree_"+speedsIQtree[s]+"_Ns_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/IQtree_"+speedsIQtree[s]+"_Ns_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_repeat\"$i\"_"+str(j)+"samples_Ns.phy -st DNA -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsIQtree[s]+"_Ns_iqtree -m GTR -quiet -nt 1 "+speedString[s]+" \n"+"done\n\n")# -redo
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitIQtree_simulationsEvaluation.sh","w")
#for j in [100, 200, 500, 1000, 2000, 5000, 10000, 20000]:
for j in [1000, 2000]:
	for s in range(len(speedsIQtree)):
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*1.7))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/IQtree_"+speedsIQtree[s]+"_evaluation_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/IQtree_"+speedsIQtree[s]+"_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_repeat\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/"+str(j)+"subsamples_repl\"$i\"_"+speedsIQtree[s]+"_iqtree.treefile -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsIQtree[s]+"_IQtreeEvaluation -m GTR -quiet -redo -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*3.0))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/IQtree_"+speedsIQtree[s]+"_4cat_evaluation_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/IQtree_"+speedsIQtree[s]+"_4cat_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_4cat_repeat\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/"+str(j)+"subsamples_repl\"$i\"_4cat_"+speedsIQtree[s]+"_iqtree.treefile -pre "+str(j)+"subsamples_repl\"$i\"_4cat_"+speedsIQtree[s]+"_IQtreeEvaluation -m GTR+G -quiet -redo -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitIQtree_simulationsEvaluation_Ns.sh","w")
for j in [100, 200, 500, 1000, 2000, 5000, 10000, 20000]:
	for s in range(len(speedsIQtree)):
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*1.7))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/IQtree_"+speedsIQtree[s]+"_Ns_evaluation_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/IQtree_"+speedsIQtree[s]+"_Ns_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_repeat\"$i\"_"+str(j)+"samples_Ns.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/"+str(j)+"subsamples_repl\"$i\"_"+speedsIQtree[s]+"_Ns_iqtree.treefile -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsIQtree[s]+"_Ns_IQtreeEvaluation -m GTR -quiet -redo -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
	file.write("\n\n")
file.close()






#FastTree2 on simulated data
speedsFastTree=["slow","medium","fast","fastest"]
speedString=["-spr 4 -mlacc 2 -slownni","-spr 4","","-fastest"]
file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitFastTree_simulations.sh","w")
#for j in [100, 200, 500, 1000, 2000, 5000, 10000, 20000]:
for j in [1000, 2000]:
	for s in range(len(speedsFastTree)):
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(1.7*j+500))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastTree_"+speedsFastTree[s]+"_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastTree_"+speedsFastTree[s]+"_console_error.txt /hps/software/users/goldman/FastTree2/FastTree -quiet -nosupport -nt -gtr -nocat "+speedString[s]+" "+folderNameSimu+str(j)+"subsamples/fastaFile_repeat\"$i\"_"+str(j)+"samples.fa \n"+"done\n\n")
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(1.7*j+500))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastTree_"+speedsFastTree[s]+"_4cat_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastTree_"+speedsFastTree[s]+"_4cat_console_error.txt /hps/software/users/goldman/FastTree2/FastTree -quiet -nosupport -nt -gtr -cat 4 "+speedString[s]+" "+folderNameSimu+str(j)+"subsamples/fastaFile_4cat_repeat\"$i\"_"+str(j)+"samples.fa \n"+"done\n\n")
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitFastTree_simulationsEvaluation.sh","w")
#for j in [100, 200, 500, 1000, 2000, 5000, 10000, 20000]:
for j in [1000, 2000]:
	for s in range(len(speedsFastTree)):
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(1.7*j+500))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastTree_"+speedsFastTree[s]+"_evaluation_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastTree_"+speedsFastTree[s]+"_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_repeat\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastTree_"+speedsFastTree[s]+"_tree.tree -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsFastTree[s]+"_fastTreeEvaluation -m GTR -quiet -redo -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(3.0*j+500))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastTree_"+speedsFastTree[s]+"_4cat_evaluation_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastTree_"+speedsFastTree[s]+"_4cat_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_4cat_repeat\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastTree_"+speedsFastTree[s]+"_4cat_tree.tree -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsFastTree[s]+"_4cat_fastTreeEvaluation -m GTR+G -quiet -redo -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitFastTree_simulations_Ns.sh","w")
for j in [100, 200, 500, 1000, 2000, 5000, 10000, 20000]:
	for s in range(len(speedsFastTree)):
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(1.7*j+500))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastTree_"+speedsFastTree[s]+"_Ns_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastTree_"+speedsFastTree[s]+"_Ns_console_error.txt /hps/software/users/goldman/FastTree2/FastTree -quiet -nosupport -nt -gtr -nocat "+speedString[s]+" "+folderNameSimu+str(j)+"subsamples/fastaFile_repeat\"$i\"_"+str(j)+"samples_Ns.fa \n"+"done\n\n")
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitFastTree_simulationsEvaluation_Ns.sh","w")
for j in [100, 200, 500, 1000, 2000, 5000, 10000, 20000]:
#for j in [20000]:
	for s in range(len(speedsFastTree)):
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(1.7*j+500))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastTree_"+speedsFastTree[s]+"_Ns_evaluation_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastTree_"+speedsFastTree[s]+"_Ns_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_repeat\"$i\"_"+str(j)+"samples_Ns.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/fastTree_"+speedsFastTree[s]+"_Ns_tree.tree -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsFastTree[s]+"_Ns_fastTreeEvaluation -m GTR -quiet -redo -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
	file.write("\n\n")
file.close()

#FastTree2
speedsFastTree=["slow","medium","fast","fastest"]
speedString=["-spr 4 -mlacc 2 -slownni","-spr 4","","-fastest"]
file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitFastTree.sh","w")
for j in [100, 200, 500, 1000, 2000, 5000, 10000, 20000]:
	for s in range(len(speedsFastTree)):
		file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"cd "+folderName+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(1.7*j+500))+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastTree_"+speedsFastTree[s]+"_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastTree_"+speedsFastTree[s]+"_console_error.txt /hps/software/users/goldman/FastTree2/FastTree -quiet -nosupport -nt -gtr -nocat "+speedString[s]+" "+folderName+str(j)+"subsamples/fastaFile_seed\"$i\"_"+str(j)+"samples.fa \n"+"done\n\n")
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitFastTreeEvaluation.sh","w")
for j in [100, 200, 500, 1000, 2000, 5000, 10000, 20000]:
#for j in [20000]:
	for s in range(len(speedsFastTree)):
		file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"cd "+folderName+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(j)+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastTree_"+speedsFastTree[s]+"_evaluation_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastTree_"+speedsFastTree[s]+"_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderName+str(j)+"subsamples/phylipFile_seed\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastTree_"+speedsFastTree[s]+"_tree.tree -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsFastTree[s]+"_fastTreeEvaluation -m GTR -quiet -redo -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
	file.write("\n\n")
file.close()





#FastTree2 starting from fastLK tree
speedsFastTree=["slow","medium","fast","fastest"]
speedString=["-spr 4 -mlacc 2 -slownni","-spr 4","","-fastest"]
file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitFastTree_startFromFastLK.sh","w")
for j in [100, 200, 500, 1000, 2000, 5000, 10000]:
#for j in [100, 200, 500, 1000, 2000]:
	for s in range(len(speedsFastTree)):
		file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"cd "+folderName+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(1.7*j+500))+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastTree_"+speedsFastTree[s]+"_startFromFastLK_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastTree_"+speedsFastTree[s]+"_startFromFastLK_console_error.txt /hps/software/users/goldman/FastTree2/FastTree -quiet -nosupport -nt -gtr -nocat -intree "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastLK_slow2_tree.tree "+speedString[s]+" "+folderName+str(j)+"subsamples/fastaFile_seed\"$i\"_"+str(j)+"samples.fa \n"+"done\n\n")
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitFastTree_startFromFastLKEvaluation.sh","w")
for j in [100, 200, 500, 1000, 2000, 5000, 10000]:
	for s in range(len(speedsFastTree)):
		file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"cd "+folderName+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(j)+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastTree_"+speedsFastTree[s]+"_startFromFastLK_evaluation_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastTree_"+speedsFastTree[s]+"_startFromFastLK_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderName+str(j)+"subsamples/phylipFile_seed\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastTree_"+speedsFastTree[s]+"_startFromFastLK_tree.tree -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsFastTree[s]+"_startFromFastLK_fastTreeEvaluation -m GTR -quiet -redo -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
	file.write("\n\n")
file.close()











#raxmlHPC
speedsRaxml=["slow","fast"]
speedString=["","-D"]
file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitRaxml.sh","w")
file.write("module load raxml-8.2.11-gcc-9.3.0-mjwrm3x \n")
for j in [100, 200, 500, 1000, 2000, 5000]:
	for s in range(len(speedsRaxml)):
		file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"bsub -M "+str(int(j+500))+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/raxml_"+speedsRaxml[s]+"_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/raxml_"+speedsRaxml[s]+"_console_error.txt raxmlHPC "+speedString[s]+" -s "+folderName+str(j)+"subsamples/phylipFile_seed\"$i\"_"+str(j)+"samples.phy -p 1 -F -c 1 -m GTRCAT -V -n raxml_"+speedsRaxml[s]+"_output -w "+folderName+str(j)+"subsamples/output_repl\"$i\"/ \n"+"done\n\n")
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitRaxmlEvaluation.sh","w")
for j in [100, 200, 500, 1000, 2000, 5000]:
	for s in range(len(speedsRaxml)):
		file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"cd "+folderName+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(j)+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/raxml_"+speedsRaxml[s]+"_evaluation_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/raxml_"+speedsRaxml[s]+"_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderName+str(j)+"subsamples/phylipFile_seed\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderName+str(j)+"subsamples/output_repl\"$i\"/RAxML_result.raxml_"+speedsRaxml[s]+"_output -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsRaxml[s]+"_raxmlEvaluation -m GTR -quiet -redo -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
	file.write("\n\n")
file.close()




#raxmlHPC starting from fastLK tree
speedsRaxml=["slow","fast"]
speedString=["","-D"]
file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitRaxml_startFromFastLK.sh","w")
file.write("module load raxml-8.2.11-gcc-9.3.0-mjwrm3x \n")
for j in [100, 200, 500, 1000, 2000, 5000]:
#for j in [100, 200, 500]:
	for s in range(len(speedsRaxml)):
		file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"bsub -M "+str(int(j+500))+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/raxml_"+speedsRaxml[s]+"_startFromFastLK_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/raxml_"+speedsRaxml[s]+"_startFromFastLK_console_error.txt raxmlHPC "+speedString[s]+" -s "+folderName+str(j)+"subsamples/phylipFile_seed\"$i\"_"+str(j)+"samples.phy -p 1 -F -c 1 -m GTRCAT -V -t "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastLK_slow2_tree.tree -n raxml_"+speedsRaxml[s]+"_startFromFastLK_output2 -w "+folderName+str(j)+"subsamples/output_repl\"$i\"/ \n"+"done\n\n")
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitRaxml_startFromFastLKEvaluation.sh","w")
for j in [100, 200, 500, 1000, 2000, 5000]:
	for s in range(len(speedsRaxml)):
		file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"cd "+folderName+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(j)+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/raxml_"+speedsRaxml[s]+"_startFromFastLK_evaluation_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/raxml_"+speedsRaxml[s]+"_startFromFastLK_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderName+str(j)+"subsamples/phylipFile_seed\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderName+str(j)+"subsamples/output_repl\"$i\"/RAxML_result.raxml_"+speedsRaxml[s]+"_startFromFastLK_output2 -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsRaxml[s]+"_startFromFastLK_raxmlEvaluation -m GTR -quiet -redo -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
	file.write("\n\n")
file.close()



#raxmlHPC on simulated data
speedsRaxml=["slow","fast"]
speedString=["","-D"]
file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitRaxml_simulations.sh","w")
file.write("module load raxml-8.2.11-gcc-9.3.0-mjwrm3x \n")
#for j in [100, 200, 500, 1000, 2000, 5000]:
for j in [1000, 2000]:
	for s in range(len(speedsRaxml)):
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"rm "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/*raxml_"+speedsRaxml[s]+"_output* \n\t"+"bsub -M "+str(int(2*j+500))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxml_"+speedsRaxml[s]+"_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxml_"+speedsRaxml[s]+"_console_error.txt raxmlHPC "+speedString[s]+" -s "+folderNameSimu+str(j)+"subsamples/phylipFile_repeat\"$i\"_"+str(j)+"samples.phy -p 1 -c 1 -m GTRCAT -V -n raxml_"+speedsRaxml[s]+"_output -w "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/ \n"+"done\n\n")
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"rm "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/*raxml_"+speedsRaxml[s]+"_4cat_output \n\t"+"bsub -M "+str(int(4*j+500))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxml_"+speedsRaxml[s]+"_4cat_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxml_"+speedsRaxml[s]+"_4cat_console_error.txt raxmlHPC "+speedString[s]+" -s "+folderNameSimu+str(j)+"subsamples/phylipFile_4cat_repeat\"$i\"_"+str(j)+"samples.phy -p 1 -c 4 -m GTRCAT -n raxml_"+speedsRaxml[s]+"_4cat_output -w "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/ \n"+"done\n\n")
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitRaxml_simulationsEvaluation.sh","w")
#for j in [100, 200, 500, 1000, 2000, 5000]:
for j in [1000, 2000]:
	for s in range(len(speedsRaxml)):
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*1.7))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxml_"+speedsRaxml[s]+"_evaluation_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxml_"+speedsRaxml[s]+"_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_repeat\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/RAxML_result.raxml_"+speedsRaxml[s]+"_output -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsRaxml[s]+"_raxmlEvaluation -m GTR -quiet -redo -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*3.0))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxml_"+speedsRaxml[s]+"_4cat_evaluation_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxml_"+speedsRaxml[s]+"_4cat_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_4cat_repeat\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/RAxML_result.raxml_"+speedsRaxml[s]+"_4cat_output -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsRaxml[s]+"_4cat_raxmlEvaluation -m GTR+G -quiet -redo -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitRaxml_simulations_Ns.sh","w")
file.write("module load raxml-8.2.11-gcc-9.3.0-mjwrm3x \n")
#for j in [100, 200, 500, 1000, 2000, 5000]:
for j in [ 5000]:
	for s in range(len(speedsRaxml)):
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"rm "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/*raxml_"+speedsRaxml[s]+"_Ns_output* \n\t"+"bsub -M "+str(int(2*j+500))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxml_"+speedsRaxml[s]+"_Ns_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxml_"+speedsRaxml[s]+"_Ns_console_error.txt raxmlHPC "+speedString[s]+" -s "+folderNameSimu+str(j)+"subsamples/phylipFile_repeat\"$i\"_"+str(j)+"samples_Ns.phy -p 1 -c 1 -m GTRCAT -V -n raxml_"+speedsRaxml[s]+"_Ns_output -w "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/ \n"+"done\n\n")
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitRaxml_simulationsEvaluation_Ns.sh","w")
#for j in [100, 200, 500, 1000, 2000, 5000]:
for j in [ 5000]:
	for s in range(len(speedsRaxml)):
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*1.7))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxml_"+speedsRaxml[s]+"_Ns_evaluation_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxml_"+speedsRaxml[s]+"_Ns_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_repeat\"$i\"_"+str(j)+"samples_Ns.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/RAxML_result.raxml_"+speedsRaxml[s]+"_Ns_output -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsRaxml[s]+"_Ns_raxmlEvaluation -m GTR -quiet -redo -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
	file.write("\n\n")
file.close()





#raxmlNG on simulated data
speedsRaxmlNG=["slow","fast"]
speedString=["--blmin 0.000000005 --tree pars\{3\}","--blmin 0.000005 --tree pars\{1\}"]
file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitRaxmlNG_simulations.sh","w")
file.write("module load raxml-ng-1.0.2-gcc-9.3.0-uicuzej \n")
#for j in [100, 200, 500, 1000, 2000, 5000]:
for j in [1000, 2000]:
	for s in range(len(speedsRaxmlNG)):
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(j+500))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_console_error.txt raxml-ng --search "+speedString[s]+" --msa "+folderNameSimu+str(j)+"subsamples/fastaFile_repeat\"$i\"_"+str(j)+"samples.fa --msa-format FASTA --data-type DNA  --redo --prefix "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_output --model GTR --threads 1 \n"+"done\n\n")
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(j+500))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_4cat_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_4cat_console_error.txt raxml-ng --search "+speedString[s]+" --msa "+folderNameSimu+str(j)+"subsamples/fastaFile_4cat_repeat\"$i\"_"+str(j)+"samples.fa --msa-format FASTA --data-type DNA  --redo --prefix "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_4cat_output --model GTR+G --threads 1 \n"+"done\n\n")
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitRaxmlNG_simulationsEvaluation.sh","w")
#for j in [100, 200, 500, 1000, 2000, 5000]:
for j in [1000, 2000]:
	for s in range(len(speedsRaxmlNG)):
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*1.7))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_evaluation_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_repeat\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_output.raxml.bestTree -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsRaxmlNG[s]+"_raxmlNGEvaluation -m GTR -quiet -redo -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*4.0))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_4cat_evaluation_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_4cat_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_4cat_repeat\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_4cat_output.raxml.bestTree -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsRaxmlNG[s]+"_4cat_raxmlNGEvaluation -m GTR+G -quiet -redo -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitRaxmlNG_simulations_Ns.sh","w")
file.write("module load raxml-ng-1.0.2-gcc-9.3.0-uicuzej \n")
for j in [100, 200, 500, 1000, 2000, 5000]:
	for s in range(len(speedsRaxmlNG)):
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"bsub -M "+str(int(j+500))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_Ns_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_Ns_console_error.txt raxml-ng --search "+speedString[s]+" --msa "+folderNameSimu+str(j)+"subsamples/fastaFile_repeat\"$i\"_"+str(j)+"samples_Ns.fa --msa-format FASTA --data-type DNA  --redo --prefix "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_Ns_output --model GTR --threads 1 \n"+"done\n\n")
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitRaxmlNG_simulationsEvaluation_Ns.sh","w")
for j in [100, 200, 500, 1000, 2000, 5000]:
	for s in range(len(speedsRaxmlNG)):
		file.write("for i in $(seq 1 "+str(nRepeatsSimu)+")\n"+"do \n\t"+"cd "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(int(500+j*1.7))+" -o "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_Ns_evaluation_console_output.txt -e "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_Ns_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderNameSimu+str(j)+"subsamples/phylipFile_repeat\"$i\"_"+str(j)+"samples_Ns.phy -st DNA -te "+folderNameSimu+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_Ns_output.raxml.bestTree -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsRaxmlNG[s]+"_Ns_raxmlNGEvaluation -m GTR -quiet -redo -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
	file.write("\n\n")
file.close()


#raxmlNG
speedsRaxmlNG=["slow","fast"]
#speedsRaxmlNG=["slow","medium","fast"]
#speedsRaxmlNG=["fast"]
#speedString=["--blmin 0.000000005","--tree pars\{3\},rand\{3\} --blmin 0.000000005","--tree pars\{1\},rand\{1\}"]
speedString=["--blmin 0.000000005 --tree pars\{3\}","--blmin 0.000005 --tree pars\{1\}"]
#speedString=["--tree pars\{1\},rand\{1\}"]
file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitRaxmlNG.sh","w")
file.write("module load raxml-ng-1.0.2-gcc-9.3.0-uicuzej \n")
for j in [100, 200, 500, 1000, 2000, 5000]:
	for s in range(len(speedsRaxmlNG)):
		file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"bsub -M "+str(int(j+500))+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_console_error.txt raxml-ng --search "+speedString[s]+" --msa "+folderName+str(j)+"subsamples/fastaFile_seed\"$i\"_"+str(j)+"samples.fa --msa-format FASTA --data-type DNA  --redo --prefix "+folderName+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_output --model GTR --threads 1 \n"+"done\n\n")
	file.write("\n\n")
file.close()

#speedString=["--blmin","--tree","fast_output"]
file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitRaxmlNGEvaluation.sh","w")
for j in [100, 200, 500, 1000, 2000, 5000]:
#for j in [2000, 5000]:
	for s in range(len(speedsRaxmlNG)):
		file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"cd "+folderName+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(j)+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_evaluation_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderName+str(j)+"subsamples/phylipFile_seed\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderName+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_output.raxml.bestTree -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsRaxmlNG[s]+"_raxmlNGEvaluation -m GTR -quiet -redo -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
	file.write("\n\n")
file.close()




#raxmlNG starting from fastLK tree
speedsRaxmlNG=["slow","fast"]
#speedsRaxmlNG=["fast"]
speedString=["--blmin 0.000000005 ","--blmin 0.000005 "]
#speedString=["--tree pars\{1\},rand\{1\}"]
file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitRaxmlNG_startFromFastLK.sh","w")
file.write("module load raxml-ng-1.0.2-gcc-9.3.0-uicuzej \n")
for j in [100, 200, 500, 1000, 2000, 5000]:
#for j in [100, 200, 500]:
	for s in range(len(speedsRaxmlNG)):
		file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"bsub -M "+str(int(j+500))+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_startFromFastLK_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_startFromFastLK_console_error.txt raxml-ng --search "+speedString[s]+" --msa "+folderName+str(j)+"subsamples/fastaFile_seed\"$i\"_"+str(j)+"samples.fa --msa-format FASTA --data-type DNA --redo --tree "+folderName+str(j)+"subsamples/output_repl\"$i\"/fastLK_slow2_tree.tree --prefix "+folderName+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_startFromFastLK_output --model GTR --threads 1 \n"+"done\n\n")
	file.write("\n\n")
file.close()

file=open("/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/submitRaxmlNG_startFromFastLKEvaluation.sh","w")
for j in [100, 200, 500, 1000, 2000, 5000]:
	for s in range(len(speedsRaxmlNG)):
		file.write("for i in $(seq 1 "+str(nRepeats)+")\n"+"do \n\t"+"cd "+folderName+str(j)+"subsamples/output_repl\"$i\"\n\t"+"bsub -M "+str(j)+" -o "+folderName+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_startFromFastLK_evaluation_console_output.txt -e "+folderName+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_startFromFastLK_evaluation_console_error.txt /hps/software/users/goldman/iqtree-2.1.3-Linux/bin/iqtree2 -s "+folderName+str(j)+"subsamples/phylipFile_seed\"$i\"_"+str(j)+"samples.phy -st DNA -te "+folderName+str(j)+"subsamples/output_repl\"$i\"/raxmlNG_"+speedsRaxmlNG[s]+"_startFromFastLK_output.raxml.bestTree -pre "+str(j)+"subsamples_repl\"$i\"_"+speedsRaxmlNG[s]+"_raxmlNG_startFromFastLKEvaluation -m GTR -quiet -redo -nt 1 -keep-ident -blmin 0.000000005 \n"+"done\n\n")
	file.write("\n\n")
file.close()

exit()



















