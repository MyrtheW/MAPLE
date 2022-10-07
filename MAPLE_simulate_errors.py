#import numpy as np
# import os
import random
import argparse
from sys import exit


def simulateErrors(file, errorRate=0.0005, siteSpecific=None, trim_seqlength=False, outputFile=None):
    """
    SIMULATING WITH ERROR RATES
    :param file: data location as a FASTA input
    :param errorRate: will not be used, other than for saving the output file, when siteSpecific lis is proived.
    :param siteSpecific: list of site specific errors, with lenght of lRef
    :param trim_seqlength: interval (tuple) to which the sequences have to be trimmed.
    :return: 
    """
    alphabet = {"A", "C", "T", "G"}
    fileI=open(file)
    trimmed = "_len"+ str(trim_seqlength[1]-trim_seqlength[0]) if trim_seqlength else ""
    if not outputFile:
        outputFile = file[:-3]+trimmed+"_errors"+str(errorRate)+".fa"
    fileO=open(outputFile,"w")
    line = fileI.readline()
    nSeqs = 0
    while line != "":
        while line == "\n":
            line = fileI.readline()
        if line == "":
            break
        nSeqs += 1
        seq = ""
        name = line.replace(">", "").replace("\n", "")
        line = fileI.readline()
        while line != "" and line != "\n" and line[0] != ">":
            seq += line.replace("\n", "")
            line = fileI.readline()

        seq = seq.lower()
        seq_errors = []
        for i in range(len(seq)):
            if siteSpecific:
                errorRate=siteSpecific[i]
            if seq[i] != 'N' and seq[i] != 'n' and random.random() < errorRate:  # error is True with prob = errorRate[pos],
                seq_errors.append(random.sample(population=alphabet - {seq[i]}, k=1)[0])
                #population = alphabet - {seq[i]}
                #seq_errors.append(np.random.sample(population= alphabet -{ seq[i]}, k=1)[0]) # sample from other nucleotides { }- current nucleotide. Question: sample all with equal probabilities, as would be with a homogenous model?
            else:
                seq_errors.append(seq[i])
        if trim_seqlength:
            seq_errors = seq_errors[trim_seqlength[0]:trim_seqlength[1]]
        fileO.write( "\n>" + name + "\n" + "".join(seq_errors))
    fileI.close()

    fileO.close()

def siteSpecificErrors(errorRate, lRef, seed=1):
    """
    position specific error rate: sample a distribution to find position specific errorRates with an average of errorRate.
    :param seed: seed to sample error rates
    :param lRef: length of reference
    :return: array with position specific error rates
    """
    siteSpecificErrorRates = [random.expovariate(1 / errorRate) for i in range(lRef)]
    # import numpy as np
    # np.random.seed(seed)
    # siteSpecificErrorRates = [np.random.exponential(scale=errorRate, size=None) for i in range(lRef)] #scale --> errorRate, [Kozlov]
    scaling_factor = errorRate/sum(siteSpecificErrorRates)*len(siteSpecificErrorRates)
    siteSpecificErrorRates = [item*scaling_factor for item in siteSpecificErrorRates] # making sure the average equals the errorRate
    return siteSpecificErrorRates


def getLRef(file):
    """
    :return: length of the genome in the file
    """
    fileI = open(file)
    seq, line = "", fileI.readline()
    while line == "\n":
        line = fileI.readline()
    line = fileI.readline()
    while line != "" and line != "\n" and line[0] != ">":
        seq += line.replace("\n", "")
        line = fileI.readline()
    fileI.close()
    return(len(seq))

parser = argparse.ArgumentParser(description='Estimate a tree from a diff format and using iterative approximate maximum likelihood sample placement.')
parser.add_argument('--input',default="/Users/demaio/Desktop/GISAID-hCoV-19-phylogeny-2021-03-12/phylogenetic_inference/2021-03-31_unmasked_differences_reduced.txt_consensus-based.txt", help='input MAPLE file name; should contain first the reference genome and then the difference of all samples with respet to the reference.')
parser.add_argument('--output',default=None)
parser.add_argument('--errorRate',type=float, default=0.0005)
parser.add_argument('--siteSpecific',default=False, action="store_true")
args = parser.parse_args()
file = args.input
outputFile = args.output
errorRate = args.errorRate

if args.siteSpecific:
    if errorRate ==0: # possible extensions is to allow for an input file of sitespecific errors
        raise Exception("please provide a non-zero error rate when using a site specific error rate")
    siteSpecific = siteSpecificErrors(errorRate, lRef=getLRef(file), seed=1)
    fileO=open(outputFile[:-3] + "_siteSpecificErrors.txt","w") #expecting a fasta extension to outputfile
    fileO.write( ", ".join([str(item) for item in siteSpecific]) )
    fileO.close()
else:
    siteSpecific = None
simulateErrors(file, errorRate=errorRate, siteSpecific=siteSpecific, outputFile=outputFile)
exit()
