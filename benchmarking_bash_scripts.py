nRepeats = 10
nRepeatsSimu = 10
folderName = "/nfs/research/goldman/demaio/fastLK/realData/subsamples/"
folderNameSimu = "/nfs/research/goldman/demaio/fastLK/simulations/subsamples/"  # folderNameSimu = args.pathToSimulationFolder + folder + '/' + str(j) + "subsamples/output_repl\"$i\"/"
speeds = ["slowest", "slow", "medium", "fast", "fastest"]
allowedFails = [5, 5, 5, 4, 3]
thresholdLogLKs = [120.0, 100.0, 80.0, 60.0, 40.0]
bLenAdjustments = [10, 10, 10, 1, 1]
numTopologyImprovements = [5, 3, 2, 1, 0]
allowedFailsTopology = [6, 4, 3, 2, 1]
thresholdLogLKtopology = [150.0, 100.0, 80.0, 60.0, 40.0]
thresholdTopologyPlacement = [-0.1, -0.2, -0.5, -1.0, -2.0]
bLenFactors = [4.0, 4.0, 3.0, 2.0, 1.0]
sampleSizes = [100, 1000, 2000, 5000]  # , 10000, 20000, 50000, 100000, 200000, 500000]
errorRates = [0, 0.0001, 0.0005]
readFileUser = "demaio"
writeFileUser = "myrthe"
folderNameCode = "/nfs/research/goldman/" + writeFileUser + "/fastLK/code/"  # folderNameSimu = args.pathToSimulationFolder + folder + '/' + str(j) + "subsamples/output_repl\"$i\"/"
siteSpecifics = [True, False]
console_output = folderNameCode + "results/console_output.txt"
console_errors = folderNameCode + "results/console_errors.txt"

"""SIMULATE ERRORS""" #also site specific ones.
file = open("simulateErrors.sh","w")
# file.write("mkdir " + folderNameCode + "\n")
file.write("rm -r " + console_output + "\n")  # + "; mkdir " + console_output + "\n")
file.write("rm -r " + console_errors + "\n")  # +"; mkdir " + console_errors + "\n")
file.write(" \n#SIMULATE ERRORS \n")
for errorRate in errorRates:
	for siteSpecific in siteSpecifics:
		if siteSpecific and not errorRate:
			continue
		for j in sampleSizes:
			pathRead = "/nfs/research/goldman/" + readFileUser + "/fastLK/simulations/subsamples/" + str(j) + 'subsamples/'  # +"subsamples/"
			pathWrite = "/nfs/research/goldman/" + writeFileUser + "/fastLK/simulations/subsamples/" + str(j) + 'subsamples/'  # +"subsamples/"
			fileNameIn = "fastaFile_repeat\"$i\"_" + str(j) + "samples_Ns.fa"
			fileNameOut = "fastaFile_repeat\"$i\"_" + str(j) + "samples_Ns_" + ("sitespecific_" if siteSpecific else '')  + "errors" + str(errorRate) + ".fa"
			file.write("for i in $(seq 1 " + str(nRepeatsSimu) + ")\n" + "do \n\t"
			           + "mkdir -p " + pathWrite + "\n\t"  # perhaps make sure we only do so if it doesn't exist already.
			           + "bsub -M " + str(int(400 + j / 2)) + " -o " + console_output + " -e " + console_errors  # job creation ; output file and error file. ,
			           + " /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 "
			           + folderNameCode + "MAPLE_simulate_errors.py "  # when using python instead?
			           + " --input " + pathRead + fileNameIn
			           + " --output " + pathWrite + fileNameOut
			           + (" --siteSpecific " if siteSpecific else '') # afterwards an error rate file is stored. outputfile[:-3] + "_siteSpecificErrors.txt","w")
			           + " --errorRate " + str(errorRate) + "\n done\n\n")

"""CALL MAPLE FILE"""
file.write(" \n#CREATE MAPLE FILES \n")
for errorRate in errorRates:
	for siteSpecific in siteSpecifics:
		if siteSpecific and not errorRate:  # No site specific errors can obviously be simulated for an error rate of 0
			continue
		for j in sampleSizes:
			pathWrite = "/nfs/research/goldman/" + writeFileUser + "/fastLK/simulations/subsamples/" + str(j) + 'subsamples/'
			fileNameIn = "fastaFile_repeat\"$i\"_" + str(j) + "samples_Ns_" + ("sitespecific_" if siteSpecific else '') + "errors" + str(errorRate) + ".fa"
			fileNameOut = "diffFile_repeat\"$i\"_" + str(j) + "samples_Ns_" + ("sitespecific_" if siteSpecific else '') + "errors" + str(errorRate) + ".txt"
			file.write("for i in $(seq 1 " + str(nRepeatsSimu) + ")\n" + "do \n\t"
			           + "bsub -M " + str(int(400 + j / 2)) + " -o " + console_output
			           + " -e " + console_errors
			           + " /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 "
			           + "/nfs/research/goldman/demaio/fastLK/code/createMapleFile.py "
			           + " --path " + pathWrite
			           + " --fasta " + fileNameIn
			           + " --output " + fileNameOut
			           + " --overwrite\n"
			           + "done\n\n")
file.close()

"""Check if files exist"""
j = 5000
errorRate = 0.0001
siteSpecific = True
filename = "diffFile_repeat\"$i\"_" + str(j) + "samples_Ns_errors" + str(errorRate) + ".txt"
filename = "fastaFile_repeat\"$i\"_" + str(j) + "samples_Ns_" + ("sitespecific_" if siteSpecific else '') + "errors" + str(errorRate) +  "_siteSpecificErrors.txt"

pathWrite = "/nfs/research/goldman/" + writeFileUser + "/fastLK/simulations/subsamples/" + str(j) + 'subsamples/'  # +"subsamples/"

print("for i in $(seq 1 10)\n" + "do \n\t",
"FILE=" + pathWrite + filename,
      """
if [ -f "$FILE" ]; then
  echo "$FILE exists."
  FILESIZE=$(stat -c%s "$FILE")
  MINSIZE=100
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

"""RUN MAPLE """ #also for sitespecific erros
version = "0.1.9_error_site_specific"
file = open("submitMAPLEv" + version + ".sh", "w")
file.write("rm -r " + console_output + "\n")  # + "; mkdir " + console_output + "\n")
file.write("rm -r " + console_errors + "\n")  # +"; mkdir " + console_errors + "\n")
file.write(" \n # creating bash scripts for MAPLE \n")
for siteSpecific in siteSpecifics:  # normal maple (error rate =0, in total 4 extra runs)
	for simulationErrorRate in errorRates:
		if not(siteSpecific and not simulationErrorRate): # No site specific errors can possibly have been simulated for an error rate of 0.
			for inferenceErrorRate in errorRates:
				if not(siteSpecific and (inferenceErrorRate != 0 and inferenceErrorRate != simulationErrorRate)): # At the moment, I only want to run site specific error inference for site specific error data
					for j in sampleSizes:
						pathRead = "/nfs/research/goldman/" + readFileUser + "/fastLK/simulations/subsamples/" + str(j) + 'subsamples/'  # +"subsamples/"
						pathWrite = "/nfs/research/goldman/" + writeFileUser + "/fastLK/simulations/subsamples/" + str(j) + 'subsamples/'
						fileName = "diffFile_repeat\"$i\"_" + str(j) + "samples_Ns_" + ("sitespecific_" if siteSpecific else '') + "errors" + str(simulationErrorRate) + ".txt"  # if errorRate else "diffFile_repeat\"$i\"_"+str(j)+"samples_Ns"  + ".txt"
						treeFile = "treeFile_repeat\"$i\"_" + str(j) + "samples.nw"
						siteSpecificFile = "fastaFile_repeat\"$i\"_" + str(j) + "samples_Ns_" + ("sitespecific_" if siteSpecific else '') + "errors" + str(errorRate) +  "_siteSpecificErrors.txt"
						benchmarkingFile = "benchmarkingFile4.tsv"	# possible to do remove existing result if already exists - this helps in case of re-running an analysis, and in case the new version fails, to wrongly use the old tree as estimate of the new version.
						file.write("for i in $(seq 1 10)\n do\n\t") # run MAPLE for various repeats.
						file.write("bsub -M " + str(int(400 + j / 20))
						           + " -o " + console_output
						           + " -e " + console_errors
						           + " /hps/software/users/goldman/pypy3/pypy3.7-v7.3.5-linux64/bin/pypy3 " + folderNameCode + "MAPLEv" + version + ".py "
						           + " --input " + pathWrite + fileName
						           + (" --errorRate " + str(inferenceErrorRate) if not siteSpecific else "")#when siteSpecific errorRates are used, this does not matter.
						           + (" --errorRateSiteSpecific " + pathWrite + siteSpecificFile if (siteSpecific and inferenceErrorRate) else "")
						           + " --benchmarkingFile " + folderNameCode + "results/" + benchmarkingFile
						           + " --trueTree " + pathRead + treeFile
						           + " --calculateLKfinalTree "
						           + "--overwrite"
						           + " --output " + pathWrite + "outputFile_infError_" + str(inferenceErrorRate) + fileName[8:]
						           + "\n")
						file.write("done\n\n")
file.close()
print("Created MAPLE bash script submitMAPLEv" + version + ".sh")




