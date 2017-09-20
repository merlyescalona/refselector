import argparse,datetime,logging,os,sys
from MEOutputFormatter import MEOutputFormatter as mof
from MELoggingFormatter import MELoggingFormatter as mlf
from MELogLevels import MELogLevels as mll
import numpy as np
import random as rnd
from select import select

APPLOGGER=logging.getLogger("ngsphy-refselector")

class LociReferenceSelection:
	"""
		LociReferenceSelection:
		This programs selects references loci would have been used to design the probes
		for a capture experiment.

	"""
	startTime=None
	endTime=None

	method=0
	inputprefix=""
	outputprefix=""
	nsize=-1

	numLoci=0
	numLociDigits=0
	numReplicates=0
	numReplicatesDigits=0

	projectName=""
	output=""

	def __init__(self, args):
		"""
			Args:
				args: Command-line parameters.

			Returns:
				The whole processing for the selection of reference loci that would
				have been used to design the probes for a capture experiment.
		"""
		# ----------------------------------------------------------------------
		self.startTime=datetime.datetime.now()
		self.endTime=None
		APPLOGGER.info(\
			"{0}".format(\
			"RefSelector started",\
		))
		# Output folder
		self.method=args.method
		self.inputprefix=args.input_prefix
		self.outputprefix=args.output_prefix
		self.nsize=args.nsize
		self.numLoci=0
		self.numLociDigits=0
		self.numReplicates=0
		self.numReplicatesDigits=0
		self.projectName=""
		if (args.path[-1]=="/"): self.projectName=os.path.basename(args.path[0:-1])
		else: self.projectName=os.path.basename(args.path)
		self.path=os.path.abspath(args.path)

		output=os.path.abspath(args.output)
		outputFolderName=""
		if (args.output[-1]=="/"): outputFolderName=os.path.basename(args.output[0:-1])
		else: outputFolderName=os.path.basename(output)
		if (os.path.exists(output)):
			listdir=os.listdir("{}".format(os.path.dirname(output)))
			counter=0
			for item in listdir:
				if outputFolderName in item:
					counter+=1
			if not counter == 0: outputFolderName+="_{0}".format(counter+1)
		self.output=os.path.join(\
			os.path.dirname(output),\
			outputFolderName
		)
		try:
			os.mkdir(self.output)
			APPLOGGER.info("Generating output folder:\t{}".format(self.output))
		except:
			APPLOGGER.info("Output folder ({0}) exists. ".format(self.output))

	def checkArgs(self):
		"""
			Checks the correctness of the input parameters.
			Returns:
				Nothing
		"""
		APPLOGGER.info("Checking arguments...")
		simphydir=os.path.exists(self.path)
		if simphydir:
			APPLOGGER.log(mll.CONFIG_LEVEL,"SimPhy folder exists:\t{0}".format(simphydir))
		else:
			self.ending(False, "SimPhy folder does not exist.")

		fileList=os.listdir(os.path.abspath(self.path))
		for index in range(0,len(fileList)):
			fileList[index]=os.path.abspath(os.path.join(self.path,fileList[index]))

		db = os.path.join(self.path,"{0}.db".format(self.projectName))
		if not db in fileList:
			self.ending(False, "SimPhy required file do not exist.")
		APPLOGGER.info("SimPhy data base exist:\t{0} ()".format(\
			os.path.basename(db),db in fileList)
		))

		# check how many of them are dirs
		for item in fileList:
			baseitem=os.path.basename(item)
			if (os.path.isdir(os.path.abspath(item)) and  baseitem.isdigit()):
				self.numReplicates=self.numReplicates+1
		self.numReplicatesDigits=len(str(self.numReplicates))
		# check if at least one
		if not (self.numReplicates>0):
			self.ending(False, "Number of replicates/folders:\t{0} [Required at least 1]".format(self.numReplicates>0))
		APPLOGGER.info("Number of replicates:\t{0}".format(self.numReplicates))

		if self.nsize > -1:
			APPLOGGER.info("Reference sequences will be concatenated.")


	def iterateOverReplicate(self):
		APPLOGGER.debug("IterateOverReplicate")
		for repID in range(1, self.numReplicates+1):
			curReplicatePath=os.path.join(\
				self.path,\
				"{0:0{1}d}".format(repID, self.numReplicatesDigits),\
			)
			APPLOGGER.debug("Replicate {1} - Iteration {0}/{2} ".format(\
				repID, curReplicatePath,self.numReplicates))
			prefixLoci=0
			fileList=os.listdir(curReplicatePath)
			for item in fileList:
				if ("{0}_".format(self.inputprefix) in item) and (".fasta" in item):
					prefixLoci+=1

			self.numLoci=prefixLoci
			self.numLociDigits=len(str(self.numLoci))
			APPLOGGER.info("Method chosen: {0}")
			if self.method==0:
				self.methodOutgroup(repID)
			elif self.method==1:
				self.methodRandomIngroup(repID)
			elif self.method==2:
				self.methodConsensusRandomSpecies(repID)
			elif self.method==3:
				self.methodConsensusAll(repID)
			else:
				self.seqPerSpecies(repID)

	def seqPerSpecies(self,repID):
		"""
			Method 4 for the selection of reference loci.

			This method selects a random single sequence from each species
			as a reference locus.

			Args:
				repID: Index of the species tree that is being used.

			Returns:
				Nothing
		"""
		for indexLOC in range(1,self.numLoci+1):
			APPLOGGER.info("Loci {0}/{1}".format(indexLOC,self.numLoci))
			lociData=self.parseLocFile(repID,indexLOC)
			sequences=[] # [(description,seq)]
			keys=lociData.keys()
			for k in keys:
				subkeys=lociData[k].keys() # Getting subkeys
				rndKey1=0;
				try:
					rndKey1=rnd.sample(subkeys,1)[0]
				except:
					rndKey1=0

				selected=lociData[k][rndKey1]
				sequences+=[(selected["description"], selected["sequence"])]
			self.writeSelectedLociMultipleSpecies(repID,indexLOC,sequences)
		APPLOGGER.info("Done Seq Per Species")


	def methodOutgroup(self,repID):
		"""
			Method 0 for the selection of reference loci.

			This method selects the outgroup as a reference locus.

			Args:
				repID: Index of the species tree that is being used.

			Returns:
				Nothing

		"""
		for indexLOC in range(1,self.numLoci+1):
			APPLOGGER.info("Loci {0}/{1}".format(indexLOC,self.numLoci))
			lociData=self.parseLocFile(repID,indexLOC)
			selected=lociData["0"]["0"]
			self.writeSelectedLoci(repID,indexLOC,selected["description"],selected["sequence"])
		APPLOGGER.info("Done outgroup sequence")


	def methodRandomIngroup(self,repID):
		"""
			Method 1 for the selection of reference loci.

			This method selects a random sequence from the ingroup species as
			a reference and from all the loci.

			Args:
				repID: Index of the species tree that is being used.

			Returns:
				Nothing

		"""
		for indexLOC in range(1,self.numLoci+1):
			APPLOGGER.info("Loci {0}/{1}".format(indexLOC,self.numLoci))
			lociData=self.parseLocFile(repID,indexLOC)
			keys=lociData.keys()
			rndKey1="0"; rndKey2="0"
			try:
				rndKey1=rnd.sample(set(keys)-set("0"),1)[0]
			except:
				rndKey1="0"
			subkeys=lociData[rndKey1]
			try:
				rndKey2=rnd.sample(len(subkeys),1)[0]
			except:
				rndKey2="0"
			selected=lociData[rndKey1][rndKey2]
			self.writeSelectedLoci(repID,indexLOC,selected["description"],selected["sequence"])
		APPLOGGER.info("Done random ingroup sequence")

	def methodConsensusRandomSpecies(self,repID):
		"""
			Method 2 for the selection of reference loci.

			This method selects a consensus sequence, obtained from the random selection
			of a species, and then computing a consensus from all the sequences
			within the selected species.

			Args:
				repID: Index of the species tree that is being used.
			Returns:
				Nothing

		"""
		for indexLOC in range(1,self.numLoci+1):
			APPLOGGER.info("Loci {0}/{1}".format(indexLOC,self.numLoci))
			lociData=self.parseLocFile(repID,indexLOC)
			keys=lociData.keys()
			rndKey1=0; rndKey2=0
			try:
				rndKey1=rnd.sample(set(keys)-set("0"),1)[0]
			except:
				rndKey1=0
			subkeys=lociData[rndKey1]
			sequences=[]
			for sk in subkeys:
				sequences+=[lociData[rndKey1][sk]["sequence"]]
			selected=self.computeConsensus(sequences)
			selectedDes=">consensus_sp_{0}".format(rndKey1)
			self.writeSelectedLoci(repID,indexLOC,selectedDes,selected)
		APPLOGGER.info("Done random ingroup consensus")



	def methodConsensusAll(self,repID):
		"""
			Method 3 for the selection of reference loci.

			Computes the consensus from all the sequences of a gene tree file,
			and uses this sequence as reference loci.

			Args:
				repID: Index of the species tree that is being used.
			Returns:
				Nothing
		"""
		for indexLOC in range(1,self.numLoci+1):
			APPLOGGER.info("Loci {0}/{1}".format(indexLOC,self.numLoci))
			lociData=self.parseLocFile(repID,indexLOC)
			keys=set(lociData.keys())-set("0")
			sequences=[]
			for mk in keys:
				subkeys=lociData[mk].keys()
				for sk in subkeys:
					sequences+=[lociData[mk][sk]["sequence"]]
			selected=self.computeConsensus(sequences)
			self.writeSelectedLoci(repID,indexLOC,">consensus_all",selected)
		APPLOGGER.info("Done all ingroups consensus")


	def computeConsensus(self, sequences):
		"""
			Method for the computation of a consensus sequence from a set
			of sequences.

			Args:
				sequences: list of the sequences

			Returns:
				Single consensus sequence
		"""
		conseq=""
		# Assume that all the sequences in a file have the same length
		seqSize=len(sequences[0])
		# Need to know how many sequences I have
		numSeqs=len(sequences)
		# Seqs for all nucleotides
		A=np.zeros(seqSize);
		C=np.zeros(seqSize);
		G=np.zeros(seqSize);
		T=np.zeros(seqSize);
		N=np.zeros(seqSize);
		for indexCol in range(0, seqSize):
			for indexRow in range(0,numSeqs):
				if sequences[indexRow][indexCol]=="A": A[indexCol]+=1
				if sequences[indexRow][indexCol]=="C": C[indexCol]+=1
				if sequences[indexRow][indexCol]=="G": G[indexCol]+=1
				if sequences[indexRow][indexCol]=="T": T[indexCol]+=1
				if sequences[indexRow][indexCol]=="N": N[indexCol]+=1

		for indexCol in range(0,seqSize):
			if A[indexCol] > C[indexCol] and A[indexCol] > G[indexCol] and A[indexCol] > T[indexCol] and A[indexCol] > N[indexCol]:
				conseq+="A"
			elif C[indexCol] > A[indexCol] and C[indexCol] > G[indexCol] and C[indexCol] > T[indexCol] and C[indexCol] > N[indexCol]:
				conseq+="C"
			elif G[indexCol] > A[indexCol] and G[indexCol] > C[indexCol] and G[indexCol] > T[indexCol] and G[indexCol] > N[indexCol]:
				conseq+="G"
			elif T[indexCol] > A[indexCol] and T[indexCol] > C[indexCol] and T[indexCol] > G[indexCol] and T[indexCol] > N[indexCol]:
				conseq+="T"
			else: conseq+="N"

		APPLOGGER.info("Consensus computed")
		return conseq


	def writeSelectedLociMultipleSpecies(self,repID,indexLOC,seqs):
		"""
			Writes multiple sequences in a single file.

			Args:
				repID: Index of the species tree that is being used.
				indexLOC: Index of the locus being used.
				seqs: sequence of the locus to be written.

			Returns:
				Nothing

			Generates a file for all the selected loci.
		"""
		# seqs=[(description,seq), ..., (description,seq)]
		APPLOGGER.info("Writing selected loci {1} from ST: {0}", repID,indexLOC)
		outname="{0}/REF_{5}_{1:0{2}d}_{3:0{4}d}.fasta".format(\
			self.output,\
			repID,\
			self.numReplicatesDigits,\
			indexLOC,\
			self.numLociDigits,\
			self.outputprefix
		)
		outfile=open(outname,'a')
		for item in range(0,len(seqs)):
			des=seqs[item][0]
			nucSeq=seqs[item][1]
			newDes=">{0}:{1:0{2}d}:REF:{7}:{6}:{3:0{4}d}:{5}".format(\
				self.projectName,\
				repID,\
				self.numReplicatesDigits,\
				indexLOC,\
				self.numLociDigits,\
				des[1:len(des)],\
				self.inputprefix,\
				self.outputprefix
			)
			outfile.write("{0}\n{1}\n".format(newDes,nucSeq))

		outfile.close()

	def writeSelectedLoci(self,repID,indexLOC,des,seq):
		"""
			Writes a single sequence per file.

			Args:
				repID: Index of the species tree that is being used.
				indexLOC: Index of the locus being used.
				des: description of the sequence to be written.
				seq: sequence of the locus to be written.

			Returns:
				Nothing

			Generates a file per selected locus.
		"""
		APPLOGGER.info("Writing selected loci {1} from ST: {0}", repID,indexLOC)
		outname="{0}/REF_{4}_{3}_{1:0{2}d}_{5:0{6}d}.fasta".format(\
			self.output,\
			repID,\
			self.numReplicatesDigits,\
			self.inputprefix,\
			self.outputprefix,\
			indexLOC,\
			self.numLociDigits
		)
		newDes=">{0}:{1:0{2}d}:REF:{7}:{6}:{3:0{4}d}:{5}".format(\
			self.projectName,\
			repID,\
			self.numReplicatesDigits,\
			indexLOC,\
			self.numLociDigits,\
			des[1:len(des)],\
			self.inputprefix,\
			self.outputprefix
		)
		# I'm assuming that if the file does not exist it will be created
		outfile=open(outname,'a')
		outfile.write("{0}\n{1}\n".format(newDes,seq))
		outfile.close()

	def parseLocFile(self,repID,indexLOC):
		"""
			Parses the file of the corresponding locus, from the corresponding
			species tree into a dictionary.

			Args:
				repID: Index of the species tree that is being used.
				indexLOC: Index of the locus being used.
			Returns:
				A dictionary with all the sequences and their corresponding
				description.
		"""
		APPLOGGER.info("Parsing...")
		# 1. Reading the sequences of the file into a dictionary
		fastapath="{0}/{1:0{2}d}/{3}_{4:0{5}d}.fasta".format(\
			self.path,\
			repID,\
			self.numReplicatesDigits,\
			self.inputprefix,\
			indexLOC,\
			self.numLociDigits)
		fastafile=open(fastapath, 'r')
		lines=fastafile.readlines()
		fastafile.close()
		seqDict=dict()
		description=""; seq=""; tmp="";count=1
		for line in lines:
			if not (line.strip()==''):
				if (count%2==0):
					seq=line[0:-1].strip()
					try:
						test=seqDict[tmp[0]]
					except:
						seqDict[tmp[0]]={}
					try:
						seqDict[tmp[0]].update({tmp[2]:{\
    						'description':description,\
    						'sequence':seq\
							}})
					except:
						seqDict[tmp[0]][tmp[2]]={}
						seqDict[tmp[0]].update({tmp[2]:{\
							'description':description,\
							'sequence':seq\
							}})
					seq=None
					description=None
					tmp=None
				else:
					description=line[0:-1].strip()
					tmp=description[1:len(description)].split("_")
				count=count+1
		return seqDict

	def ending(self, good, message):
		"""
			Ending  (method to properly end the program)
			Args:
				good: Boolean, to represent whether is a correct or wrong ending.
				message: message to be shown.
			Returns:
				Nothing

			Exits the program.
		"""
		# good: Whether there's a good ending or not (error)
		if not good:
			APPLOGGER.error(message)
		self.endTime=datetime.datetime.now()
		APPLOGGER.info("--------------------------------------------------------------------------------")
		APPLOGGER.info("Elapsed time:\t{0}".format(self.endTime-self.startTime))
		APPLOGGER.info("Ending at:\t{0}".format(self.endTime.strftime("%a, %b %d %Y. %I:%M:%S %p")))
		sys.exit()

	def run(self):
		"""
			Run process of the program.
		"""
		self.checkArgs()
		self.iterateOverReplicate()
		self.ending(True,"Loci Reference Selection finished")
