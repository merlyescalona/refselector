import datetime,itertools,logging,os,sys, msatools, sqlite3,re, glob, gzip
import numpy as np
import random as rnd

APPLOGGER=logging.getLogger("ngsphy-refselector")
class NRSException(Exception):
	def __init__(self, expression, message, time):
		self.expression = expression
		self.message = message
		self.time= time

class ReferenceSelection:
	startTime=None
	endTime=None
	path=""
	projectName=""
	output=""
	inputprefix=""
	outputprefix=""
	method=0
	nsize=0
	numLociPerReplicate=[]
	numLociPerReplicateDigits=[]
	numReplicates=0
	numReplicatesDigits=0
	referenceIndexFile=""
	referenceLabels=[]
	offtargetloci=[]
	offtargetlocifile=""

	def __init__(self, args):
		self.startTime=datetime.datetime.now()
		self.endTime=None
		APPLOGGER.info(\
			"{0}".format(\
			"RefSelector started",\
		))
		# Variable initialization
		self.method=args.method
		self.inputprefix=args.input_prefix
		self.outputprefix=args.output_prefix
		self.nsize=args.nsize
		if args.seq_desc_file:
			self.referenceIndexFile=os.path.abspath(args.seq_desc_file)

		########################################################################
		# Checking correctness of the given paths
		if (args.simphy_path[-1]=="/"):
			self.projectName=os.path.basename(args.simphy_path[0:-1])
		else:
			self.projectName=os.path.basename(args.simphy_path)
		self.path=os.path.abspath(args.simphy_path)
		output=os.path.abspath(args.output)
		outputFolderName=""
		if (args.output[-1]=="/"):
			outputFolderName=os.path.basename(args.output[0:-1])
		else:
			outputFolderName=os.path.basename(output)
		if (os.path.exists(output)):
			listdir=os.listdir("{}".format(os.path.dirname(output)))
			counter=0
			for item in listdir:
				if outputFolderName in item:
					counter+=1
			if not counter == 0: outputFolderName+="_{0}".format(counter+1)
		self.output=os.path.join(os.path.dirname(output),outputFolderName)
		if (args.off_target_loci_file):
			if (os.path.exists(os.path.abspath(args.off_target_loci_file))):
				self.offtargetlocifile=os.path.abspath(args.off_target_loci_file)
				self.readOffTargetLociFile()
				APPLOGGER.info("Using off-target loci")
		else:
			APPLOGGER.info("NOT using off-target loci")
			
		########################################################################
		# Generation of the output folder
		try:
			os.mkdir(self.output)
			APPLOGGER.info("Generating output folder:\t{}".format(self.output))
		except:
			APPLOGGER.info("Output folder ({0}) exists. ".format(self.output))

	def checkArgs(self):
		APPLOGGER.info("Checking arguments...")
		APPLOGGER.info("\tSimPhy...")
		simphydir=os.path.exists(self.path)
		########################################################################
		if simphydir:
			APPLOGGER.info("SimPhy folder exists:\t{0}".format(simphydir))
		else:
			exc_type, _, exc_tb = sys.exc_info()
			fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
			ex="SimPhy folder does not exist."
			message="{1} | File: {2} - Line:{3}\n\t{0}".format(\
				ex,exc_type,fname, exc_tb.tb_lineno,\
				"Please verify. Exiting.")
			raise NRSException(False, message, datetime.datetime.now()-self.startTime)
		fileList=os.listdir(os.path.abspath(self.path))
		for index in range(0,len(fileList)):
			fileList[index]=os.path.abspath(os.path.join(self.path,fileList[index]))
		self.db = os.path.join(self.path,"{0}.db".format(self.projectName))
		if not self.db in fileList:
			ex="SimPhy required file do not exist."
			exc_type, _, exc_tb = sys.exc_info()
			fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
			message="{1} | File: {2} - Line:{3}\n\t{0}".format(\
				ex,exc_type,fname, exc_tb.tb_lineno,\
				"Please verify. Exiting."\
			)
			raise NRSException(False, message, datetime.datetime.now()-self.startTime)
		APPLOGGER.info("SimPhy data base exist:\t{0} ({1})".format(\
			os.path.basename(self.db),self.db in fileList)
		)
		APPLOGGER.info("\tIdentifying replicates...")
		# check how many of them are dirs
		for item in fileList:
			baseitem=os.path.basename(item)
			if (os.path.isdir(os.path.abspath(item)) and  baseitem.isdigit()):
				self.numReplicates=self.numReplicates+1
		self.numReplicatesDigits=len(str(self.numReplicates))
		########################################################################
		# check if at least one
		if not (self.numReplicates>0):
			ex="Number of replicates/folders:\t{0} [Required at least 1]".format(self.numReplicates>0)
			exc_type, _, exc_tb = sys.exc_info()
			fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
			message="{1} | File: {2} - Line:{3}\n\t{0}".format(\
				ex,exc_type,fname, exc_tb.tb_lineno,\
				"Please verify. Exiting."\
			)
			raise NRSException(False, message, datetime.datetime.now()-self.startTime)
		APPLOGGER.info("\tDone!")
		APPLOGGER.info("Number of replicates:\t{0}".format(self.numReplicates))
		########################################################################
		self.referenceLabels=["1_0_0" for item in range(0,self.numReplicates)]
		self.numLociPerReplicate=[0]*self.numReplicates
		self.numLociPerReplicateDigits=[0]*self.numReplicates
		self.getNumLociPerReplicate()
		APPLOGGER.info("Number of loci per replicates: {}".format(self.numLociPerReplicate))
		APPLOGGER.info("Number of loci digits per replicates: {}".format(self.numLociPerReplicateDigits))
		########################################################################
		if self.nsize > -1:
			APPLOGGER.info("Reference sequences will be concatenated.")
		########################################################################
		if self.method==1:
			if not (os.path.exists(self.referenceIndexFile) and os.path.isfile(self.referenceIndexFile)):
				message="{0}\n\t{1}".format(\
					"Sequence description file does not exist.",\
					"Please verify. Exiting."\
				)
				raise NRSException(False, message, datetime.datetime.now()-self.startTime)
			else:
				self.readReferenceIndices()
		########################################################################

	def readReferenceIndices(self):
		APPLOGGER.info("Reading the reference index file")
		with open(self.referenceIndexFile, "rb") as f:
			referenceLabels=[line.strip() for line in f if not line.strip() ==""]
		for item in range(0,len(referenceLabels)):
			self.referenceLabels[item]=referenceLabels[item]


	def filterSpeciesTreeReplicates(self):
		self.filtered=[(idx+1) for idx in range(0,self.numReplicates) \
		 	if not self.numLociPerReplicate[idx]==0]

	def getNumLociPerReplicate(self):
		for index in range(0, self.numReplicates):
			repID=index+1
			APPLOGGER.debug("Replicate {0}/{1} ".format(repID, self.numReplicates))
			filename=os.path.join(\
				self.path,\
				"{0:0{1}d}".format(repID, self.numReplicatesDigits),\
				"{0}_{1:0{2}d}_TRUE.fasta.gz".format(self.inputprefix,repID, self.numReplicatesDigits)
			)
			num=0
			if os.path.exists(filename):
				seqDict=self.readGzippedFastaFile(filename)
				key=seqDict.keys()[0]
				seqs=seqDict[key]
				seqs=seqs.split("N")
				seqs=[item for item in seqs if not item == ""]
				num=len(seqs)
			self.numLociPerReplicate[index]=num
			self.numLociPerReplicateDigits[index]=len(str(num))

	def readGzippedFastaFile(self, filename):
		with gzip.open(filename, 'rb') as f:
		    file_content = f.read()
		file_content=file_content.split("\n")
		file_content=[line for line in file_content if not line ==""]
		descriptionIndices=range(0,len(file_content),2)
		seqDict=dict()
		for item in descriptionIndices:
		    key=file_content[item][1:]
		    seqDict[key]=file_content[item+1]
		return seqDict
	
	def readOffTargetLociFile(self):
		APPLOGGER.info("Reading OFF-target loci file")
		with open(self.offtargetlocifile, "rb") as f: 
			self.offtargetloci=f.read().strip().split(",")			
		self.offtargetloci=[ int(item) for item in self.offtargetloci if item != ""]
		self.offtargetloci.sort()
		

	def iterateOverReplicate(self):
		APPLOGGER.debug("IterateOverReplicate")
		for index in range(0, self.numReplicates):
			repID=index+1
			if repID in self.filtered:
				APPLOGGER.debug("Replicate {0}/{1} ".format(repID, self.numReplicates))
				APPLOGGER.info("Method chosen: {0}".format(self.method))
				####################################################################
				if self.method==0:
					self.methodOutgroup(index)
				####################################################################
				if self.method==1:
					self.methodSelectFromFile(index)
				####################################################################
				if self.method==2:
					self.methodRandomIngroup(index)
				####################################################################
				if self.method==3:
					self.methodConsensusRandomSpecies(index)
				####################################################################
				if self.method==4:
					self.methodConsensusAll(index)
				####################################################################

	def writeSequence(self,index, description, sequence):
		APPLOGGER.info("Writing reference sequence for replicate {0:0{1}d}".format(index+1,self.numReplicatesDigits))
		filename=os.path.join(\
			self.output,\
			"{0}_{1:0{2}d}.fasta".format(\
				self.outputprefix,
				index+1,\
				self.numReplicatesDigits
			)
		)
		nsequence="".join(["N"]*self.nsize)
		splitSequence=sequence.split("N")
		locsequences=[splitSequence[s] for s in range(0,len(splitSequence)) if not splitSequence[s] == "" and not s in self.offtargetloci]
		seq=["{}{}".format(locsequences[s],nsequence) for s in range(0,len(locsequences))]
		seq="{}{}".format("".join(seq),locsequences[-1])
		with open(filename,"wb") as f:
			f.write(">{}\n{}\n".format(description,seq))
		locsequences=[splitSequence[s] for s in range(0,len(splitSequence)) if not splitSequence[s] == ""]
		seq=["{}{}".format(locsequences[s],nsequence) for s in range(0,len(locsequences))]
		seq="{}{}".format("".join(seq),locsequences[-1])
		self.generateBEDFile(index, description, seq)

	def methodSelectFromFile(self,index):
		APPLOGGER.debug("method rndingroup")
		repID=index+1
		fastapath=os.path.join(\
			self.path,\
			"{0:0{1}d}".format(repID, self.numReplicatesDigits),\
			"{0}_{1:0{2}d}_TRUE.fasta.gz".format(self.inputprefix,(index+1), self.numReplicatesDigits)\
		)
		seqDict=self.readGzippedFastaFile(fastapath)
		description=self.referenceLabels[index]
		selectedSequence=seqDict[description]
		self.writeSequence(index,description,selectedSequence)
		APPLOGGER.info("Done selected sequence")


	def methodRandomIngroup(self,index):
		"""
		Method 0 for the selection of reference loci.
		This method selects the outgroup as a reference locus.
		"""
		APPLOGGER.debug("method rndingroup")
		repID=index+1
		description0="0_0_0"
		fastapath=os.path.join(\
			self.path,\
			"{0:0{1}d}".format(repID, self.numReplicatesDigits),\
			"{0}_{1:0{2}d}_TRUE.fasta.gz".format(self.inputprefix,(index+1), self.numReplicatesDigits)\
		)
		seqDict=self.readGzippedFastaFile(fastapath)
		description="1_0_0"
		keys=seqDict.keys()
		try:
			k=rnd.choice(keys)[0]
		except:
			pass
		description=keys[int(k)]
		if description==description0:
			description="1_0_0"
		selectedSequence=seqDict[description]
		self.writeSequence(index,description,selectedSequence)
		APPLOGGER.info("Done random ingroup sequence")

	def methodOutgroup(self,index):
		"""
		Method 0 for the selection of reference loci.
		This method selects the outgroup as a reference locus.
		------------------------------------------------------------------------
		attributes: repID: Index of the species tree that is being used.
		returns: Nothing
		"""
		repID=index+1
		APPLOGGER.debug("method outgroup")
		description="0_0_0"
		fastapath=os.path.join(\
			self.path,\
			"{0:0{1}d}".format(repID, self.numReplicatesDigits),\
			"{0}_{1:0{2}d}_TRUE.fasta.gz".format(self.inputprefix,(index+1), self.numReplicatesDigits)\
		)
		seqDict=self.readGzippedFastaFile(fastapath)
		selectedSequence=seqDict[description]
		self.writeSequence(index,description,selectedSequence)
		APPLOGGER.info("Done outgroup sequence")

	def methodConsensusRandomSpecies(self,index):
		"""
		Method 2 for the selection of reference loci.
		------------------------------------------------------------------------
		This method selects a consensus sequence, obtained from the random selection
		of a species, and then computing a consensus from all the sequences
		within the selected species.
		Args: repID: Index of the species tree that is being used.
		Returns: Nothing
		"""
		repID=index+1
		fastapath=os.path.join(\
			self.path,\
			"{0:0{1}d}".format(repID, self.numReplicatesDigits),\
			"{0}_{1:0{2}d}_TRUE.fasta.gz".format(self.inputprefix,(index+1), self.numReplicatesDigits)\
		)
		seqDict=self.readGzippedFastaFile(fastapath)
		speciesKeys=np.unique([item.split("_")[0] for item in seqDict.keys()])
		rndKey1=1
		try:
			rndKey1=rnd.sample(set(speciesKeys),1)[0]
		except:
			rndKey1=1
		if rndKey1==0: rndKey1=1
		subkeys=[item for item in seqDict.keys() if item.split("_")[0]==rndKey1]
		sequences=[seqDict[item] for item in seqDict.keys() if item in subkeys]
		selected=self.computeConsensus(sequences)
		selectedDes="consensus_sp_{0}".format(rndKey1)
		self.writeSequence(index,selectedDes,selected)
		APPLOGGER.info("Done random ingroup consensus")

	def methodConsensusAll(self,index):
		"""
		Method 3 for the selection of reference loci.
		------------------------------------------------------------------------
		Computes the consensus from all the sequences of a gene tree file,
		and uses this sequence as reference loci.
		Args: repID: Index of the species tree that is being used.
		Returns: Nothing
		"""
		repID=index+1
		fastapath=os.path.join(\
			self.path,\
			"{0:0{1}d}".format(repID, self.numReplicatesDigits),\
			"{0}_{1:0{2}d}_TRUE.fasta.gz".format(self.inputprefix,(index+1), self.numReplicatesDigits)\
		)
		seqDict=self.readGzippedFastaFile(fastapath)
		sequences=[seqDict[k] for k in seqDict.keys()]
		selected=self.computeConsensus(sequences)
		self.writeSequence(index,"consensus_all",selected)
		APPLOGGER.info("Done all ingroups consensus")


	def computeConsensus(self, sequences):
		"""
		Method for the computation of a consensus sequence from a set
		of sequences.
		------------------------------------------------------------------------
		Args: sequences: list of the sequences
		Returns: Single consensus sequence
		"""
		conseq=""
		# Assume that all the sequences in a file have the same length
		seqSize=len(sequences[0])
		# Need to know how many sequences I have
		numSeqs=len(sequences)
		# Seqs for all nucleotides
		A=np.zeros(seqSize)
		C=np.zeros(seqSize)
		G=np.zeros(seqSize)
		T=np.zeros(seqSize)
		N=np.zeros(seqSize)
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
			else: conseq+=rnd.choice(["A","C","G","T"])

		APPLOGGER.info("Consensus computed")
		return conseq

	def generateBEDFile(self,index,description, sequence):
		APPLOGGER.info("Writing BED file for replicate {0:0{1}d}".format(index+1,self.numReplicatesDigits))
		splitSequence=sequence.split("N")
		seqlist=[splitSequence[s] for s in range(0,len(splitSequence)) if not splitSequence[s] == "" and not s in self.offtargetloci]
		sizelist=[len(item) for item in seqlist]
		repID=index+1
		bed=dict()
		startpos=0; endpos=0
		for i in range(0, len(seqlist)):
			if not (i+1) in self.offtargetloci:
				endpos=startpos+sizelist[i]
				locID="LOCUS_{0:0{1}d}".format(\
					i+1,\
					self.numLociPerReplicateDigits[index]\
				)
				bed[locID]={"start":startpos,"end":endpos}
				startpos=endpos+(self.nsize-1)

		bedfile=os.path.join(\
			self.output,\
			"{0}_{1:0{2}d}.bed".format(\
				self.outputprefix,\
				repID,self.numReplicatesDigits\
			)
		)
		totalConcatSizeDigits=len(str(len(sequence)))
		with open(bedfile,'a') as outfile:
			for i in range(0, len(seqlist)):
				if not (i+1) in self.offtargetloci:
					locID="LOCUS_{0:0{1}d}".format(i+1,self.numLociPerReplicateDigits[index])
					positions="{startPOS:{align}{posSIZE}}\t{endPOS:{align}{posSIZE}}".format(\
						align=">",\
						startPOS=bed[locID]["start"],\
						endPOS=bed[locID]["end"],\
						posSIZE=totalConcatSizeDigits\
					)
					outfile.write("{0}\t{1}\t{2}\n".format(\
						locID,\
						positions,\
						description\
					))

	def run(self):
		"""
		Run process of the program.
		"""
		self.checkArgs()
		self.filterSpeciesTreeReplicates()
		self.iterateOverReplicate()
		raise NRSException(True,"",datetime.datetime.now()-self.startTime)
