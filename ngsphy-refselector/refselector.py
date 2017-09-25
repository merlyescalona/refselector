import argparse,datetime,logging,os,sys, msatools, sqlite3
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
	seqDescriptionFile=""
	method=0
	nsize=-1
	numLociPerReplicate=[]
	numLociPerReplicateDigits=[]
	numReplicates=0
	numReplicatesDigits=0
	db="" # Path of the SimPhy DB

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
		self.seqDescriptionFile=args.seq_desc_file
		self.nsize=args.nsize

		# Checking correctness of the given paths
		if (args.path[-1]=="/"):
			self.projectName=os.path.basename(args.path[0:-1])
		else:
			self.projectName=os.path.basename(args.path)
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
		# Generation of the output folder
		try:
			os.mkdir(self.output)
			APPLOGGER.info("Generating output folder:\t{}".format(self.output))
		except:
			APPLOGGER.info("Output folder ({0}) exists. ".format(self.output))
		if  self.nsize > -1:
			try:
				os.mkdir(os.path.join(self.output, "bed"))
				APPLOGGER.info("Generating output folder:\t{}".format(os.path.join(self.output, "bed")))
			except:
				APPLOGGER.info("Output folder ({0}) exists. ".format(os.path.join(self.output, "bed")))

	def checkArgs(self):
		APPLOGGER.info("Checking arguments...")
		APPLOGGER.info("\tSimPhy...")
		simphydir=os.path.exists(self.path)
		if simphydir:
			APPLOGGER.info("SimPhy folder exists:\t{0}".format(simphydir))
		else:
			exc_type, exc_obj, exc_tb = sys.exc_info()
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
		print(self.db)
		if not self.db in fileList:
			ex="SimPhy required file do not exist."
			exc_type, exc_obj, exc_tb = sys.exc_info()
			fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
			message="{1} | File: {2} - Line:{3}\n\t{0}".format(\
				ex,exc_type,fname, exc_tb.tb_lineno,\
				"Please verify. Exiting."\
			)
			raise NRSException(False, message, datetime.datetime.now()-self.startTime)
		APPLOGGER.info("SimPhy data base exist:\t{0} ()".format(\
			os.path.basename(self.db),self.db in fileList)
		)
		APPLOGGER.info("\tIdentifying replicates...")
		# check how many of them are dirs
		for item in fileList:
			baseitem=os.path.basename(item)
			if (os.path.isdir(os.path.abspath(item)) and  baseitem.isdigit()):
				self.numReplicates=self.numReplicates+1
		self.numReplicatesDigits=len(str(self.numReplicates))
		# check if at least one
		if not (self.numReplicates>0):
			ex="Number of replicates/folders:\t{0} [Required at least 1]".format(self.numReplicates>0)
			exc_type, exc_obj, exc_tb = sys.exc_info()
			fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
			message="{1} | File: {2} - Line:{3}\n\t{0}".format(\
				ex,exc_type,fname, exc_tb.tb_lineno,\
				"Please verify. Exiting."\
			)
			raise NRSException(False, message, datetime.datetime.now()-self.startTime)
		APPLOGGER.info("\tDone!")
		APPLOGGER.info("Number of replicates:\t{0}".format(self.numReplicates))
		if self.nsize > -1:
			APPLOGGER.info("Reference sequences will be concatenated.")
		########################################################################
		self.getNumLociPerReplicate()
		########################################################################

	def getNumLociPerReplicate(self):
		query="select N_Loci from Species_Trees"
		con = sqlite3.connect(self.db)
		res=con.execute(query).fetchall()
		con.close()
		self.numLociPerReplicate=[item for sublist in res for item in sublist]
		for item in range(0,len(self.numLociPerReplicate)):
			self.numLociPerReplicateDigits+=[len(str(self.numLociPerReplicate[item]))]

	def getBEDfile(self, repID):
		control=os.path.join(\
			self.path,\
			"{0:0{1}d}".format(repID, self.numReplicatesDigits),\
			"control.txt"\
		)
		f=open(control,"r")
		lines=f.readlines()
		f.close()
		sizes=[ int(line.replace("[","").replace("]","").strip().split()[4]) for line in lines if "PARTITION" in line]
		startpos=0
		bedfile=os.path.join(\
			self.output,\
			"bed",\
			"{0}_{1:0{2}d}.bed".format(\
				self.outputprefix,\
				repID,self.numReplicatesDigits\
			)
		)
		outfile=open(bedfile,'a')
		totalConcatSize=len(str(sum(sizes)+(self.nsize*len(sizes))))
		for locID in range(1,self.numLociPerReplicate[repID-1]+1):
			endpos=startpos+sizes[locID-1]
			chromField="{chrom:{align}0{fillChrom}d}".format(\
				align="<",\
				chrom=repID,\
				fillChrom=self.numReplicatesDigits)
			positions="{startPOS:{align}{posSIZE}}\t{endPOS:{align}{posSIZE}}".format(\
				align=">",\
				startPOS=startpos,\
				endPOS=(endpos-1),\
				posSIZE=int(str(totalConcatSize))\
			)
			nameField="{name:0{fillName}d}".format(\
				name=locID,\
				fillName=self.numLociPerReplicateDigits[repID-1]
			)
			outfile.write("Replicate{0}\t{1}\tLocus{2}\n".format(\
				chromField,\
				positions,\
				nameField\
			))
			startpos=endpos+self.nsize
		outfile.close()

	def iterateOverReplicate(self):
		APPLOGGER.debug("IterateOverReplicate")
		for repID in range(1, self.numReplicates+1):
			APPLOGGER.debug("Replicate {0}/{1} ".format(repID, self.numReplicates))
			curReplicatePath=os.path.join(\
				self.path,\
				"{0:0{1}d}".format(repID, self.numReplicatesDigits),\
			)
			prefixLoci=0
			fileList=os.listdir(curReplicatePath)
			for item in fileList:
				if ("{0}_".format(self.inputprefix) in item) and (".fasta" in item):
					prefixLoci+=1
			APPLOGGER.info("Method chosen: {0}")
			if (self.nsize>-1): # only way I know sequences are concatenated
				self.getBEDfile(repID)
			####################################################################
			if self.method==0:
				self.methodOutgroup(repID)
			####################################################################
			if self.method==1:
				if not self.seqDescriptionFile == "":
					self.seqPerLocus(repID)
				else:
					message="{0}\n\t{1}".format(\
						"Required file for method == 4.",\
						"Please verify. Exiting."
					)
					raise NRSException(False, message, datetime.datetime.now()-self.startTime)
			####################################################################
			if self.method==2:
				self.methodRandomIngroup(repID)
			####################################################################
			if self.method==3:
				self.methodConsensusRandomSpecies(repID)
			####################################################################
			if self.method==4:
				self.methodConsensusAll(repID)


	def writeLocus(self, repID,locID,description,sequence):
		APPLOGGER.debug("Write Locis()")
		if self.nsize > -1:
			APPLOGGER.debug("CONCAT")
			self.concatSelectedLoci(repID,locID,description,sequence)
		else:
			APPLOGGER.debug("SEPARATE FILE")
			self.writeSelectedLoci(repID,locID,description,sequence)

	def methodOutgroup(self,repID):
		"""
		Method 0 for the selection of reference loci.
		This method selects the outgroup as a reference locus.
		------------------------------------------------------------------------
		attributes: repID: Index of the species tree that is being used.
		returns: Nothing
		"""
		APPLOGGER.debug("method outgroup")
		description="0_0_0"
		for locID in range(1,self.numLociPerReplicate[repID-1]+1):
			APPLOGGER.info("Loci {0}/{1}".format(locID,self.numLociPerReplicate[repID-1]))
			fastapath=os.path.join(\
				self.path,\
				"{0:0{1}d}".format(repID, self.numReplicatesDigits),\
				"{0}_{1:0{2}d}.fasta".format(self.inputprefix,locID, self.numLociPerReplicateDigits[repID-1])\
			)
			lociData=msatools.parseMSAFileWithDescriptions(fastapath)
			selectedSequence=lociData[description]
			self.writeLocus(repID,locID,description,selectedSequence)
		APPLOGGER.info("Done outgroup sequence")


	def seqPerLocus(self,repID):
		"""
		Method 4 for the selection of reference loci.
		------------------------------------------------------------------------
		This method selects a sequence per locus as indicated in the seq_desc_file.
		Args: repID: Index of the species tree that is being used.
		Returns: Nothing
		"""
		for locID in range(1,self.numLociPerReplicate[repID-1]+1):
			APPLOGGER.info("Loci {0}/{1}".format(locID,self.numLociPerReplicate[repID-1]))
			lociData=self.parseLocFile(repID,locID)
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
			self.writeSelectedLociMultipleSpecies(repID,locID,sequences)
		APPLOGGER.info("Done Seq Per Species")

	def seqPerSpecies(self,repID):
		for locID in range(1,self.numLociPerReplicate[repID-1]+1):
			APPLOGGER.info("Loci {0}/{1}".format(locID,self.numLociPerReplicate[repID-1]))
			lociData=self.parseLocFile(repID,locID)
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
			self.writeSelectedLociMultipleSpecies(repID,locID,sequences)
		APPLOGGER.info("Done Seq Per Species")


	def methodRandomIngroup(self,repID):
		"""
		Method 1 for the selection of reference loci.
		------------------------------------------------------------------------
		This method selects a random sequence from the ingroup species as
		a reference and from all the loci.
		Args: repID: Index of the species tree that is being used.
		Returns: Nothing
		"""
		for locID in range(1,self.numLociPerReplicate[repID-1]+1):
			APPLOGGER.info("Loci {0}/{1}".format(locID,self.numLociPerReplicate[repID-1]))
			lociData=self.parseLocFile(repID,locID)
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
			self.writeSelectedLoci(repID,locID,selected["description"],selected["sequence"])
		APPLOGGER.info("Done random ingroup sequence")

	def methodConsensusRandomSpecies(self,repID):
		"""
		Method 2 for the selection of reference loci.
		------------------------------------------------------------------------
		This method selects a consensus sequence, obtained from the random selection
		of a species, and then computing a consensus from all the sequences
		within the selected species.
		Args: repID: Index of the species tree that is being used.
		Returns: Nothing
		"""
		for locID in range(1,self.numLociPerReplicate[repID-1]+1):
			APPLOGGER.info("Loci {0}/{1}".format(locID,self.numLociPerReplicate[repID-1]))
			lociData=self.parseLocFile(repID,locID)
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
			self.writeSelectedLoci(repID,locID,selectedDes,selected)
		APPLOGGER.info("Done random ingroup consensus")



	def methodConsensusAll(self,repID):
		"""
		Method 3 for the selection of reference loci.
		------------------------------------------------------------------------
		Computes the consensus from all the sequences of a gene tree file,
		and uses this sequence as reference loci.
		Args: repID: Index of the species tree that is being used.
		Returns: Nothing
		"""
		for locID in range(1,self.numLociPerReplicate[repID-1]+1):
			APPLOGGER.info("Loci {0}/{1}".format(locID,self.numLociPerReplicate[repID-1]))
			lociData=self.parseLocFile(repID,locID)
			keys=set(lociData.keys())-set("0")
			sequences=[]
			for mk in keys:
				subkeys=lociData[mk].keys()
				for sk in subkeys:
					sequences+=[lociData[mk][sk]["sequence"]]
			selected=self.computeConsensus(sequences)
			self.writeSelectedLoci(repID,locID,">consensus_all",selected)
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


	def writeSelectedLociMultipleSpecies(self,repID,locID,seqs):
		"""
		Writes multiple sequences in a single file.
		------------------------------------------------------------------------
		Args:
		- repID: Index of the species tree that is being used.
		- locID: Index of the locus being used.
		- seqs: sequence of the locus to be written.
		------------------------------------------------------------------------
		Returns: Nothing
		------------------------------------------------------------------------
		Generates a file for all the selected loci.
		"""
		# seqs=[(description,seq), ..., (description,seq)]
		APPLOGGER.info("Writing selected loci {1} from ST: {0}", repID,locID)
		outname=os.path.join(\
			self.output,\
			"{0}_{1:0{2}d}_{3:0{4}d}.fasta".format(\
				self.outputprefix,\
				repID,self.numReplicatesDigits,\
				locID,self.numLociPerReplicateDigits[repID-1]\
			)\
		)
		outfile=open(outname,'a')
		for item in range(0,len(seqs)):
			des=seqs[item][0]
			nucSeq=seqs[item][1]
			newDes=">{0}:{1:0{2}d}:REF:{7}:{6}:{3:0{4}d}:{5}".format(\
				self.projectName,\
				repID,\
				self.numReplicatesDigits,\
				locID,\
				self.numLociPerReplicateDigits[repID-1],\
				des[1:len(des)],\
				self.inputprefix,\
				self.outputprefix
			)
			outfile.write("{0}\n{1}\n".format(newDes,nucSeq))
		outfile.close()


	def concatSelectedLoci(self,repID,locID,description,sequence):
		"""
		BEDFILE: replicateID startPOS endPOS locID
		"""
		APPLOGGER.info("Writing selected loci {1} from ST: {0}".format(repID,locID))
		outname=os.path.join(\
			self.output,\
			"{0}_{1:0{2}d}.fasta".format(\
				self.outputprefix,\
				repID,self.numReplicatesDigits
			)
		)
		newDes=">{0}:{1}:{2:0{3}d}".format(\
			self.projectName,\
			self.outputprefix,\
			repID, self.numReplicatesDigits
		)
		nsequence="N"*self.nsize
		# I'm assuming that if the file does not exist it will be created
		if os.path.exists(outname):
			with open(outname, 'a+') as f:
				f.seek(-1,2)
				if locID == self.numLociPerReplicate[repID-1]:
					f.write('\n'.encode())
				else:
					f.write('{0}{1}'.format(sequence,nsequence).encode())
		else:
			f=open(outname,"w")
			f.write(">{0}\n{1}".format(description, sequence))


	def writeSelectedLoci(self,repID,locID,des,seq):
		"""
		Writes a single sequence per file.
		-----------------------------------------------------------------------
		Input:
			- repID: Index of the species tree that is being used.
			- locID: Index of the locus being used.
			- des: description of the sequence to be written.
			- seq: sequence of the locus to be written.
		Returns: Nothing
			- Generates a file per selected locus.
		"""
		APPLOGGER.info("Writing selected loci {1} from ST: {0}", repID,locID)
		outname=os.path.join(\
			self.output,\
			"{0}_{1:0{2}d}_{3:0{4}d}.fasta".format(\
				self.outputprefix,\
				repID,self.numReplicatesDigits,\
				locID, self.numLociPerReplicateDigits[repID-1]
			)
		)
		newDes=">{0}:{1:0{2}d}:REF:{7}:{6}:{3:0{4}d}:{5}".format(\
			self.projectName,\
			repID,\
			self.numReplicatesDigits,\
			locID,\
			self.numLociPerReplicateDigits[repID-1],\
			des[1:len(des)],\
			self.inputprefix,\
			self.outputprefix
		)
		# I'm assuming that if the file does not exist it will be created
		outfile=open(outname,'a')
		outfile.write("{0}\n{1}\n".format(newDes,seq))
		outfile.close()


	def run(self):
		"""
		Run process of the program.
		"""
		self.checkArgs()
		self.iterateOverReplicate()
		raise NRSException(True,"", datetime.datetime.now()-self.startTime )
