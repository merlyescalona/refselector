import argparse,datetime,logging,os,sys, refselector, msatools
import loggingformatter as lf
import numpy as np
import random as rnd
################################################################################
# CONSTANTS
VERSION=1
MIN_VERSION=0
FIX_VERSION=1
PROGRAM_NAME="ngsphy-refselector"
AUTHOR="Merly Escalona <merlyescalona@uvigo.es>"
INSTITUTION="University of Vigo, Spain."
LOG_LEVEL_CHOICES=["DEBUG","INFO","WARNING","ERROR"]
LINE="--------------------------------------------------------------------------------"
################################################################################
# python LociReferenceSelection.py -p <prefix> -SF <simphy_path> -o outout -m method
################################################################################
# Logger init
################################################################################
ch = logging.StreamHandler()
loggerFormatter=lf.MELoggingFormatter(\
	fmt="%(asctime)s - %(levelname)s:\t%(message)s",\
	datefmt="%d/%m/%Y %I:%M:%S %p")
ch.setFormatter(loggerFormatter)
ch.setLevel(logging.NOTSET)
APPLOGGER=logging.getLogger("ngsphy-refselector")
APPLOGGER.addHandler(ch)
################################################################################
def createLogFile():
	formatString=""
	if platform.system()=="Darwin":
		formatString="%(asctime)s - %(levelname)s (%(module)s:%(lineno)d):\t%(message)s"
	else:
		formatString="%(asctime)s - %(levelname)s (%(module)s|%(funcName)s:%(lineno)d):\t%(message)s"
	fh=logging.FileHandler(\
		"{0}/{2}.{1:%Y}{1:%m}{1:%d}-{1:%H}:{1:%M}:{1:%S}.log".format(\
			os.getcwd(),\
			datetime.datetime.now(),\
			PROGRAM_NAME[0:-3].upper()\
			)\
		)
	fh.setLevel(logging.DEBUG)
	formatter=logging.Formatter(formatString)
	fh.setFormatter(formatter)
	APPLOGGER.addHandler(fh)
################################################################################
# Handling parameters
################################################################################
def handlingParameters():
	parser = argparse.ArgumentParser(\
        prog="{0} (v.{1}.{2}.{3})".format(PROGRAM_NAME,VERSION,MIN_VERSION,FIX_VERSION),\
        formatter_class=argparse.RawDescriptionHelpFormatter,\
        description=\
'''
\033[1m
================================================================================
NGSphy - RefSelector
================================================================================
\033[0m
Description:
============

For simulations of targeted-sequencing experiments under a known species/gene
tree distribution, extracts the reference sequences that would have been used as
target in the probe desing.

Assumptions:
============

- We are working under a SimPhy - NGSphy simulation pipeline scenario.
Following the same hierarchical folder structure.

- Also, it is assumed that the SimPhy folder project has been compressed using simphycompress(https://github.com/merlyescalona/simphycompress) and the
length of the concatenation N sequence used is known. To know more about the simulation pipeline scenario go to:
    - SimPhy's repository (https://github.com/adamallo/simphy), and/or check:
        - Mallo D, de Oliveira Martins L, Posada D (2016) SimPhy: Phylogenomic Simulation of Gene, Locus and Species Trees. *Syst. Biol.* **65**(2) 334-344. doi: (http://dx.doi.org/10.1093/sysbio/syv082
    - INDELible's site (http://abacus.gene.ucl.ac.uk/software/indelible/), and/or check:
        - Fletcher, W. and Yang, Z.( 2009) INDELible: a flexible simulator of biological sequence evolution. *Mol. Biol. Evol.* **26**(8):1879-1888. doi: https://doi.org/10.1093/molbev/msp098
    - NGSphy'S repository (https://github.com/merlyescalona/ngsphy), and/or check:
        - Escalona M, Rocha S, Posada D (2018) NGSphy: phylogenomic simulation of next-generation sequencing data. *Bioinformatics.* bty146. doi: https://doi.org/10.1093/bioinformatics/bty146
    - simphycompress repository (https://github.com/merlyescalona/simphycompress)

- Species tree replicates are filtered based on the number of loci (number of sequences/ FASTA files) existing
in each folder.

- True sequences from the SimPhy/INDELible simulation process do not contain N's.

General:
========
Specified method to obtain the reference loci used for the design of probes.
Values range from 0-4 ( Default: (0)), where:

	(0): 	Considers the outgroup sequence as the reference loci.
	(1): 	Extracts a specific sequence per locus
	(2): 	Selects a random sequence from the ingroups. Same sequence throughout the loci.
	(3): 	Selects randomly a specie and generates a consensus sequence of the
		sequences belonging to that species.
	(4): 	Generates a consensus sequences from all the sequences involved.
		(will need parameter -sdf/--seq-desc-file)
NOTE: 	The higher the method number, the longer it will take to generate the
	reference loci.

Output:
=======

- 	The output will be a directory of FASTA files.
- 	There should be as many FASTA files as replicates have been generated for
	the current SimPhy project.
- 	Each file will contain all the selected loci.
		''',\
			epilog="Version {0}.{1}.{2} (Still under development)".format(VERSION,MIN_VERSION,FIX_VERSION),\
			add_help=False
		)
	requiredGroup= parser.add_argument_group('Required arguments')
	requiredGroup.add_argument('-s','--simphy-path',metavar='<path>', type=str,\
		help='Path of the SimPhy folder.', required=True)
	requiredGroup.add_argument('-ip','--input-prefix', metavar='<input_prefix>', type=str,\
		help='Prefix of the FASTA filenames.', required=True)
	requiredGroup.add_argument('-op','--output-prefix', metavar='<output_prefix>', type=str,\
		help='Prefix for the output filename.', required=True)
	requiredGroup.add_argument('-m','--method', metavar="<method_code>",type=int,\
		choices=range(0,6), default=0,\
		help="Specified method to obtain the reference loci used for the design of probes. Values range from 0-4.",\
        required=True)
	requiredGroup.add_argument('-o','--output', metavar="<output_path>",type=str,\
		help="Path where output will be written. ",\
		required=True)
	optionalGroup= parser.add_argument_group('Optional arguments')
	optionalGroup.add_argument('-sdf','--seq-desc-file',metavar='<sequence_descriptions_file_path>', type=str,\
		help="When method = 4 has been selected, it is required to identify a file with the sequence descriptions, "+\
		"used as identifiers, corresponding to the sequence that will be used as reference per replicate.")
	optionalGroup.add_argument('-n','--nsize',metavar='<N_seq_size>', type=int,\
		default=-1,\
		help="Number of N's that will be introduced to separate the reference sequences selected. "+\
            "If the parameter is not set, the output file per replicate will be a multiple alignment sequence file, "+\
            "otherwise, the output will be a single sequence file per replicate consisting of a concatenation  "+\
            "of the reference sequences selected separated with as many N's as set for this parameter.")
	optionalGroup.add_argument('-l','--log',metavar='<log_level>', type=str,\
		choices=LOG_LEVEL_CHOICES, default="INFO",\
		help='Specified level of log that will be shown through the standard output. Entire log will be stored in a separate file. Values:{0}. Default: {1}. '.format(LOG_LEVEL_CHOICES,LOG_LEVEL_CHOICES[1]))
	informationGroup= parser.add_argument_group('Information arguments')
	informationGroup.add_argument('-v', '--version',\
		action='version',\
		version='Version {0}.{1}.{2}'.format(VERSION,MIN_VERSION,FIX_VERSION),\
		help="Show program's version number and exit")
	informationGroup.add_argument('-h', '--help',\
		action='store_true',\
		help="Show this help message and exit")
	try:
		tmpArgs = parser.parse_args()
	except:
		sys.stdout.write("\n\033[1m{}\033[0m\n".format(LINE))
		APPLOGGER.error("Something happened while parsing the arguments.")
		APPLOGGER.error("Please verify. Exiting.\n{}".format(LINE))

		parser.print_help()
		sys.exit(-1)
	return tmpArgs

################################################################################
# MAIN
################################################################################
def main():
	try:
		cmdArgs = handlingParameters()
		print(cmdArgs)
		APPLOGGER.setLevel(cmdArgs.log.upper())
		APPLOGGER.debug("Args. introduced: {}".format(cmdArgs))
		prog = refselector.ReferenceSelection(cmdArgs)
		prog.run()
	except refselector.NRSException as ex:
	    if ex.expression:
	        APPLOGGER.info("REFSELECTOR finished properly.")
	        APPLOGGER.info("Elapsed time (ETA):\t{0}".format(ex.time))
	        APPLOGGER.info("Ending at:\t{0}".format(datetime.datetime.now().strftime("%a, %b %d %Y. %I:%M:%S %p")))
	        sys.exit()
	    else:
	        APPLOGGER.error(ex.message)
	        APPLOGGER.error("Elapsed time (ETA):\t{0}".format(ex.time))
	        APPLOGGER.error("Ending at:\t{0}".format(datetime.datetime.now().strftime("%a, %b %d %Y. %I:%M:%S %p")))
	        sys.exit(-1)
	except KeyboardInterrupt:
	    APPLOGGER.error("{0}{1}\nProgram has been interrupted.{2}\nPlease run again for the expected outcome.\n{3}\n".format("\033[91m","\033[1m","\033[0m",LINE))
	    sys.exit(-1)

if __name__=="__main__":
	main()
