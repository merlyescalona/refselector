import argparse,datetime,logging,os,sys
import LociReferenceSelection as lrf
from MEOutputFormatter import MEOutputFormatter as mof
from MELoggingFormatter import MELoggingFormatter as mlf
from MELogLevels import MELogLevels as mll
import numpy as np
import random as rnd
from select import select
################################################################################
# CONSTANTS
VERSION=1
MIN_VERSION=0
FIX_VERSION=0
PROGRAM_NAME="ngsphy-refselector"
AUTHOR="Merly Escalona <merlyescalona@uvigo.es>"
INSTITUTION="University of Vigo, Spain."
LOG_LEVEL_CHOICES=["DEBUG","INFO","WARNING","ERROR"]
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
def handlingParameters():
	"""
		Handling parameters.
		This functiuon takes care of the cmd-line input.
	"""

	parser = argparse.ArgumentParser(\
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
- We are workning under a SimPhy - NGSphy simulation pipeline scenario.
Meaning, it follows hierarchical SimPhy's folder structure and sequence
labeling.

General:
========
Specified method to obtain the reference loci used for the design of probes.
Values range from 0-4 ( Default: (0)), where:

	(0): Considers the outgroup sequence as the reference loci.
    (1): Extracts a specific sequence per replicate (will need parameter -sdf/--seq-desc-file)
	(2): Selects a random sequence from the ingroups.
	(3): Selects randomly a specie and a consensus sequence of the sequences belonging to that species.
	(4): Generates a consensus sequences from all the sequences involved.

NOTE: The higher the method number, the longer it will take to generate the reference loci.

Output:
=======

- The output will be a directory of FASTA files.
- There should be as many FASTA files as replicates have been generated for the current SimPhy project.
- Each file will contain all the selected loci.

		''',\
			epilog="Version {0}.{1}.{2} (Still under development)".format(VERSION,MIN_VERSION,FIX_VERSION),\
			add_help=False
		)
	requiredGroup= parser.add_argument_group('Required arguments')
	requiredGroup.add_argument('-p','--path',metavar='<path>', type=str,\
		help='Path of the SimPhy folder.', required=True)
	requiredGroup.add_argument('-ip','--input-prefix', metavar='<input_prefix>', type=str,\
		help='Prefix of the FASTA filenames.', required=True)
	requiredGroup.add_argument('-op','--output-prefix', metavar='<output_prefix>', type=str,\
		help='Prefix for the output filename.', required=True)
	requiredGroup.add_argument('-m','--method', metavar="<method_code>",type=int,\
		choices=range(0,5), default=0,\
		help="Specified method to obtain the reference loci used for the design of probes. Values range from 0-4.",\
        required=True)
	requiredGroup.add_argument('-o','--output', metavar="<output_path>",type=str,\
		help="Path where output will be written. ",\
		required=True)
	optionalGroup= parser.add_argument_group('Optional arguments')
    optionalGroup.add_argument('-sdf','--seq-desc-file',metavar='<sequence_descriptions_file_path>', type=int,\
    default="",\
		help="When method = 4 has been selected, it is required to identify a file with the sequence descriptions, "+
        "used as identifiers, corresponding to the sequence that will be used as reference per replicate.")
	optionalGroup.add_argument('-n','--nsize',metavar='<N_seq_size>', type=int,\
		default=-1,\
		help="Number of N's that will be introduced to separate the reference sequences selected. "+\
            "If the parameter is not set, the output file per replicate will be a multiple alignment sequence file, "+\
            "otherwise, the output will be a single sequence file per replicate consisting of a concatenation  "+\
            "of the reference sequences selected separated with as many N's as set for this parameter.")
	optionalGroup.add_argument('-l','--log',metavar='<log_level>', type=str,\
		choices=mll.LOG_LEVEL_CHOICES, default="INFO",\
		help='Specified level of log that will be shown through the standard output. Entire log will be stored in a separate file. Values:{0}. Default: {1}. '.format(mll.LOG_LEVEL_CHOICES,mll.LOG_LEVEL_CHOICES[1]))
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
		parser.print_help()
		sys.exit(-1)
	return tmpArgs

"""
MAIN
"""
if __name__=="__main__":
	try:
		cmdArgs = handlingParameters()
		print(cmdArgs)
		prog = lrf.LociReferenceSelection(cmdArgs)
		prog.run()
        sys.exit(0)
	except KeyboardInterrupt:
		sys.stdout.write("{0}{1}\nInterrupted!{2}\nPlease run again for the expected outcome.\n".format(mof.BOLD,mof.DARKCYAN, mof.END))
		sys.exit(-1)
_
