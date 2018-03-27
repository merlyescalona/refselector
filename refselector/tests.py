import argparse,datetime,logging,os,sys, refselector, msatools
import loggingformatter as lf
import numpy as np
import random as rnd

def method1():
    try:
        cmdArgs=argparse.Namespace(help=False, input_prefix='data', log='INFO', method=0, nsize=300, output='reference', output_prefix='reference.300', seq_desc_file=None, simphy_path='testwsimphy/')
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