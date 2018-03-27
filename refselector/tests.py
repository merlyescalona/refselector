import argparse,datetime,logging,os,sys, refselector, msatools
import loggingformatter as lf
import numpy as np
import random as rnd

def method0():
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

def method1():
    try:
        cmdArgs=argparse.Namespace(help=False, input_prefix='data', log='INFO', method=1, nsize=300, output='reference.1', output_prefix='reference.300', seq_desc_file='sequence.description.txt', simphy_path='testwsimphy/')
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

def method2():
    try:
        cmdArgs=argparse.Namespace(help=False, input_prefix='data', log='INFO', method=2, nsize=300, output='reference', output_prefix='reference.300', seq_desc_file=None, simphy_path='testwsimphy/')
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

def method3():
    try:
        cmdArgs=argparse.Namespace(help=False, input_prefix='data', log='INFO', method=3, nsize=300, output='reference', output_prefix='reference.300', seq_desc_file=None, simphy_path='testwsimphy/')
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


def method4():
    try:
        cmdArgs=argparse.Namespace(help=False, input_prefix='data', log='INFO', method=4, nsize=300, output='reference', output_prefix='reference.300', seq_desc_file=None, simphy_path='testwsimphy/')
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

def allMethods():
    try:
        cmdArgs0=argparse.Namespace(help=False, input_prefix='data', log='INFO', method=0, nsize=300, output='reference', output_prefix='reference.300', seq_desc_file=None, simphy_path='testwsimphy/')
        prog0 = refselector.ReferenceSelection(cmdArgs0)
        prog0.run()
        cmdArgs1=argparse.Namespace(help=False, input_prefix='data', log='INFO', method=1, nsize=300, output='reference.1', output_prefix='reference.300', seq_desc_file='sequence.description.txt', simphy_path='testwsimphy/')
        prog1 = refselector.ReferenceSelection(cmdArgs1)
        prog1.run()
        cmdArgs2=argparse.Namespace(help=False, input_prefix='data', log='INFO', method=2, nsize=300, output='reference', output_prefix='reference.300', seq_desc_file=None, simphy_path='testwsimphy/')
        prog2 = refselector.ReferenceSelection(cmdArgs2)
        prog2.run()
        cmdArgs3=argparse.Namespace(help=False, input_prefix='data', log='INFO', method=3, nsize=300, output='reference', output_prefix='reference.300', seq_desc_file=None, simphy_path='testwsimphy/')
        prog3 = refselector.ReferenceSelection(cmdArgs3)
        prog3.run()
        cmdArgs4=argparse.Namespace(help=False, input_prefix='data', log='INFO', method=4, nsize=300, output='reference', output_prefix='reference.300', seq_desc_file=None, simphy_path='testwsimphy/')
        prog4 = refselector.ReferenceSelection(cmdArgs4)
        prog4.run()
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
