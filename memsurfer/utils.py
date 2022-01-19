"""
Copyright (c) 2019, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.
Written by Harsh Bhatia (hbhatia@llnl.gov) and Peer-Timo Bremer (bremer5@llnl.gov)
LLNL-CODE-763493. All rights reserved.

This file is part of MemSurfer, Version 1.0.
Released under GNU General Public License 3.0.
For details, see https://github.com/LLNL/MemSurfer.
"""

# ------------------------------------------------------------------------------
# This file provides some simple utilities used by MemSurfer
# ------------------------------------------------------------------------------

import os
import numpy as np
import timeit
import inspect
import logging

# ------------------------------------------------------------------------------
# logging utils
# ------------------------------------------------------------------------------
LOGGER = logging.getLogger(__name__)

# Configuration globals
LFORMAT = '%(asctime)s [%(levelname)s] %(name)s:%(funcName)s - %(message)s'
ACCEPTED_INPUT = set(["yes", "y"])


# Create an interface to initialize a logger
def create_logger(loglevel, logstdout, logfile, logpath, logfname):
    """
        logstdout -- enable logging to stdout!
        logfile   -- enable logging to file!
        logfname  -- logfile name!
        logpath   -- logpath name!
        loglevel  -- log level!
    """
    # log level!
    assert(loglevel >= 1 and loglevel <= 5)
    if loglevel == 1:       loglevel = logging.DEBUG
    elif loglevel == 2:     loglevel = logging.INFO
    elif loglevel == 3:     loglevel = logging.WARN
    elif loglevel == 4:     loglevel = logging.ERROR
    elif loglevel == 5:     loglevel = logging.CRITICAL

    ROOTLOGGER = logging.getLogger(inspect.getmodule(__name__))
    formatter = logging.Formatter(LFORMAT)
    ROOTLOGGER.setLevel(loglevel)

    # Add the StreamHandler
    if logstdout == 1:
        sh = logging.StreamHandler()
        sh.setLevel(loglevel)
        sh.setFormatter(formatter)
        ROOTLOGGER.addHandler(sh)
        LOGGER.debug('stdout Enabled')

    # Add the FileHandler
    if logfile == 1:
        logpath = os.path.expanduser(logpath)
        if not os.path.exists(logpath):
            os.makedirs(logpath)

        if logfname[-4:] != '.log':
            logfname += '.log'

        logfname = os.path.join(logpath, logfname)
        fh = logging.FileHandler(logfname)
        fh.setLevel(loglevel)
        fh.setFormatter(formatter)
        ROOTLOGGER.addHandler(fh)
        LOGGER.debug(f'file {logfname} Enabled')

    return LOGGER

    # Print the level of logging.
    LOGGER.debug('Enabled')
    LOGGER.info('Enabled')
    LOGGER.warning('Enabled')
    LOGGER.error('Enabled')
    LOGGER.critical('Enabled')
    return LOGGER


# ------------------------------------------------------------------------------
# timing utils
# ------------------------------------------------------------------------------
class Timer(object):

    def __init__(self):
        self.stime = None
        self.etime = None
        self.start()

    def start(self):
        self.stime = timeit.default_timer()

    def end(self, print_time=False):
        self.etime = timeit.default_timer()
        if print_time:
            print(self)

    def __str__(self):

        tseconds = self.etime - self.stime
        if tseconds < np.power(10.,-6.):
            return "%.3f micro-sec." % (tseconds*np.power(10,6))

        elif tseconds < np.power(10.,-3.):
            return "%.3f milli-sec." % (tseconds*np.power(10.,3.))

        elif tseconds < 60.0:
            return "%.3f sec." % (tseconds)

        else:
            m = int(tseconds/ 60.0)
            return "%d min. %.3f sec." % (m, (tseconds - 60*m))

    def __repr__(self):
        return self.__str__()

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
