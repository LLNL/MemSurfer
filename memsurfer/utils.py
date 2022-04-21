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
import scipy.spatial.distance


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
def test_overlapping_points(points, dthreshold = 0.00001, fix = False):

    # --------------------------------------------------------------------------
    # https://stackoverflow.com/questions/13079563/how-does-condensed-distance-matrix-work-pdist
    def calc_row_idx(k, n):
        return int(np.ceil((1 / 2.) * (- (-8 * k + 4 * n ** 2 - 4 * n - 7) ** 0.5 + 2 * n - 1) - 1))

    def elem_in_i_rows(i, n):
        return i * (n - 1 - i) + (i * (i + 1)) // 2

    def calc_col_idx(k, i, n):
        return int(n - elem_in_i_rows(i + 1, n) + k)

    def condensed_to_square(k, n):
        i = calc_row_idx(k, n)
        j = calc_col_idx(k, i, n)
        return i, j

    # --------------------------------------------------------------------------
    LOGGER.debug(f'Testing for overlapping points ({points.shape}, threshold={dthreshold}, fix={fix})')

    pdist = scipy.spatial.distance.pdist(points)
    is_zero = np.where(pdist < dthreshold)[0]
    nzero = is_zero.shape[0]

    if nzero == 0:
        return True

    is_zero = [condensed_to_square(zidx, points.shape[0]) for zidx in is_zero]
    LOGGER.warning(f'Found {nzero} overlapping points: {is_zero}')
    if not fix:
        return False

    # now, try to fix them
    # we want to move the fewest points possible
    dupls, dupls_cnts = np.unique([x for sublist in is_zero for x in sublist], return_counts=True)
    dupls = {dupls[i]: dupls_cnts[i] for i in range(len(dupls))}

    # try to separate each pair
    max_move = 0.5 * pdist[pdist > dthreshold].min()
    for x,y in is_zero:

        dist0 = np.linalg.norm(points[x]-points[y])

        # we have moved either of these points before
        if not np.isclose(dist0, 0.):
            continue
        if dupls[x] == -1 or dupls[y] == -1:
            continue

        # move the one that impacts more points
        if dupls[x] > dupls[y]:
            a, b = x, y
        else:
            a, b = y, x

        # move a
        pa = f'{points[a]}'

        rad = np.random.rand() * max_move
        theta = np.random.rand() * np.pi
        points[a,0] += rad * np.cos(theta)
        points[a,1] += rad * np.sin(theta)

        dist1 = np.linalg.norm(points[x] - points[y])
        LOGGER.debug(f'Moved point {a}: {pa} --> {points[a]} : dist = {dist0} --> {dist1}')

        # mark this as moved
        dupls[a] = -1


# ------------------------------------------------------------------------------
def get_order_param(a2,a3, b2,b3, z=np.array([0,0,1])):
    """Get lipid "order" parameter
    # (a2, a3) are the second and third beads in tail a
    # (b2, b3) are the second and third beads in tail b
    # z represents the normal at the head

    # P = 0.5 * (3*cos^2(theta) - 1) order param at pos3 and average over tails
    #    where "theta" is the angle between the bond and the bilayer normal.
    #    P2 = 1      perfect alignment with the bilayer normal
    #    P2 = -0.5   anti-alignment
    #    P2 = 0      random orientation

    # returns (n,2) array with a column for each tail
    """

    if z.ndim == 1:
        z = np.array([z])

    LOGGER.info(f'Computing order parameter'
                f'({a2.shape}, {a3.shape}, {b2.shape}, {b3.shape}; {z.shape})')

    assert a2.shape == a3.shape and \
           a2.shape == b2.shape and \
           a2.shape == b3.shape, f'Tail\'s shapes must match'

    assert z.shape == (1,3) or z.shape == a2.shape, \
        'Should provide a single normal or a normal per lipid'

    # compute the norm of z once
    norm_z = np.linalg.norm(z, axis=1)

    # smaller function that implements the formula for order parameter
    def _compute(dx):
        dx_dot_z = np.einsum('ij, ij->i', dx, z)
        costht = dx_dot_z / (np.linalg.norm(dx, axis=1) * norm_z)
        return 0.5 * (3.0 * costht * costht - 1)

    # compute for both vectors
    orda = _compute(a2-a3)
    ordb = _compute(b2-b3)
    return np.swapaxes(np.vstack((orda, ordb)), 0,1)


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
