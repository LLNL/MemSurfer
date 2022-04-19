"""
Copyright (c) 2019, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.
Written by Harsh Bhatia (hbhatia@llnl.gov) and Peer-Timo Bremer (bhatia4@llnl.gov)
LLNL-CODE-763493. All rights reserved.

This file is part of MemSurfer, Version 1.0.
Released under GNU General Public License 3.0.
For details, see https://github.com/LLNL/MemSurfer.
"""

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
import os
import sys
import numpy as np
import argparse
import logging
LOGGER = logging.getLogger(__name__)

import memsurfer
from memsurfer import utils

import MDAnalysis
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.analysis.leaflet import LeafletFinder, optimize_cutoff

import mdreader
from lipidType import *

# ------------------------------------------------------------------------------
# This example script demonstrates how MemSurfer can be used to compute
# membrane surfaces.
# ------------------------------------------------------------------------------
if __name__ == '__main__':

    print(f'using memsurfer from ({os.path.basename(memsurfer.__file__)})')

    # --------------------------------------------------------------------------
    # Path to input data.
    # Currently, we're using hardcoded data, but it can be passed as argument
    '''
    parser = argparse.ArgumentParser(description='Example script for MemSurfer')
    parser.add_argument('--traj', required=True, nargs=1, help='Trajectory file')
    parser.add_argument('--topo', required=True, nargs=1, help='Topology file')
    args = vars(parser.parse_args())
    '''
    args = dict()
    ddir = './data'
    args['traj'] = [os.path.join(ddir, '10us.35fs-DPPC.40-DIPC.30-CHOL.30.gro')]
    args['topo'] = [os.path.join(ddir, '10us.35fs-DPPC.40-DIPC.30-CHOL.30.tpr')]

    # this data has three lipids
    lipids = ['DIPC', 'DPPC', 'CHOL']

    # --------------------------------------------------------------------------
    # Create a logger
    memsurfer.utils.create_logger(2,1,0,'','')

    # The prefix we will use to name the output files
    outprefix = args['traj'][0]
    outprefix = outprefix[:outprefix.rfind('.')]

    # --------------------------------------------------------------------------
    # Use MDAnalysis to read the data
    fargs = ['-f', args['traj'][0], '-s', args['topo'][0]]
    LOGGER.info(f'arguments = {fargs}')

    syst = mdreader.MDreader(fargs)
    syst.add_argument('-sframe', metavar='SELFRAME', type=int, dest='selframe', default=-1, help='int \tframe to save')
    syst.do_parse()

    LOGGER.info(f'Number of frames: {len(syst)}')
    LOGGER.info(f'System dimensions: {syst.dimensions[0:3]}')

    # In this example, we only look at a single frame (the last one)
    selFrame = [len(syst)-1]

    # Our domain is periodic (in xy)
    periodic = True

    # And this is the bounding box
    bbox=np.zeros((2,3))
    bbox[1,:] = syst.dimensions[:3]

    # --------------------------------------------------------------------------
    # Use a lipid master list to identify headgrups
    lipidTypeList = lipid_masterlist()

    # Get list of all suported lipid types
    lipidTypeNames = np.array([l[0] for l in lipidTypeList])

    # Get all resnames in this simulation
    resnames = np.unique(syst.atoms.resnames)

    # Make a lipid type dictionary
    mresnames = resnames[np.logical_not(np.in1d(resnames, lipidTypeNames))]
    LOGGER.warning(f'Following resnames either not lipids or not supported {mresnames}')

    lipidTypeDic = {}
    for i in np.where(np.in1d(lipidTypeNames, resnames))[0]:
        lipidTypeDic[lipidTypeNames[i]] = LipidType(lipidTypeList[i], syst)

    # Get beads that define the heads of all flip-flopping none leaflet defining lipids
    defFlipFlopHeadgroups = MDAnalysis.core.groups.AtomGroup([], syst)

    for i in list(lipidTypeDic.values()):
        defFlipFlopHeadgroups += i.getNoneLeafletSelection(syst)

    # --------------------------------------------------------------------------
    # Define top/bottom leaflets

    # Get beads that define the bilayer leaflets
    # (not in bilayer middle, but closly packed, use linker + first tail bead of lipids that do not flip-flop)
    defHeadgroups = MDAnalysis.core.groups.AtomGroup([], syst)
    for i in list(lipidTypeDic.values()):
        defHeadgroups += i.getLeafletSelection(syst)

    # get leaflets
    rcutoff, n = optimize_cutoff(syst, defHeadgroups)
    lfls = LeafletFinder(syst, defHeadgroups, cutoff=rcutoff, pbc=True)

    grps = lfls.groups()

    if len(grps) == 2:
        LOGGER.info(f'Found {len(grps)} leaflet groups')

    # check if they're even
    top_head = grps[0]
    bot_head = grps[1]

    rt = float(len(top_head))/len(bot_head)
    if rt > 1.3 or rt < 0.77:
        raise ValueError(f'Found uneven leaflets. '
                         f'top = {len(top_head)}, bot = {len(bot_head)}')

    LOGGER.info(f'Found leaflets of size {len(top_head)} and {len(bot_head)}')

    # --------------------------------------------------------------------------
    # Compute membranes for all frames
    for i in selFrame:

        print('\n')
        # Select frame 0 to len(syst)
        ts = syst.trajectory[i]
        LOGGER.info(f'Frame: {ts.frame:5d}, Time: {syst.trajectory.time:8.3f} ps')

        # Get all lipids in top/bot leaflets (including flip-flop lipids - therefore has to be done for each frame)
        tp = top_head + defFlipFlopHeadgroups.select_atoms('around 12 global group topsel', topsel=top_head)
        bt = bot_head + defFlipFlopHeadgroups.select_atoms('around 12 global group botsel', botsel=bot_head)

        if len(tp.select_atoms('group bt', bt=bt)):
            errmsg = f'Frame {syst.trajectory.frame}: ' \
                     f'{len(tp.select_atoms("group bt", bt=bt))} common atoms between leaflets'
            LOGGER.warning(errmsg)
            raise ValueError(errmsg)

        # Select one bead per lipid (could do this better by getting the residues - and selecting e.g. first of each
        LOGGER.info(f'We have {len(tp)} lipids in upper and {len(bt)} in lower leaflets')

        # ----------------------------------------------------------------------
        # This is where we use MemSurfer
        # ----------------------------------------------------------------------
        mt = memsurfer.Membrane.compute(tp.positions, tp.resnames.astype('U'), bbox, periodic)
        mb = memsurfer.Membrane.compute(bt.positions, bt.resnames.astype('U'), bbox, periodic)

        # compute total density
        sigmas = [10,20,30,40]
        get_nlipdis = True
        memsurfer.Membrane.compute_densities([mt, mb], [2], sigmas, get_nlipdis, 'all')

        # compute density of each type of lipid
        for l in lipids:
            memsurfer.Membrane.compute_densities([mt, mb], [2], sigmas, get_nlipdis, 'all')

        if True:
            # compute density of each type of lipid
            mt.write_all(f'{outprefix}_f{ts.frame}-top',
                         {'frame': ts.frame, 'time': syst.trajectory.time})
            mb.write_all(f'{outprefix}_f{ts.frame}-bot',
                         {'frame': ts.frame, 'time': syst.trajectory.time})

        # --------------------------------------------------------------------------
        # --------------------------------------------------------------------------
