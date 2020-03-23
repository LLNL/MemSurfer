'''
Copyright (c) 2019, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.
Written by Harsh Bhatia (hbhatia@llnl.gov) and Peer-Timo Bremer (bhatia4@llnl.gov)
LLNL-CODE-763493. All rights reserved.

This file is part of MemSurfer, Version 1.0.
Released under GNU General Public License 3.0.
For details, see https://github.com/LLNL/MemSurfer.
'''

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
import os, sys
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

    print ('using memsurfer from ({})'.format(memsurfer.__file__))

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
    LOGGER.info('arguments = {}'.format(fargs))

    syst = mdreader.MDreader(fargs)
    syst.add_argument('-sframe', metavar='SELFRAME', type=int, dest='selframe', default=-1, help='int \tframe to save')
    syst.do_parse()

    LOGGER.info('Number of frames in sim {}'.format(len(syst)))
    LOGGER.info('System dimensions: {}'.format(syst.dimensions[0:3]))

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
    LOGGER.warning('Following resnames either not lipids or not supported {}'
                    .format(resnames[np.logical_not(np.in1d(resnames, lipidTypeNames))]))

    lipidTypeDic = {}
    for i in np.where(np.in1d(lipidTypeNames, resnames))[0]:
        lipidTypeDic[lipidTypeNames[i]] = LipidType(lipidTypeList[i], syst)

    # Get beads that define the heads of all flip-flopping none leaflet defining lipids
    defFlipFlopHeadgroups = MDAnalysis.core.groups.AtomGroup([], syst)

    for i in list(lipidTypeDic.values()):
        defFlipFlopHeadgroups += i.getNoneLeafletSelection(syst)

    #LOGGER.info('FlipFlopHeadgroups = {}'.format(defFlipFlopHeadgroups))

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
        LOGGER.info('Found {} leaflet groups'.format(len(grps)))

    # check if they're even
    top_head = grps[0]
    bot_head = grps[1]

    rt = float(len(top_head))/len(bot_head)
    if rt > 1.3 or rt < 0.77:
        raise ValueError('Found uneven leaflets. top = {}, bot = {}'.format(len(top_head), len(bot_head)))

    LOGGER.info('Found leaflets of size {}, {}'.format(len(top_head), len(bot_head)))

    # --------------------------------------------------------------------------
    # Compute membranes for all frames
    for i in selFrame:

        print('\n')
        # Select frame 0 to len(syst)
        ts = syst.trajectory[i]
        LOGGER.info('Frame: %5d, Time: %8.3f ps' % (ts.frame, syst.trajectory.time))

        # Get all lipids in top/bot leaflets (including flip-flop lipids - therefore has to be done for each frame)
        tp = top_head + defFlipFlopHeadgroups.select_atoms('around 12 fullgroup topsel', topsel=top_head)
        bt = bot_head + defFlipFlopHeadgroups.select_atoms('around 12 fullgroup botsel', botsel=bot_head)

        if len(tp.select_atoms('group bt', bt=bt)):
            errmsg = 'Frame {}: {} common atoms between leaflets.'.format(syst.trajectory.frame, len(tp.select_atoms('group bt', bt=bt)))
            LOGGER.warning(errmsg)
            raise ValueError(errmsg)

        ## Select one bead per lipid (could do this better by getting the residues - and selecting e.g. first of each
        LOGGER.info('We have {} lipids in upper and {} in lower leaflets'.format(len(tp), len(bt)))

        # ----------------------------------------------------------------------
        # This is where we use MemSurfer
        # ----------------------------------------------------------------------
        mt = memsurfer.Membrane.compute(tp.positions, tp.resnames, bbox, periodic)
        mb = memsurfer.Membrane.compute(bt.positions, bt.resnames, bbox, periodic)

        # compute total density
        sigmas = [10,20,30,40]
        get_nlipdis = True
        memsurfer.Membrane.compute_densities([mt, mb], [2], sigmas, get_nlipdis, 'all')

        # compute density of each type of lipid
        for l in lipids:
            memsurfer.Membrane.compute_densities([mt, mb], [2], sigmas, get_nlipdis, 'all')

        if True:
            # compute density of each type of lipid
            mt.write_all(outprefix+'_f{}-top'.format(ts.frame),
                    {'frame': ts.frame, 'time': syst.trajectory.time})

            mb.write_all(outprefix+'_f{}-bot'.format(ts.frame),
                    {'frame': ts.frame, 'time': syst.trajectory.time})

        # --------------------------------------------------------------------------
        # --------------------------------------------------------------------------
