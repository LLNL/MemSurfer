#!/usr/bin/env python3
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
import glob
import numpy as np
import argparse

import logging
from datetime import datetime
LOGGER = logging.getLogger(__name__)

import memsurfer
from memsurfer import membrane
from memsurfer import utils

import MDAnalysis as mda
import pandas as pd

# Our domain is periodic (in xy)
periodic = True


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
def read_fatslim_to_df(filename):
    LOGGER.info(f"Reading ({filename})")
    df = pd.read_csv(filename)
    df = df.rename(columns={"X coords": "pos_x", "Y coords": "pos_y",
                            "Z coords": "pos_z", "Area per lipid": "apl(fatslim)"})
    df['leaflet'] = df['leaflet'].apply(lambda x: x.split(' ')[0])
    LOGGER.info(f"Read {df.shape} data")
    return df


# ------------------------------------------------------------------------------
def fetch_helices(u):

    # --------------------------------------------------------------------------
    # info provided by Tim
    nhelices = 7
    upper_helix_residues = [[7, 8, 9], [66, 67, 68], [78, 79, 80],
                            [144, 145, 146], [172, 173, 174], [258, 259, 260],
                            [267, 268, 269]]
    lower_helix_residues = [[27, 28, 29], [47, 48, 49], [103, 104, 105],
                            [126, 127, 128], [195, 196, 197], [233, 234, 235],
                            [285, 285, 287]]

    # --------------------------------------------------------------------------
    resids = {'upper': upper_helix_residues, 'lower': lower_helix_residues}
    resnames = {'upper': [], 'lower': []}
    coords = {'upper': np.zeros([nhelices, 3]), 'lower': np.zeros([nhelices, 3])}

    for nm in ['upper', 'lower']:
        for i in range(nhelices):
            sel = 'name BB and resid ' + ' '.join([str(h) for h in resids[nm][i]])
            coords[nm][i] = (u.select_atoms(sel).center_of_geometry() / 10)
            resnames[nm].append(f'helix_{nm[0]}{i}')
        resnames[nm] = np.array(resnames[nm])

    return resids, resnames, coords


# ------------------------------------------------------------------------------
def process_frame(trajfile, fslmfile, outfile):

    if not os.path.isfile(trajfile):
        LOGGER.warning(f'File ({trajfile}) not found!')
        return False

    if not os.path.isfile(fslmfile):
        LOGGER.warning(f'File ({fslmfile}) not found!')
        return False

    if os.path.isfile(outfile):
        LOGGER.info(f'File ({outfile}) already exists!')
        return False

    # --------------------------------------------------------------------------
    LOGGER.info(f"Reading ({trajfile})")

    u = mda.Universe(trajfile)
    LOGGER.info(f'Number of frames: {len(u.trajectory)}')
    LOGGER.info(f'System dimensions: {u.dimensions[0:3]}')
    bbox = np.zeros((2, 3))
    bbox[1, :] = u.dimensions[:3] * 0.1

    # to compute order parameter, select all pos2 and pos3 atoms
    ag_pos2a = u.select_atoms("name C2A D2A")
    ag_pos2b = u.select_atoms("name C2B D2B")
    ag_pos3a = u.select_atoms("name C3A D3A")
    ag_pos3b = u.select_atoms("name C3B D3B")

    # get positions of protein's helices (from mdanalysis data)
    prot_resids, prot_resnames, prot_coords = fetch_helices(u)

    nhelices = 7
    for nm in ['upper','lower']:
        assert nhelices == prot_coords[nm].shape[0]
        assert nhelices == prot_resnames[nm].shape[0]

    helix_nans = np.empty(nhelices)
    helix_nans.fill(np.nan)
    helix_resids = -1 * np.ones(nhelices, dtype=int)
    helix_idxs = np.arange(nhelices)

    LOGGER.info(f'Extracted {nhelices} helices')

    # --------------------------------------------------------------------------
    df = read_fatslim_to_df(fslmfile)

    # --------------------------------------------------------------------------
    # This is where we use MemSurfer
    # --------------------------------------------------------------------------
    mdf = []
    for nm in ['upper','lower']:

        LOGGER.info(f'Processing Bilayer {nm}')

        # extract the df with headgroup info for this leaflet from fatslim data
        hg_df = df[df['leaflet'] == nm]
        LOGGER.info(f'leaflet {nm}: for {hg_df.shape} dataframe')

        # get resids and positions from the dataframe
        hg_resids = hg_df['resid'].to_numpy()
        hg_pos = np.dstack((hg_df['pos_x'].to_numpy(),
                            hg_df['pos_y'].to_numpy(),
                            hg_df['pos_z'].to_numpy()))[0]
        hg_apl = hg_df['apl(fatslim)'].to_numpy()

        # find residue names for all these residues using mdanalysis
        hg = u.select_atoms('resid '+' '.join([str(h) for h in hg_resids]))
        hg_resnames = hg.residues.resnames

        nresidues = hg_resids.shape[0]
        nchol = sum(hg_resnames == 'CHOL')

        LOGGER.info(f'leaflet {nm}: '
                    f'nresidues = {nresidues}, nCHOL = {nchol} ',
                    f'unique values = {np.unique(hg_resnames)}')

        # ----------------------------------------------------------------------
        # let's sort on residue names to put all of them together
        # this is useful for managing tails, which are not defined for CHOL
        resnames_sidx = np.argsort(hg_resnames)

        hg_pos = hg_pos[resnames_sidx]
        hg_resids = hg_resids[resnames_sidx]
        hg_resnames = hg_resnames[resnames_sidx]
        hg_apl = hg_apl[resnames_sidx]

        # make sure we dont use this again, since we have reordered the data
        hg_df = None

        # ----------------------------------------------------------------------
        # now, extract tail for all non-CHOL residues
        tsearch = 'resid ' + ' '.join([str(h) for h in hg_resids[nchol:]])

        tail_2a = ag_pos2a.select_atoms(tsearch)
        tail_2b = ag_pos2b.select_atoms(tsearch)
        tail_3a = ag_pos3a.select_atoms(tsearch)
        tail_3b = ag_pos3b.select_atoms(tsearch)

        # data extraction and wrangling for the computation for order parameters
        # when selecting tails, filter out CHOL to track the tails wrt heads
        LOGGER.info(f'leaflet {nm}: found tail atoms for non-CHOL residues: '
                    f'({tail_2a.n_atoms}, {tail_2b.n_atoms},'
                    f' {tail_3a.n_atoms}, {tail_3b.n_atoms})')

        assert tail_2a.n_atoms == (nresidues-nchol) and \
               tail_2a.n_atoms == tail_2b.n_atoms and \
               tail_2a.n_atoms == tail_2b.n_atoms and \
               tail_2a.n_atoms == tail_2b.n_atoms, \
               f'Found mismatching number of atoms for tails: ({nresidues-nchol} ' \
               f'{tail_2a.n_atoms}, {tail_2b.n_atoms}, ' \
               f'{tail_2b.n_atoms}, {tail_2b.n_atoms}'

        # ----------------------------------------------------------------------
        # extract protein positions and resnames
        LOGGER.info(f'Bilayer {nm}: '
                    f'#headgroups = {hg_pos.shape[0]}, '
                    f'#helices = {prot_coords[nm].shape[0]}')

        # concatenate positions and resnames for headgroups and proteins
        pos = np.concatenate((prot_coords[nm], hg_pos))
        res = np.concatenate((prot_resnames[nm], hg_resnames)).astype('U')

        # ----------------------------------------------------------------------
        # compute membrane properties
        m = memsurfer.membrane.Membrane.compute(pos, res, bbox, periodic)
        m.memb_smooth.compute_normals()
        m.memb_smooth.compute_pointareas()
        m.memb_smooth.compute_curvatures()
        m.memb_smooth.compute_shells(helix_idxs)

        # order parameters
        z = m.memb_smooth.attributes['pnormal']
        ord_param_chol = np.ones((nchol, 2))
        ord_param_helices = -1*np.ones((nhelices, 2))
        ord_param_lips = utils.get_order_param(tail_2a.atoms.positions,
                                               tail_3a.atoms.positions,
                                               tail_2b.atoms.positions,
                                               tail_3b.atoms.positions,
                                               z[nhelices+nchol:])
        order_params = np.concatenate((ord_param_helices, ord_param_chol, ord_param_lips))

        # grab densities and labels from the membrane!
        params = {key: val for key, val in m.attributes.items()}
        params['order_param_a'] = order_params[:,0]
        params['order_param_b'] = order_params[:,1]
        params['resnames'] = res
        params['respos'] = pos

        # for comparison and reference, we want to look at fatslim apl
        if 1:
            params['apl(fatslim)'] = np.concatenate((helix_nans, hg_apl))
            params['resids'] = np.concatenate((helix_resids, hg_resids))

        if 1:
            m.memb_smooth.write_vtp(os.path.join(path_out, f'{frame}-{nm}.vtp'), params)
            #m.write_all(os.path.join(path_out, f'all-{frame}-{nm}'), params)

        # ----------------------------------------------------------------------
        tdf = m.memb_smooth.as_pddataframe(additional_attributes=params)
        tdf['leaflet'] = nm
        mdf.append(tdf)

    # --------------------------------------------------------------------------
    mdf = pd.concat(mdf).reset_index(drop=True)

    # reorder and write
    if 1:
        mdf = mdf[['leaflet', 'resids', 'resnames',
                   'respos_x', 'respos_y', 'respos_z',
                   'pos_x', 'pos_y', 'pos_z',
                   'pnormal_x', 'pnormal_y', 'pnormal_z',
                   'shell', 'parea', 'mean_curv', 'gaus_curv',
                   'apl(fatslim)','order_param_a','order_param_b'
                   ]]

    LOGGER.info(f'Writing {mdf.shape} dataframe to ({outfile})')
    mdf.to_csv(outfile, sep=',', index=False, float_format='%.6f', na_rep='-1')
    return True


# ------------------------------------------------------------------------------
# This example script demonstrates how MemSurfer can be used to compute
# membrane surfaces.
# ------------------------------------------------------------------------------
if __name__ == '__main__':

    # --------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description='Process a (set of) frames using memsurfer')

    parser.add_argument('--outpath', dest='outpath', required=True,
                        help='Output path')

    # either provide a directory to process all files
    parser.add_argument('--inpath', dest='inpath',
                        help='Input path')

    # or a single frame (both files)
    parser.add_argument('--filegro', dest='ingro',
                        help='Path to input trajectory file (*.gro)')

    parser.add_argument('--filefat', dest='infat',
                        help='Path to input fatslim file (*.csv)')

    args = parser.parse_args()

    mode_frame = (args.ingro is not None) and (args.infat is not None)
    mode_traj = args.inpath is not None

    if mode_traj == mode_frame:
        print ('Need either a frame (gro + fatslim) or an input directory!\n')
        parser.print_help()
        exit(1)

    # --------------------------------------------------------------------------
    #path_out = os.path.join(args.outpath, 'shellanalysis', 'memsurfer')
    path_out = args.outpath
    os.makedirs(path_out, exist_ok=True)

    # --------------------------------------------------------------------------
    print(f'using memsurfer from ({os.path.dirname(memsurfer.__file__)})')

    if 1:
        memsurfer.utils.create_logger(1, 1,0, '','')
    else:
        ts = datetime.now().strftime('%Y%m%d-%H%M%S')
        memsurfer.utils.create_logger(2, 0,1, '.',f'memsurfer-{ts}.log')

    # --------------------------------------------------------------------------
    # find files to read
    # --------------------------------------------------------------------------
    gro2frame = lambda f: int(os.path.basename(f)[5:-4])
    fat2frame = lambda f: int(os.path.basename(f)[:-4].split('_')[-1])

    files = {}

    if mode_frame:
        print (args.ingro)
        t1 = gro2frame(args.ingro)
        print (args.infat)
        t2 = fat2frame(args.infat)
        assert t1 == t2, f'Input files have mismatch in frame number: {t1} != {t2}'

        files[t1] = [args.ingro, args.infat]

    else:
        path_in = args.inpath
        path_gro = os.path.join(path_in, 'frames')
        path_fslm = os.path.join(path_in, 'output')

        files_gro = sorted(glob.glob(os.path.join(path_gro, 'frame*.gro')))
        files_fslm = sorted(glob.glob(os.path.join(path_fslm, 'fatslim-output.frame_*.csv')))

        print (f'found {len(files_gro)} gro files in ({path_gro})')
        print (f'found {len(files_fslm)} fatslim files in ({path_fslm})')

        for f in files_gro:
            t = gro2frame(f)
            files[t] = [f, '']

        for f in files_fslm:
            t = fat2frame(f)
            if t in files:
                files[t][1] = f
            else:
                files[t] = ['', f]

    # --------------------------------------------------------------------------
    # Compute membranes for all frames
    # --------------------------------------------------------------------------
    for frame, files in files.items():

        trajfile = files[0]
        fslmfile = files[1]
        trajname = os.path.basename(files[0])[:-4]
        outfile = os.path.join(path_out, f'{trajname}-memsurfer.csv')

        try:
            process_frame(trajfile, fslmfile, outfile)

        except Exception as e:
            LOGGER.error(f'Failed for ({trajname}): {e}')

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
