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

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if __name__ == '__main__':

    print ('using memsurfer from ({})'.format(memsurfer.__file__))

    ddir = './data'
    filename = os.path.join(ddir, 'noisy.off')

    # --------------------------------------------------------------------------
    # Create a logger
    memsurfer.utils.create_logger(2,1,0,'','')

    # The prefix we will use to name the output files
    outprefix = filename
    outprefix = outprefix[:outprefix.rfind('.')]

    # --------------------------------------------------------------------------
    # read a mesh (discard the faces, we only need vertices)
    verts, faces = memsurfer.utils.read_off(filename)
    verts = np.array(verts)

    LOGGER.info('{} : vertices = {}'.format(filename, verts.shape))

    # whether a periodic mesh
    periodic = True

    # bounding box (used only in the periodic case)
    box = np.array([-1.01,-1.01,-1.01,1.01,1.01,1.01]).reshape(2,3)

    # --------------------------------------------------------------------------
    # assign random labels to the vertices (optional use)
        # labels can be any datatype, e.g., strings or ints
    labels = np.random.randint(0, 5, size=(verts.shape[0]))
    labels = np.array([str(l) for l in labels])

    # --------------------------------------------------------------------------
    # create the membrane object
    knbrs = 18

    m = memsurfer.Membrane(verts, periodic=periodic, bbox=box, labels=labels)

    # fit the points to the box (optional)
    fit_needed = False
    if periodic and fit_needed:
        m.fit_points_to_box_xy()

    # estimate point normals (using k nearest neighbors)
    m.compute_pnormals(knbrs)

    # fit an approximating surface to the points
    m.compute_approx_surface()

    # compute the final membrane surfaces
    m.compute_membrane_surface()

    # compute properties on the final mesh
    m.compute_properties('exact')
    m.compute_properties('smooth')

    # or, you can directly access the triangulation to compute what you want
    '''
    self.memb_smooth.compute_normals()
    self.memb_smooth.compute_pointareas()
    self.memb_smooth.compute_curvatures()
    '''

    # --------------------------------------------------------------------------
    # estimate density of the points (all points, or selected on labels)
    sigmas = [0.2]
    get_nlipdis = True

    memsurfer.Membrane.compute_densities([m], [2], sigmas, get_nlipdis, 'all')
    memsurfer.Membrane.compute_densities([m], [2], sigmas, get_nlipdis, '1')

    # write the output
        # (required) a prefix for output files
        # (optional) a dictionary for any other parameters you may want to save
    m.write_all(outprefix, {'frame': 0, 'time': 0.0})

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
