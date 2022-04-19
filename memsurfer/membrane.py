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
# ------------------------------------------------------------------------------
import sys
import numpy as np
import logging
LOGGER = logging.getLogger(__name__)

from pypoisson import poisson_reconstruction
from . import memsurfer_cmod
from .trimesh import TriMesh
from .utils import Timer, test_overlapping_points


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
class Membrane(object):
    """
       Class to create and manipulate membrane surfaces
    """

    # --------------------------------------------------------------------------
    # constructor
    def __init__(self, points, **kwargs):
        """
        points: ndarray of shape (npoints, 3)
        kwargs:
                labels:         label for each point
                periodic:       boolean flag for periodicity (default: False)
                bbox:           ndarray of shape (nverts, 3)
                                    required for periodic domain
                boundary_layer: float value for boundary layer thickess (default: 0.2)
                                    used only for periodic domain
        """
        # 3d points
        if points.shape[1] !=  3:
            raise ValueError(f'Membrane needs 3D points: got {points.shape}')

        self.points = points
        self.npoints = self.points.shape[0]

        # periodicity information
        self.periodic = kwargs.get('periodic', False)
        if self.periodic:
            self.blayer = kwargs.get('boundary_layer', 0.2)
            if 'bbox' not in list(kwargs.keys()):
                raise ValueError('Periodic membrane needs 3D bounding box: ndarray (2 ,3)')

        # bounding box may be not given, but is needed for periodicity
        if 'bbox' in list(kwargs.keys()):
            self.bbox = kwargs['bbox']
            if self.bbox.shape[0] not in [2, 3]:
                raise ValueError('Membrane needs 3D bounding box: ndarray (2 ,3)')

        # labels for points
        if 'labels' in list(kwargs.keys()):
            self.labels = kwargs['labels']
            if self.labels.shape[0] != self.npoints or len(self.labels.shape) > 1:
                raise ValueError('Membrane expects one label per point')
        else:
            self.labels = None

        LOGGER.info(f'Initializing Membrane with {self.points.shape} points')
        LOGGER.info(f'\t actual bbox = {self.points.min(axis=0)} {self.points.max(axis=0)}')
        if self.periodic and 'bbox' in list(kwargs.keys()):
            LOGGER.info(f'\t given periodic bbox = {self.bbox[0]} {self.bbox[1]}')

        # create point set object
        self.points = self.points.astype(np.float32)
        self.bbox = self.bbox.astype(np.float32)
        self.pset = memsurfer_cmod.PointSet(self.points)

        # check for overlapping points, but no fix!
        test_overlapping_points(self.points, fix=False)

        # different constructs we will compute
        self.pnormals = None
        self.spoints = None
        self.ppoints = None

        self.surf_poisson = None
        self.memb_planar = None
        self.memb_smooth = None
        self.memb_exact = None

        # other properties
        self.attributes = {}

    # --------------------------------------------------------------------------
    def fit_points_to_box_xy(self):
        """
            Fit the points in a periodic domain to the given bounding box
        """
        nadjusted = 0
        boxw = self.bbox[1,:] - self.bbox[0,:]
        for d in range(2):

            l = np.where(self.points[:,d] < self.bbox[0,d])[0]
            if l.shape[0] > 0:
                self.points[l,d] = self.points[l,d] + boxw[d]
                nadjusted += l.shape[0]

            r = np.where(self.points[:,d] > self.bbox[1,d])[0]
            if r.shape[0] > 0:
                self.points[r,d] = self.points[r,d] - boxw[d]
                nadjusted += r.shape[0]

        if nadjusted > 0:
            LOGGER.info(f'\t adjusted bbox = {self.points.min(axis=0)} {self.points.max(axis=0)}')

        return nadjusted

    # --------------------------------------------------------------------------
    def compute_pnormals(self, knbrs=18, ndir_hint=+1):
        """
            Estimate normals for the point set, using k nearest neighhbors
        """
        if self.pnormals is not None:
            return self.pnormals

        cverbose = LOGGER.isEnabledFor(logging.DEBUG)

        LOGGER.info('Computing normals')
        mtimer = Timer()

        # this function will duplicate points within the boundary layer
        if self.periodic:
            self.pset.set_periodic(self.bbox.reshape(-1).tolist(), self.blayer, cverbose)

        self.pset.need_normals(knbrs, cverbose)
        normals = self.pset.get_normals()
        if len(normals) != 3*self.npoints:
            raise ValueError(f'PointSet computed incorrect normals (got {len(normals)} values')

        self.pnormals = np.array(normals).reshape(self.npoints, 3)
        self.pnormals = self.pnormals.astype(np.float32)

        # if the max absolute z doesnt match ndir, then need to invert
        maxAbsZ = np.argmax(np.abs(self.pnormals[:,2]))
        if self.pnormals[maxAbsZ, 2] * ndir_hint < 0:
            LOGGER.debug('Inverting normals')
            self.pnormals *= -1.0

        mtimer.end()
        LOGGER.info(f'Computed {self.pnormals.shape} normals! took {mtimer}')
        LOGGER.info(f'\t normals = {self.pnormals.min(axis=0)}, {self.pnormals.max(axis=0)}')
        return self.pnormals

    # --------------------------------------------------------------------------
    def compute_approx_surface(self, exactness_level=10):
        """
            Compute an approximating surface using Poisson recronstruction
                exactness_level:  controls the smoothness
                        larger exactness_level will fit the points more --> less smooth
        """
        self.compute_pnormals()

        LOGGER.info(f'Computing Poisson surface for {self.npoints} points')
        mtimer = Timer()
        sfaces, sverts = poisson_reconstruction(self.points.tolist(),
                                                self.pnormals.tolist(),
                                                depth=exactness_level)
        mtimer.end()
        LOGGER.info(f'Poisson surface computed! took {mtimer}. '
                    f'created {len(sfaces)} faces and {len(sverts)} vertices.')

        # represent the poisson surface as a triangulation
        self.surf_poisson = TriMesh(sverts, faces=sfaces, label='surf_poisson')

        #self.surf_poisson.write_vtp('_temp3.vtp', {})#params)
        #self.surf_poisson.remesh()
        #self.surf_poisson.write_vtp('_temp3_remeshed.vtp', {})#params)
        LOGGER.info(self.surf_poisson)

    # --------------------------------------------------------------------------
    def compute_membrane_surface(self):
        """
            Compute the membrane surfaces that pass
                through the given set of points
        """

        # 1. parameterize the membrane surface
        self.surf_poisson.parameterize()

        # 2. project the points on the surface and 2D plane
        (self.spoints, self.ppoints) = self.surf_poisson.project_on_surface_and_plane(self.points)

        # 2b. handle overlapping points
        test_overlapping_points(self.spoints, fix=True)
        test_overlapping_points(self.ppoints, fix=True)

        # 3. create a 2D triangulation of the projected points
        self.memb_planar = TriMesh(self.ppoints, periodic=self.periodic, label='memb_planar')
        self.memb_smooth = TriMesh(self.spoints, periodic=self.periodic, label='memb_smooth')
        self.memb_exact  = TriMesh(self.points,  periodic=self.periodic, label='memb_exact')

        if self.periodic:

            # use the bbox of parameterized vertices (= bbox of poisson surface)
            # that is the absolute maximum
            bb0, bb1 = self.surf_poisson.parameterization_bbox()
            self.memb_planar.set_bbox(bb0, bb1)

            # bounding box of the points projected on the surface
            bb0 = self.spoints.min(axis=0)
            bb1 = self.spoints.max(axis=0)
            self.memb_smooth.set_bbox(bb0, bb1)

            # bounding box of the actual points
            bb0 = self.points.min(axis=0)
            bb1 = self.points.max(axis=0)
            self.memb_exact.set_bbox(bb0, bb1)

        # compute delaunay triangulation of the planar points
        self.memb_planar.compute_delaunay()
        self.memb_smooth.copy_triangulation(self.memb_planar)
        self.memb_exact.copy_triangulation(self.memb_planar)

        # ----------------------------------------------------------------------
        # change the bounding box of the planar surface
        if self.periodic:
            bb0 = self.ppoints.min(axis=0)
            bb1 = self.ppoints.max(axis=0)
            self.memb_planar.set_bbox(bb0, bb1)

    # --------------------------------------------------------------------------
    def compute_properties(self, mtype='smooth'):

        assert mtype == 'smooth' or mtype == 'exact'
        if mtype == 'exact':
            self.memb_exact.compute_normals()
            self.memb_exact.compute_pointareas()
            self.memb_exact.compute_curvatures()
        else:
            self.memb_smooth.compute_normals()
            self.memb_smooth.compute_pointareas()
            self.memb_smooth.compute_curvatures()

    # --------------------------------------------------------------------------
    # compute density of points given by plabels
        # on every vertex
    def compute_density(self, type, sigma, name, get_nlipdis, labels=[]):

        nlabels = len(labels)

        if type < 1 or type > 3:
            raise InvalidArgument(f'Invalid density type, {type}. Should be 1 (geodesic), 2 (2D) or 3 (3D)')

        # if labels are not available
        if nlabels == 0 and self.labels.shape is None:
            raise ValueError('Cannot compute density of selected labels, because point labels are not available')

        # estimate density of all points
        if nlabels == 0:
            self.attributes[name] = self.memb_smooth.compute_density(type, sigma, name, get_nlipdis, np.empty([0]))
            return

        # estimate density of a subset of points
        # still compute on all vertices
        if self.labels.shape is None:
            raise ValueError('Cannot compute density of selected labels, because point labels are not available')

        lidxs = np.where(np.in1d(self.labels, labels))[0]
        self.attributes[name] = self.memb_smooth.compute_density(type, sigma, name, get_nlipdis, lidxs)

    # --------------------------------------------------------------------------
    @staticmethod
    def compute_thickness(a, b, mtype='smooth'):

        assert mtype in ['smooth', 'exact']
        if mtype == 'exact':
            a.properties['thickness'] = a.memb_exact.compute_distance_to_surface(b.memb_exact)
            b.properties['thickness'] = b.memb_exact.compute_distance_to_surface(a.memb_exact)
        else:
            a.properties['thickness'] = a.memb_smooth.compute_distance_to_surface(b.memb_smooth)
            b.properties['thickness'] = b.memb_smooth.compute_distance_to_surface(a.memb_smooth)

    # --------------------------------------------------------------------------
    @staticmethod
    def compute_densities(membranes, types, sigmas, get_nlipdis, l):

        labels = [] if l == 'all' else [l]
        for t in types:
            for s in sigmas:
                name = 'density_type{0}_{1}_k{2:.1f}'.format(t, l, s)
                for m in membranes:
                    m.compute_density(t, s, name, get_nlipdis, labels)

    # --------------------------------------------------------------------------
    # A static method that computes and returns a membrane object
    # --------------------------------------------------------------------------
    @staticmethod
    def compute(positions, labels, bbox, periodic):

        knbrs = 18

        # initialize the membrane!
        m = Membrane(positions, labels=labels, periodic=periodic, bbox=bbox)

        # compute the membrane
        m.fit_points_to_box_xy()
        m.compute_pnormals(knbrs)
        m.compute_approx_surface()
        m.compute_membrane_surface()
        return m

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    def write_all(self, outprefix, params={}):

        pparams = {'normals': self.pnormals}
        if self.labels.shape != (0, 0):
            pparams['labels'] = self.labels

        if True:
            if 'shells' in self.memb_smooth.attributes: #self.memb_smooth.shells is not None:
                pparams['shells'] = self.memb_smooth.attributes['shells']

            # create a temporary mesh just to output the points
            TriMesh(self.points).write_vtp(outprefix+'_points.vtp', pparams)

        if True:
            self.surf_poisson.write_vtp(outprefix + '_surface_poisson.vtp')

        if True:
            for key in list(self.attributes.keys()):
                params[key] = self.attributes[key]

            if 'labels' in pparams.keys():
                params['labels'] = pparams['labels']

            self.memb_planar.write_vtp(outprefix+'_planar.vtp', params)
            self.memb_exact.write_vtp(outprefix+'_membrane_exact.vtp', params)
            self.memb_smooth.write_vtp(outprefix+'_membrane_smooth.vtp', params)

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
