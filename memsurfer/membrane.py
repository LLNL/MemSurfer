'''
Copyright (c) 2019, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.
Written by Harsh Bhatia (hbhatia@llnl.gov) and Peer-Timo Bremer (bremer5@llnl.gov)
LLNL-CODE-763493. All rights reserved.

This file is part of MemSurfer, Version 1.0.
Released under GNU General Public License 3.0.
For details, see https://github.com/LLNL/MemSurfer.
'''

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
import sys
import numpy as np

from trimesh import TriMesh
from utils import Timer

import pysurfer
from pypoisson import poisson_reconstruction

import logging
LOGGER = logging.getLogger(__name__)

# ------------------------------------------------------------------------------
class Membrane(object):

    # --------------------------------------------------------------------------
    # constructor
    def __init__(self, points, **kwargs):
        '''
        points: ndarray of shape (npoints, 3)
        kwargs:
                periodic:       boolean flag for periodicity (default: False)
                boundary_layer: float value for boundary layer thickess (default: 0.2)
                                    used only for periodic domain
                bbox:           ndarray of shape (nverts, 3)
                                    required for periodic domain
                labels:         label for each point
        '''
        # 3d points
        if points.shape[1]!= 3:
            raise ValueError('Membrane needs 3D points: ndarray (npoints, 3)')

        self.points = points
        self.npoints = self.points.shape[0]

        # periodicity information
        self.periodic = kwargs.get('periodic', False)
        if self.periodic:
            self.blayer = kwargs.get('boundary_layer', 0.2)
            if 'bbox' not in kwargs.keys():
                raise ValueError('Periodic membrane needs 3D bounding box: ndarray (2 ,3)')

        # bounding box may be not given, but is needed for periodicity
        if 'bbox' in kwargs.keys():
            self.bbox = kwargs['bbox']
            if self.bbox.shape[0]!= 2 and self.bbox.shape[1] != 3:
                raise ValueError('Membrane needs 3D bounding box: ndarray (2 ,3)')

        # labels for points
        if 'labels' in kwargs.keys():
            self.labels = kwargs['labels']
            if self.labels.shape[0] != self.npoints or len(self.labels.shape) > 1:
                raise ValueError('Membrane expects one label per point')
        else:
            self.labels = np.empty((0,0))

        LOGGER.info('Initializing Membrane with {} points'.format(self.points.shape))
        LOGGER.info('\t actual bbox   = {} {}'.format(self.points.min(axis=0),self.points.max(axis=0)))
        if self.periodic  and 'bbox' in kwargs.keys():
            LOGGER.info('\t given periodic bbox = {} {}'.format(self.bbox[0],self.bbox[1]))

        # create point set object

        self.points = self.points.astype(np.float32)
        self.bbox = self.bbox.astype(np.float32)
        self.pset = pysurfer.PointSet(self.points)

        # point normals
        self.properties = {}
        self.pnormals = np.empty((0,0))

    # --------------------------------------------------------------------------
    def fit_points_to_box_xy(self):
        '''
            Fit the points in a periodic domain to the given bounding box
        '''
        # move the points to make sure the points lie in the box
        nadjusted = 0
        boxw = self.bbox[1,:] - self.bbox[0,:]

        for d in xrange(2):

            l = np.where(self.points[:,d] < self.bbox[0,d])[0]
            if l.shape[0] > 0:
                self.points[l,d] = self.points[l,d] + boxw[d]
                nadjusted += l.shape[0]

            r = np.where(self.points[:,d] > self.bbox[1,d])[0]
            if r.shape[0] > 0:
                self.points[r,d] = self.points[r,d] - boxw[d]
                nadjusted += r.shape[0]

        if nadjusted > 0:
            LOGGER.info('\t adjusted bbox = {} {}'.format(self.points.min(axis=0), self.points.max(axis=0)))

        return nadjusted

    # --------------------------------------------------------------------------
    def need_pnormals(self, knbrs=18, ndir_hint=+1):
        '''
            Estimate normals for the point set, using k nbrs
        '''
        if self.pnormals.shape != (0,0):
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
            raise ValueError('PointSet computed incorrect normals (got {} values}'.len(normals))

        self.pnormals = np.array(normals).reshape(self.npoints, 3)
        self.pnormals = self.pnormals.astype(np.float32)

        # if the max absolute z doesnt match ndir, then need to invert
        maxAbsZ = np.argmax(np.abs(self.pnormals[:,2]))
        if self.pnormals[maxAbsZ, 2] * ndir_hint < 0:
            LOGGER.debug('Inverting normals')
            self.pnormals = -1 * self.pnormals

        mtimer.end()
        LOGGER.info('Computed {} normals! took {}'.format(self.pnormals.shape, mtimer))
        LOGGER.info('\t normals = {} {}'.format(self.pnormals.min(axis=0),self.pnormals.max(axis=0)))
        return self.pnormals

    # --------------------------------------------------------------------------
    def fit_surface(self, depth=10):
        '''
            Compute an approximating surface using Poisson recronstruction
                depth:  controls the smoothness
                        larger depth will fit the points more --> less smooth
        '''
        self.need_pnormals()

        LOGGER.info('Computing Poisson surface for {} points'.format(self.npoints))
        mtimer = Timer()
        sfaces, sverts = poisson_reconstruction(self.points.tolist(),
                                                self.pnormals.tolist(),
                                                depth=depth)
        mtimer.end(False)
        LOGGER.info('Poisson surface computed! took {} created {} faces and {} vertices.'
                    .format(mtimer, len(sfaces), len(sverts)))

        # represent the poisson surface as a triangulation
        self.mesh3poisson = TriMesh(sverts, faces=sfaces, label='mesh3poisson')
        #self.mesh3poisson.write_vtp('_temp3.vtp', {})#params)
        #self.mesh3poisson.remesh()
        #self.mesh3poisson.write_vtp('_temp3_remeshed.vtp', {})#params)
        self.mesh3poisson.display()

    # --------------------------------------------------------------------------
    # retriangulate the surface by projecting a set of points on it
    def retriangulate(self, points):

        # 1. parameterize the membrane surface
        self.mesh3poisson.parameterize()

        #self.mesh2poisson = TriMesh(self.mesh3poisson.pverts, faces=self.mesh3poisson.faces)
        #self.mesh2poisson.write_vtp('temp_mesh2poisson.vtp')

        # 2. project the points on the surface and 2D plane
        (self.spoints, self.ppoints) = self.mesh3poisson.project_on_surface_and_plane(points)

        # 3. create a 2D triangulation of the projected points
        self.mesh2 = TriMesh(self.ppoints, periodic=self.periodic, label='mesh2')
        self.mesh32 = TriMesh(self.spoints, periodic=self.periodic, label='mesh32')
        self.mesh3 = TriMesh(self.points, periodic=self.periodic, label='mesh3')

        if self.periodic:

            # use the bbox of parameterized vertices (= bbox of poisson surface)
                # that is the absolute maximum
            bb0 = self.mesh3poisson.pverts.min(axis=0)
            bb1 = self.mesh3poisson.pverts.max(axis=0)
            self.mesh2.set_bbox(bb0, bb1)

            # bounding box of the points projected on the surface
            bb0 = self.spoints.min(axis=0)
            bb1 = self.spoints.max(axis=0)
            self.mesh32.set_bbox(bb0, bb1)

            # bounding box of the actual points
            bb0 = self.points.min(axis=0)
            bb1 = self.points.max(axis=0)
            self.mesh3.set_bbox(bb0, bb1)

        # compute delaunay triangulation of the planar points
        self.mesh2.delaunay()
        self.mesh32.copy_triangulation(self.mesh2)
        self.mesh3.copy_triangulation(self.mesh2)

        # ----------------------------------------------------------------------
        # set the correct bounding box
        if self.periodic:
            bb0 = self.ppoints.min(axis=0)
            bb1 = self.ppoints.max(axis=0)
            self.mesh2.set_bbox(bb0, bb1)

    # --------------------------------------------------------------------------
    # compute density of points given by plabels
        # on every vertex
    def estimate_density(self, sigma, name, normalize, labels=[]):

        nlabels = len(labels)

        # if labels are not available
        if nlabels == 0 and self.labels.shape == (0,0):
            raise ValueError('Cannot compute density of selected labels, because point labels are not available')

        # estimate density of all points
        if nlabels == 0:
            self.properties[name] = self.mesh2.estimate_density(sigma, name, normalize, np.empty([0]))
            return

        # estimate density of a subset of points
        # still compute on all vertices
        if self.labels.shape == (0,0):
            raise ValueError('Cannot compute density of selected labels, because point labels are not available')

        lidxs = np.where(np.in1d(self.labels, labels))[0]
        self.properties[name] = self.mesh2.estimate_density(sigma, name, normalize, lidxs)

    # --------------------------------------------------------------------------
    def write_all(self, outprefix, params={}):

        from utils import write2vtkpolydata

        pparams = {'normals': self.pnormals}
        if self.labels.shape != (0,0):
            pparams['labels'] = self.labels

        write2vtkpolydata(outprefix+'_points.vtp', self.points, pparams)
        self.mesh3poisson.write_vtp(outprefix+'_mesh3poisson.vtp')

        if self.labels.shape != (0,0):
            params['labels'] = self.labels

        for key in self.properties.keys():
            params[key] = self.properties[key]

        self.mesh2.write_vtp(outprefix+'_mesh2.vtp', params)
        self.mesh3.write_vtp(outprefix+'_mesh3.vtp', params)
        self.mesh32.write_vtp(outprefix+'_mesh32.vtp', params)

        self.mesh2.tmesh.write_binary(outprefix+'_mesh2.bin')
        self.mesh3.tmesh.write_binary(outprefix+'_mesh3.bin')
        self.mesh32.tmesh.write_binary(outprefix+'_mesh32.bin')

    def write(self, outprefix, params={}):

        from utils import write2vtkpolydata

        pparams = {'normals': self.pnormals}
        if self.labels.shape != (0,0):
            pparams['labels'] = self.labels

        write2vtkpolydata(outprefix+'_points.vtp', self.points, pparams)

        if self.labels.shape != (0,0):
            params['labels'] = self.labels

        for key in self.properties.keys():
            params[key] = self.properties[key]

        self.mesh32.write_vtp(outprefix+'_membrane.vtp', params)

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
        m.need_pnormals(knbrs)
        m.fit_surface()
        m.retriangulate(positions)

        # compute properties on the membrane (mesh3 is the final mesh)
        m.mesh3.need_normals()
        m.mesh3.need_pointareas()
        m.mesh3.need_curvatures()
        return m

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
