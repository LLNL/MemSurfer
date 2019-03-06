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

import numpy as np
import logging
LOGGER = logging.getLogger(__name__)

import pymemsurfer
from utils import Timer

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
class TriMesh(object):
    '''
       Class to create and manipulate triangular meshes (in 2D and 3D)
    '''

    # --------------------------------------------------------------------------
    # constructor
    def __init__(self, vertices, **kwargs):
        '''
        vertices: ndarray of shape (nverts, dim) # dim = 2,3
        kwargs:
                  faces: ndarray of shape (nfaces,3)
                  label: label for this mesh
                  periodic: whether this mesh is periodic
        '''
        if vertices.shape[1]!= 2 and vertices.shape[1]!= 3:
            raise ValueError('TriMesh needs 2D or 3D vertices: ndarray (npoints, 2/3)')

        self.vertices = vertices.astype(np.float32)
        self.nverts = self.vertices.shape[0]

        self.periodic = kwargs.get('periodic', False)
        if self.periodic:
            self.tmesh = pymemsurfer.TriMeshPeriodic(self.vertices)
            self.label = kwargs.get('label', 'TriMeshPeriodic')
        else:
            self.tmesh = pymemsurfer.TriMesh(self.vertices)
            self.label = kwargs.get('label', 'TriMesh')

        self.nfaces = 0
        if 'faces' in kwargs.keys():
            self.faces = kwargs['faces']
            self.faces = self.faces.astype(np.uint32)
            self.nfaces = self.faces.shape[0]

            if (self.faces.shape[1] != 3):
                raise ValueError('TriMesh needs triangles: ndarray (nfaces, 3)')

            self.tmesh.set_faces(self.faces)

        LOGGER.info('Created <{}> with {} vertices and {} faces'
                    .format(self.label, self.nverts, self.nfaces))

        # properties to be computed
        self.pnormals = np.empty((0,0))
        self.pareas = np.empty(0)
        self.mean_curv = np.empty(0)
        self.gaus_curv = np.empty(0)
        self.pverts = np.empty((0,0))
        self.cverbose = LOGGER.isEnabledFor(logging.DEBUG)

    # --------------------------------------------------------------------------
    def set_bbox(self, bb0, bb1):

        if bb0.shape[0] == 2 and bb1.shape[0] == 2:
            self.bbox = np.array([bb0[0], bb0[1],  bb1[0], bb1[1]])
            self.boxw = np.array([bb1[0] - bb0[0], bb1[1] - bb0[1]])

        elif bb0.shape[0] == 3 and bb1.shape[0] == 3:
            self.bbox = np.array([bb0[0], bb0[1], bb0[2], bb1[0], bb1[1], bb1[2]])
            self.boxw = np.array([bb1[0]-bb0[0], bb1[1]-bb0[1], bb1[2]-bb0[2]])

        else:
            raise ValueError('TriMesh got incorrect bbox: {}, {}'.format(bb0, bb1))

        self.bbox = self.bbox.astype(np.float32)
        self.boxw = self.boxw.astype(np.float32)
        LOGGER.info('<{}> setting bbox = {}'.format(self.label, self.bbox))
        self.tmesh.set_bbox(self.bbox)

    # --------------------------------------------------------------------------
    def copy_triangulation(self, mesh):

        assert(self.periodic == mesh.periodic)

        self.tmesh.set_faces(mesh.tmesh)
        self.faces = mesh.faces

        if self.periodic:
            self.pfaces = mesh.pfaces
            self.tfaces = mesh.tfaces
            dverts = self.tmesh.duplicated_vertices()
            self.dverts = np.array(dverts).reshape(-1, self.vertices.shape[1]).astype(np.float32)

    # --------------------------------------------------------------------------
    def parameterize(self, xy=False):

        if self.pverts.shape != (0,0):
            return self.pverts

        LOGGER.info('Parameterizing the surface')
        mtimer = Timer()

        if xy:
            pverts = self.tmesh.parameterize_xy(self.cverbose)
        else:
            pverts = self.tmesh.parameterize(self.cverbose)

        self.pverts = np.array(pverts).reshape(-1, 2).astype(np.float32)
        assert self.pverts.shape[0] == self.nverts

        mtimer.end()
        LOGGER.info('Parameterization took {}'.format(mtimer))
        return self.pverts

    # --------------------------------------------------------------------------
    def project_on_surface_and_plane(self, points):
        '''
        Project a set of points (given as barycentric coordinates and face ids)
            onto the triangulation, and parameterized triangulation
        '''
        self.parameterize()

        npoints = points.shape[0]
        LOGGER.info('Projecting {} points on the parameterized surface'.format(npoints))

        mtimer = Timer()

        r = self.tmesh.project_on_surface(points.reshape(-1).tolist(), self.cverbose)
        r = np.array(r).reshape(npoints, 4)

        fids = r[:,0].astype(int)
        fbarys = r[:,1:4]

        # create points using the projections
        spoints = np.zeros((npoints, 3), dtype=np.float32)
        ppoints = np.zeros((npoints, 2), dtype=np.float32)

        for i in xrange(npoints):

            f = self.faces[fids[i]]
            for j in xrange(3):

                spoints[i] += fbarys[i][j] * self.vertices[f[j]]
                ppoints[i] += fbarys[i][j] * self.pverts[f[j]]

        mtimer.end()
        LOGGER.info('Projecting took {}'.format(mtimer))
        LOGGER.info('\tbbox of planar projections: {}, {}'.format(ppoints.min(axis=0), ppoints.max(axis=0)))
        LOGGER.info('\tbbox of surface projections: {}, {}'.format(spoints.min(axis=0), spoints.max(axis=0)))
        return (spoints, ppoints)

    # --------------------------------------------------------------------------
    '''
    def remesh(self):

        LOGGER.info('Remeshing the surface with {} vertices ad {} faces'.format(self.nverts, self.nfaces))
        mtimer = Timer()

        self.tmesh.remesh()
        self.vertices = np.array(self.tmesh.vertices()).reshape(-1,3)
        self.faces = np.array(self.tmesh.faces()).reshape(-1,3)
        self.nverts = self.vertices.shape[0]
        self.nfaces = self.faces.shape[0]

        mtimer.end()
        LOGGER.info('Created {} vertices and {} faces! took {}'.format(self.nverts, self.nfaces, mtimer))
    '''

    # --------------------------------------------------------------------------
    def delaunay(self):

        if self.nfaces > 0:
            return

        tag = ' periodic ' if self.periodic else ' '
        LOGGER.info('Computing 2D{}Delaunay triangulation.'.format(tag))
        mtimer = Timer()

        faces = self.tmesh.delaunay(self.cverbose)
        self.faces = np.array(faces).reshape(-1, 3).astype(np.uint32)
        self.nfaces = self.faces.shape[0]

        if self.periodic:

            pfaces = self.tmesh.periodic_faces()
            self.pfaces = np.array(pfaces).reshape(-1, 3).astype(np.uint32)

            tfaces = self.tmesh.trimmed_faces()
            self.tfaces = np.array(tfaces).reshape(-1, 3).astype(np.uint32)

            dverts = self.tmesh.duplicated_vertices()
            self.dverts = np.array(dverts).reshape(-1, self.vertices.shape[1]).astype(np.float32)

            mtimer.end()
            LOGGER.info('Periodic Delaunay triangulation took {}! created [{} {} {}] faces, and [{} {}] vertices.'
                        .format(mtimer, self.nfaces, self.pfaces.shape[0], self.tfaces.shape[0],
                                        self.nverts, self.dverts.shape))

    # --------------------------------------------------------------------------
    def compute_normals(self):

        if self.pnormals.shape != (0,0):
            return self.pnormals

        LOGGER.info('Computing normals')
        mtimer = Timer()

        rval = self.tmesh.need_normals(self.cverbose)
        if len(rval) != 3*self.nverts:
            raise ValueError('Incorrect normals!')

        self.pnormals = np.asarray(rval).reshape(-1,3).astype(np.float32)

        mtimer.end()
        LOGGER.info('Computed {} normals! took {}'.format(self.nverts, mtimer))
        return self.pnormals

    # --------------------------------------------------------------------------
    def compute_pointareas(self):

        if self.pareas.shape != (0,):
            return self.pareas

        LOGGER.info('Computing point areas')
        mtimer = Timer()

        rval = self.tmesh.need_pointareas(self.cverbose)
        if len(rval) != self.nverts:
            raise ValueError('Incorrect point areas!')

        self.pareas = np.array(rval).astype(np.float32)

        mtimer.end()
        LOGGER.info('Computed {} point areas! took {}'.format(self.nverts, mtimer))
        return self.pareas

    # --------------------------------------------------------------------------
    def compute_curvatures(self):

        if self.mean_curv.shape != (0,):
            return (self.mean_curv, self.gaus_curv)

        LOGGER.info('Computing curvatures')
        mtimer = Timer()

        rval = self.tmesh.need_curvature(self.cverbose)
        if len(rval) != 2*self.nverts:
            raise ValueError('Incorrect curvature!')

        self.mean_curv = np.zeros(self.nverts,)
        self.gaus_curv = np.zeros(self.nverts,)

        for i in xrange(self.nverts):
            self.mean_curv[i] = rval[i]
            self.gaus_curv[i] = rval[self.nverts + i]

        self.mean_curv = self.mean_curv.astype(np.float32)
        self.gaus_curv = self.gaus_curv.astype(np.float32)

        mtimer.end()
        LOGGER.info('Computed {} x2 curvatures! took {}'.format(self.nverts, mtimer))
        return (self.mean_curv, self.gaus_curv)

    # --------------------------------------------------------------------------
    def compute_distance_to_surface(self, other):
        d = self.tmesh.distance_to_other_mesh(other.tmesh)
        return np.asarray(d, dtype=np.float32)

    # --------------------------------------------------------------------------
    def compute_density(self, type, sigma, name, normalize, pidxs):

        if type < 1 or type > 3:
            raise InvalidArgument('Invalid density type, {}. Should be 1 (geodesic), 2 (2D) or 3 (3D)'.format(type))

        cnt = len(self.vertices) if len(pidxs) == 0 else len(pidxs)
        tag = 'all ({})'.format(cnt) if cnt == 0 else cnt

        LOGGER.info('Estimating density of {} points [name = {}]'.format(tag, name))
        mtimer = Timer()

        k = pymemsurfer.DensityKernel(float(sigma))
        d = self.tmesh.kde(type, k, name, pidxs.tolist(), self.cverbose)
        d = np.asarray(d, dtype=np.float32)

        if normalize:
            d *= (1.0/float(cnt))

        mtimer.end()
        LOGGER.info('Computed density! took {}'.format(mtimer))
        return d

    # --------------------------------------------------------------------------
    def display(self):
        LOGGER.info('<{}> : {} verts, {} faces. bbox = {}, {}'
                    .format(self.label, self.nverts, self.nfaces,
                    self.vertices.min(axis=0), self.vertices.max(axis=0)))

    # --------------------------------------------------------------------------
    def write_vtp(self, filename, properties={}):

        duplicate_verts = False

        if not self.periodic or not duplicate_verts:
            verts = self.vertices
            properties['faces'] = self.faces

            if self.pnormals.shape != (0,0):
                properties['pnormals'] = self.pnormals
            if self.pareas.shape != (0,):
                properties['pareas'] = self.pareas
            if self.mean_curv.shape != (0,):
                properties['mean_curv'] = self.mean_curv
            if self.gaus_curv.shape != (0,):
                properties['gaus_curv'] = self.gaus_curv

        else:
            def append(a,b):
                return np.concatenate((a,b), axis=0)

            def append_dups(a,dids):
                return append(a, a[dids])

            dids = self.tmesh.duplicate_ids()
            dids = np.array(dids, dtype=np.uint32)

            verts = append(self.vertices, self.dverts)
            properties['faces'] = append(self.faces, self.tfaces)

            if self.pnormals.shape != (0,0):
                properties['pnormals'] = append_dups(self.pnormals, dids)
            if self.pareas.shape != (0,):
                properties['pareas'] = append_dups(self.pareas, dids)
            if self.mean_curv.shape != (0,):
                properties['mean_curv'] = append_dups(self.mean_curv, dids)
            if self.gaus_curv.shape != (0,):
                properties['gaus_curv'] = append_dups(self.gaus_curv, dids)

        from utils import write2vtkpolydata
        write2vtkpolydata(filename, verts, properties)

    # --------------------------------------------------------------------------
    def write_off(self, filename):

        from utils import write_off

        verts = self.vertices
        faces = self.faces

        if self.periodic:
            faces = np.vstack((faces, self.pfaces))
            #verts = np.vstack((verts, self.dverts[:,0:verts.shape[1]]))
            #faces = np.vstack((faces, self.tfaces))

        write_off(filename, verts, faces)

    # --------------------------------------------------------------------------
    def write_ply(self, filename):

        from utils import write_off

        verts = self.vertices
        faces = self.faces

        if self.periodic:
            #verts = np.vstack((verts, self.dverts[:,0:verts.shape[1]]))
            faces = np.vstack((faces, self.tfaces))

        write_ply(filename, self.vertices, self.faces)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
