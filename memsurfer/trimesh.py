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
import numpy as np
import vtk
from vtk.util import numpy_support

import logging
LOGGER = logging.getLogger(__name__)

from . import memsurfer_cmod
from .utils import Timer


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
class TriMesh(object):
    """
       Class to create and manipulate triangular meshes (in 2D and 3D)
    """

    # --------------------------------------------------------------------------
    # constructor
    def __init__(self, vertices, **kwargs):
        """
        vertices: ndarray of shape (nverts, dim) # dim = 2,3
        kwargs:
                  faces: ndarray of shape (nfaces,3)
                  label: label for this mesh
                  periodic: whether this mesh is periodic
        """
        if vertices.shape[1] not in [2,3]:
            raise ValueError(f'TriMesh needs 2D or 3D vertices: got {vertices.shape}')

        self.vertices = vertices.astype(np.float32)
        self.tmesh = memsurfer_cmod.TriMesh(self.vertices)

        if 'faces' in list(kwargs.keys()):
            self.faces = kwargs['faces']
            self.faces = self.faces.astype(np.uint32)
            if self.faces.shape[1] != 3:
                raise ValueError(f'TriMesh needs triangles: got {faces.shape}')
            self.tmesh.set_faces(self.faces)

        else:
            self.faces = None

        self.periodic = kwargs.get('periodic', False)
        if self.periodic:
            self.tmesh.set_periodic()
            self.label = kwargs.get('label', 'TriMeshPeriodic')
        else:
            self.label = kwargs.get('label', 'TriMesh')

        nfaces = self.faces.shape[0] if self.faces is not None else 0
        LOGGER.info(f'{self.tag()} Created {self.vertices.shape} vertices and {nfaces} faces')
        self.cverbose = LOGGER.isEnabledFor(logging.DEBUG)

        # other parameters that will be created later!
        self.bbox = None            # bounding box

        self.pnormals = None        # point normals
        self.pareas = None          # point areas
        self.mean_curv = None       # mean curvature
        self.gaus_curv = None       # gaussian curvature
        self.shells = None          # shell index

        self.pverts = None          # parameterized vertices
        self.pfaces = None          # periodic faces
        self.tfaces = None          # trimmed faces
        self.dverts = None          # duplicated verts

    # --------------------------------------------------------------------------
    def tag(self):
        return f'[{self.label}, periodic={self.periodic}]'

    def __str__(self):
        return f'{self.tag()}: ' \
               f'{self.vertices.shape[0]} verts, ' \
               f'{self.faces.shape[0]} faces. ' \
               f'bbox = {self.vertices.min(axis=0)}, {self.vertices.max(axis=0)}'

    def __repr__(self):
        return self.__str__()

    # --------------------------------------------------------------------------
    def set_bbox(self, bb0, bb1):

        if bb0.shape[0] == 2 and bb1.shape[0] == 2:
            self.bboxl = np.array([bb0[0], bb0[1],  bb1[0], bb1[1]])
            self.bbox = np.array([[bb0[0], bb0[1], 0.], [bb1[0], bb1[1], 0.]])

        elif bb0.shape[0] == 3 and bb1.shape[0] == 3:
            self.bboxl = np.array([bb0[0], bb0[1], bb0[2], bb1[0], bb1[1], bb1[2]])
            self.bbox = np.array([[bb0[0], bb0[1], bb0[2]], [bb1[0], bb1[1], bb1[2]]])

        else:
            raise ValueError(f'TriMesh got incorrect bbox: {bb0}, {bb1}')

        self.bbox = self.bbox.astype(np.float32)
        LOGGER.info(f'{self.tag()} Setting bbox = {self.bbox[0]}, {self.bbox[1]}')
        self.tmesh.set_bbox(self.bboxl)

    # --------------------------------------------------------------------------
    def as_vtkpolydata(self, props_to_include = [], additional_props = {}):

        LOGGER.info(f'{self.tag()} Converting to vtk '
                    f'(include={props_to_include}, '
                    f'additional={list(additional_props.keys())})')

        polydata = vtk.vtkPolyData()
        points = vtk.vtkPoints()
        cells = vtk.vtkCellArray()

        # ----------------------------------------------------------------------
        if not self.periodic:
            verts = self.vertices
            faces = self.faces
            dids = np.array([], dtype=np.uint32)
        else:
            verts = np.concatenate((self.vertices, self.dverts), axis=0)
            faces = np.concatenate((self.faces, self.tfaces), axis=0)
            dids = np.array(self.tmesh.duplicate_ids(), dtype=np.uint32)

        # ----------------------------------------------------------------------
        # convert points
        if verts.shape[1] == 2:
            for v in verts:
                points.InsertNextPoint(v[0], v[1], 0.)
        else:
            for v in verts:
                points.InsertNextPoint(v[0], v[1], v[2])
        polydata.SetPoints(points)

        # ----------------------------------------------------------------------
        # if faces are available, create polygons
        if faces is not None:
            if not self.periodic:
                for f in faces:
                    cell = vtk.vtkTriangle()
                    for i in range(3):
                        cell.GetPointIds().SetId(i, f[i])
                    cells.InsertNextCell(cell)

            else:
                xw = 0.5 * (self.bbox[1,0] - self.bbox[0,0])
                yw = 0.5 * (self.bbox[1,1] - self.bbox[0,1])

                for f in faces:
                    is_periodic = False
                    cell = vtk.vtkTriangle()
                    for i in range(3):
                        cell.GetPointIds().SetId(i, f[i])

                        p = polydata.GetPoint(f[i])
                        q = polydata.GetPoint(f[(i + 1) % 3])
                        is_periodic |= abs(p[0]-q[0]) > xw or \
                                       abs(p[1]-q[1]) > yw

                    if not is_periodic:
                        cells.InsertNextCell(cell)
            polydata.SetPolys(cells)

        # ----------------------------------------------------------------------
        # if faces are not available, each point is a vtk vertex
        else:
            for i in range(self.vertices.shape[0]):
                cell = vtk.vtkVertex()
                cell.GetPointIds().SetId(0, i)
                cells.InsertNextCell(cell)
            polydata.SetVerts(cells)

        # ----------------------------------------------------------------------
        def _add_to_polydata(_name, _data):

            # everything must be a numpy array
            if not isinstance(_data, np.ndarray):
                _data = np.array([_data])

            # duplicate if needed, if this is point data
            if _data.shape[0] == self.vertices.shape[0] and dids.shape[0] > 0:
                _data = np.concatenate((_data, _data[dids]), axis=0)

            if _data.dtype.kind not in ['U', 'S']:
                _vtkdata = numpy_support.numpy_to_vtk(num_array=_data)

            else:
                _vtkdata = vtk.vtkStringArray()
                _vtkdata.SetNumberOfValues(_data.shape[0])
                for i, v in enumerate(_data):
                    _vtkdata.SetValue(i, str(v))

            _vtkdata.SetName(_name)

            # now, where does this need to go?
            if _data.shape[0] == 1:
                polydata.GetFieldData().AddArray(_vtkdata)

            elif _data.shape[0] == points.GetNumberOfPoints():
                polydata.GetPointData().AddArray(_vtkdata)

            elif _data.shape[0] == cells.GetNumberOfCells():
                polydata.GetCellData().AddArray(_vtkdata)

            else:
                LOGGER.warning(f'Could not add {_name} with size {_data.shape}')

        # ----------------------------------------------------------------------
        if 'pnormals' in props_to_include and self.pnormals is not None:
            _add_to_polydata('pnormals', self.pnormals)

        if 'pareas' in props_to_include and self.pareas is not None:
            _add_to_polydata('pareas', self.pareas)

        if 'mean_curv' in props_to_include and self.mean_curv is not None:
            _add_to_polydata('mean_curv', self.mean_curv)

        if 'gaus_curv' in props_to_include and self.gaus_curv is not None:
            _add_to_polydata('gaus_curv', self.gaus_curv)

        if 'shells' in props_to_include and self.shells is not None:
            _add_to_polydata('shell_idx', self.shells)

        for k, v in additional_props.items():
            _add_to_polydata(k, v)

        # ----------------------------------------------------------------------
        LOGGER.info(f'{self.tag()} Converted to vtk polydata')
        return polydata

    # --------------------------------------------------------------------------
    # mesh manupulation
    # --------------------------------------------------------------------------
    def parameterize(self, xy=False):

        if self.pverts is not None:
            return self.pverts

        LOGGER.info(f'{self.tag()} Parameterizing the surface')
        mtimer = Timer()

        if xy:
            pverts = self.tmesh.parameterize_xy(self.cverbose)
        else:
            pverts = self.tmesh.parameterize(self.cverbose)

        self.pverts = np.array(pverts).reshape(-1, 2).astype(np.float32)
        assert self.pverts.shape[0] == self.vertices.shape[0]

        mtimer.end()
        LOGGER.info(f'{self.tag()} Parameterization took {mtimer}')
        return self.pverts

    # --------------------------------------------------------------------------
    def project_on_surface_and_plane(self, points):
        """
        Project a set of points (given as barycentric coordinates and face ids)
            onto the triangulation and the parameterized triangulation
        """
        self.parameterize()

        npoints = points.shape[0]
        LOGGER.info(f'{self.tag()} Projecting {npoints} points on the parameterized surface')

        mtimer = Timer()

        r = self.tmesh.project_on_surface(points.reshape(-1).tolist(), self.cverbose)
        r = np.array(r).reshape(npoints, 4)

        fids = r[:,0].astype(int)
        fbarys = r[:,1:4]

        # create points using the projections
        spoints = np.zeros((npoints, 3), dtype=np.float32)
        ppoints = np.zeros((npoints, 2), dtype=np.float32)

        for i in range(npoints):

            f = self.faces[fids[i]]
            for j in range(3):
                spoints[i] += fbarys[i][j] * self.vertices[f[j]]
                ppoints[i] += fbarys[i][j] * self.pverts[f[j]]

        mtimer.end()
        LOGGER.info(f'{self.tag()} Projecting took {mtimer}')
        LOGGER.info(f'\tbbox of planar projections: {ppoints.min(axis=0)}, {ppoints.max(axis=0)}')
        LOGGER.info(f'\tbbox of surface projections: {spoints.min(axis=0)}, {spoints.max(axis=0)}')
        return spoints, ppoints

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

        if self.faces is not None:
            return

        tag = ' periodic ' if self.periodic else ' '
        LOGGER.info(f'{self.tag()} Computing 2D{tag}Delaunay triangulation')
        mtimer = Timer()

        faces = self.tmesh.delaunay(self.cverbose)
        self.faces = np.array(faces).reshape(-1, 3).astype(np.uint32)

        if not self.periodic:
            mtimer.end()
            LOGGER.info(f'{self.tag()} Delaunay triangulation took {mtimer}! '
                        f'created {self.faces.shape[0]} faces and {self.vertices.shape[0]} vertices')
            return

        # handle duplicated elements
        pfaces = self.tmesh.periodic_faces()
        self.pfaces = np.array(pfaces).reshape(-1, 3).astype(np.uint32)

        tfaces = self.tmesh.trimmed_faces()
        self.tfaces = np.array(tfaces).reshape(-1, 3).astype(np.uint32)

        dverts = self.tmesh.duplicated_vertices()
        self.dverts = np.array(dverts).reshape(-1, self.vertices.shape[1]).astype(np.float32)

        mtimer.end()
        LOGGER.info(f'{self.tag()} Delaunay triangulation took {mtimer}! '
                    f'created [{self.faces.shape[0]} {self.pfaces.shape[0]} {self.tfaces.shape[0]}] faces, '
                    f'and [{self.vertices.shape[0]} {self.dverts.shape[0]}] vertices')

    # --------------------------------------------------------------------------
    # computation of properties!
    # --------------------------------------------------------------------------
    def compute_normals(self):

        if self.pnormals is not None:
            return self.pnormals

        LOGGER.info(f'{self.tag()} Computing normals')
        mtimer = Timer()

        rval = self.tmesh.need_normals(self.cverbose)
        if len(rval) != 3*self.vertices.shape[0]:
            raise ValueError('Incorrect normals!')

        self.pnormals = np.asarray(rval).reshape(-1,3).astype(np.float32)

        mtimer.end()
        LOGGER.info(f'{self.tag()} Computed {self.pnormals.shape[0]} normals! took {mtimer}')
        return self.pnormals

    # --------------------------------------------------------------------------
    def compute_pointareas(self):

        if self.pareas is not None:
            return self.pareas

        LOGGER.info(f'{self.tag()} Computing point areas')
        mtimer = Timer()

        rval = self.tmesh.need_pointareas(self.cverbose)
        if len(rval) != self.vertices.shape[0]:
            raise ValueError('Incorrect point areas!')

        self.pareas = np.array(rval).astype(np.float32)

        mtimer.end()
        LOGGER.info(f'{self.tag()} Computed {self.pareas.shape[0]} areas! took {mtimer}')
        return self.pareas

    # --------------------------------------------------------------------------
    def compute_curvatures(self):

        if self.mean_curv is not None and self.gaus_curv is not None:
            return self.mean_curv, self.gaus_curv

        LOGGER.info(f'{self.tag()} Computing curvature')
        mtimer = Timer()

        # convert into vtk polydata
        polydata = self.as_vtkpolydata()

        # use vtk to compute curvature
        mc = vtk.vtkCurvatures()
        mc.SetInputData(polydata)
        mc.SetCurvatureTypeToMean()
        mc.Update()

        gc = vtk.vtkCurvatures()
        gc.SetInputData(polydata)
        gc.SetCurvatureTypeToGaussian()
        gc.Update()

        self.mean_curv = numpy_support.vtk_to_numpy(mc.GetOutput().GetPointData().GetArray('Mean_Curvature'))
        self.gaus_curv = numpy_support.vtk_to_numpy(gc.GetOutput().GetPointData().GetArray('Gauss_Curvature'))

        mtimer.end()
        LOGGER.info(f'{self.tag()} Computed {self.vertices.shape[0]} x2 curvatures! took {mtimer}')
        return self.mean_curv, self.gaus_curv

    # --------------------------------------------------------------------------
    def compute_distance_to_surface(self, other):
        d = self.tmesh.distance_to_other_mesh(other.tmesh)
        return np.asarray(d, dtype=np.float32)

    # --------------------------------------------------------------------------
    def compute_density(self, type, sigma, name, get_nlipids, pidxs):

        if type < 1 or type > 3:
            raise InvalidArgument('Invalid density type, {}. Should be 1 (geodesic), 2 (2D) or 3 (3D)'.format(type))

        cnt = f'all (of {self.vertices.shape[0]})' if len(pidxs) == 0 else \
              f'{len(pidxs)} (of {self.vertices.shape[0]})'

        LOGGER.info(f'{self.tag()} Estimating density of {cnt} points [name = {name}]')
        mtimer = Timer()

        # ----------------------------------------------------------------------
        # based on the type of density, choose the correct kernel!
        if type == 1 or type == 2:
            dens_kern = memsurfer_cmod.GaussianKernel2D(float(sigma))
        elif type == 3:
            dens_kern = memsurfer_cmod.GaussianKernel3D(float(sigma))
        else:
            assert False

        # ----------------------------------------------------------------------
        # based on periodicity, choose the correct distance kernel!
        if self.periodic:
            dist_kern = memsurfer_cmod.DistancePeriodicXYSquared(self.bboxl)
        else:
            dist_kern = memsurfer_cmod.DistanceSquared()

        # ----------------------------------------------------------------------
        d = self.tmesh.kde(name, type, get_nlipids, dens_kern, dist_kern,
                            pidxs.tolist(), self.cverbose)
        d = np.asarray(d, dtype=np.float32)

        mtimer.end()
        LOGGER.info(f'{self.tag()} Computed density! took {mtimer}')
        LOGGER.info(f'\trange = [{d.min()}, {d.max()}], sum = {d.sum()}')

        # ----------------------------------------------------------------------
        return d

    # --------------------------------------------------------------------------
    def compute_shells(self, ref_pts):

        assert isinstance(ref_pts, (list, np.ndarray))
        assert len(ref_pts) > 0

        # compute shell index with respect to the reference points!
        LOGGER.info(f'{self.tag()} Computing shells with respect to {len(ref_pts)} vertices')
        mtimer = Timer()

        nverts = self.vertices.shape[0]

        # get neighborhood graph
        ret_val = self.tmesh.need_neighbors()
        ret_val = np.array(ret_val).astype(np.int)

        # reformat the received data to create neighbor list
        nnbrs = ret_val[:nverts]        # first n are the counts

        nbrs = [[] for i in range(nverts)]
        nidx = nverts
        for vidx in range(nverts):
            nbrs[vidx] = ret_val[nidx:nidx+nnbrs[vidx]]
            nidx += nnbrs[vidx]

        # now, populate the shells using a simple bfs on this graph
        self.shells = -1 * np.ones(nverts)

        # initialize the bfs with reference points as shell = 0
        for p in ref_pts:
            self.shells[p] = 0

        vertices_to_process = list(ref_pts)
        while len(vertices_to_process) > 0:
            current_vertex = vertices_to_process.pop(0)
            for nbr in nbrs[current_vertex]:
                if self.shells[nbr] == -1:   # this vertex has not been visited yet!
                    self.shells[nbr] = 1+self.shells[current_vertex]
                    vertices_to_process.append(nbr)

        assert self.shells.min() == 0, 'Failed to assign shell id for some vertices'
        mtimer.end()
        LOGGER.info(f'{self.tag()} Computed shells for {self.vertices.shape[0]} vertices! took {mtimer}')

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    def copy_triangulation(self, mesh):

        if self.periodic != mesh.periodic:
            raise ValueError(f'Cannot copy a (periodic={self.periodic}) mesh into ({mesh.periodic})')

        self.tmesh.copy_periodicDelaunay(mesh.tmesh)
        self.faces = self.tmesh.get_faces()
        self.faces = np.array(self.faces).reshape(-1, 3).astype(np.uint32)

        if self.periodic:
            self.pfaces = self.tmesh.periodic_faces()
            self.tfaces = self.tmesh.trimmed_faces()
            self.dverts = self.tmesh.duplicated_vertices()

            self.pfaces = np.array(self.pfaces).reshape(-1, 3).astype(np.uint32)
            self.tfaces = np.array(self.tfaces).reshape(-1, 3).astype(np.uint32)
            self.dverts = np.array(self.dverts).reshape(-1, self.vertices.shape[1]).astype(np.float32)

    def copy_densities(self, mesh):
        self.tmesh.set_fields(mesh.tmesh, 'density')

    # --------------------------------------------------------------------------
    def write_vtp(self, filename, additional_props={}):

        # 4x properties are stored in the python object
        props_to_include = ['pnormals', 'pareas', 'mean_curv', 'gaus_curv', 'shells']

        polydata = self.as_vtkpolydata(props_to_include, additional_props)

        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(filename)
        writer.SetInputData(polydata)
        writer.Write()
        LOGGER.info(f'{self.tag()} Written to ({filename})')

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
