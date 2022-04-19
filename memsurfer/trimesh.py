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

    KNOWN_ATTRIBUTES = ['pnormal', 'parea', 'mean_curv', 'gaus_curv', 'shell']

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

        # geometric properties that will be created later!
        self.bbox = None       # bounding box
        self.bboxl = None      # same as bbox, but linearized

        self._periodic_faces = None         # periodic faces
        self._trimmed_faces = None          # trimmed faces
        self._duplicated_vidxs = None       # duplicated verts
        self._duplicated_verts = None       # duplicated verts
        self._parameterized_verts = None    # parameterized vertices

        self.attributes = {}

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
    # copy info from outside
    # --------------------------------------------------------------------------
    def copy_triangulation(self, mesh):

        if self.periodic != mesh.periodic:
            raise ValueError(f'Cannot copy a (periodic={self.periodic}) mesh into ({mesh.periodic})')

        self.tmesh.copy_periodicDelaunay(mesh.tmesh)
        self._fetch_faces()

    def copy_densities(self, mesh):
        self.tmesh.set_fields(mesh.tmesh, 'density')

    # --------------------------------------------------------------------------
    # mesh manipulation
    # --------------------------------------------------------------------------
    def parameterization_bbox(self):
        bb0 = self._parameterized_verts.min(axis=0)
        bb1 = self._parameterized_verts.max(axis=0)
        return bb0, bb1

    def parameterize(self, xy=False):

        if self._parameterized_verts is not None:
            return self._parameterized_verts

        LOGGER.info(f'{self.tag()} Parameterizing the surface')
        mtimer = Timer()

        if xy:
            pverts = self.tmesh.parameterize_xy(self.cverbose)
        else:
            pverts = self.tmesh.parameterize(self.cverbose)

        self._parameterized_verts = np.array(pverts, dtype=np.float32)\
            .reshape(-1, 2)
        assert self._parameterized_verts.shape[0] == self.vertices.shape[0]

        mtimer.end()
        LOGGER.info(f'{self.tag()} Parameterization took {mtimer}')
        return self._parameterized_verts

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
                ppoints[i] += fbarys[i][j] * self._parameterized_verts[f[j]]

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
    def compute_delaunay(self):

        if self.faces is not None:
            return

        tag = ' periodic ' if self.periodic else ' '
        LOGGER.info(f'{self.tag()} Computing 2D{tag}Delaunay triangulation')
        mtimer = Timer()

        self.tmesh.delaunay(self.cverbose)
        self._fetch_faces()
        mtimer.end()

        if not self.periodic:
            LOGGER.info(f'{self.tag()} Delaunay triangulation took {mtimer}! '
                        f'created {self.faces.shape[0]} faces and {self.vertices.shape[0]} vertices')
        else:
            LOGGER.info(f'{self.tag()} Delaunay triangulation took {mtimer}! '
                        f'created [{self.faces.shape[0]} {self._periodic_faces.shape[0]} {self._trimmed_faces.shape[0]}] faces '
                        f'and [{self.vertices.shape[0]} {self._duplicated_verts.shape[0]}] vertices')

        # check vertex incidence
        vertex_degree = np.zeros(self.vertices.shape[0], dtype=int)
        for f in self.faces:
            vertex_degree[f] += 1

        incorrect = np.where(vertex_degree < 2)[0]
        if incorrect.shape[0] > 0:
            LOGGER.warning(f'Found vertices with small incidence: {incorrect} = {vertex_degree[incorrect]}')
            print (self.vertices[incorrect])

    # --------------------------------------------------------------------------
    # computation of attributes!
    # --------------------------------------------------------------------------
    def compute_normals(self):

        key = 'pnormal'
        if key in self.attributes:
            return self.attributes[key]

        LOGGER.info(f'{self.tag()} Computing vertex normals')
        mtimer = Timer()

        rval = self.tmesh.need_normals(self.cverbose)
        if len(rval) != 3*self.vertices.shape[0]:
            raise ValueError('Incorrect vertex normals!')

        self.attributes[key] = np.asarray(rval).reshape(-1,3).astype(np.float32)

        mtimer.end()
        LOGGER.info(f"{self.tag()} Computed {self.attributes[key].shape[0]} normals! took {mtimer}")
        return self.attributes[key]

    # --------------------------------------------------------------------------
    def compute_pointareas(self):

        key = 'parea'
        if key in self.attributes:
            return self.attributes[key]

        LOGGER.info(f'{self.tag()} Computing vertex areas')
        mtimer = Timer()

        rval = self.tmesh.need_pointareas(self.cverbose)
        if len(rval) != self.vertices.shape[0]:
            raise ValueError('Incorrect vertex areas!')

        self.attributes[key] = np.array(rval).astype(np.float32)

        mtimer.end()
        LOGGER.info(f"{self.tag()} Computed {self.attributes[key].shape[0]} areas! took {mtimer}")
        return self.attributes[key]

    # --------------------------------------------------------------------------
    def compute_curvatures(self):

        keys = ['mean_curv', 'gaus_curv']
        if keys[0] in self.attributes and keys[1] in self.attributes:
            return self.attributes[keys[0]], self.attributes[keys[1]]

        # convert into vtk polydata
        polydata = self.as_vtkpolydata(include_attributes=[])

        LOGGER.info(f'{self.tag()} Computing curvature')
        mtimer = Timer()

        # use vtk to compute curvature
        mc = vtk.vtkCurvatures()
        mc.SetInputData(polydata)
        mc.SetCurvatureTypeToMean()
        mc.Update()

        gc = vtk.vtkCurvatures()
        gc.SetInputData(polydata)
        gc.SetCurvatureTypeToGaussian()
        gc.Update()

        m = numpy_support.vtk_to_numpy(mc.GetOutput().GetPointData().GetArray('Mean_Curvature'))
        g = numpy_support.vtk_to_numpy(gc.GetOutput().GetPointData().GetArray('Gauss_Curvature'))

        nverts = self.vertices.shape[0]
        self.attributes[keys[0]] = m[:nverts]
        self.attributes[keys[1]] = g[:nverts]

        mtimer.end()
        LOGGER.info(f'{self.tag()} Computed {nverts} x2 curvatures! took {mtimer}')
        return self.attributes[keys[0]], self.attributes[keys[1]]

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

        key = 'shell'
        if key in self.attributes:
            return self.attributes[key]

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
        shells = -1 * np.ones(nverts, dtype=int)

        # initialize the bfs with reference points as shell = 0
        for p in ref_pts:
            shells[p] = 0

        vertices_to_process = list(ref_pts)
        while len(vertices_to_process) > 0:
            current_vertex = vertices_to_process.pop(0)
            for nbr in nbrs[current_vertex]:
                if shells[nbr] == -1:   # this vertex has not been visited yet!
                    shells[nbr] = 1+shells[current_vertex]
                    vertices_to_process.append(nbr)

        if shells.min() < 0:
            wrong = np.where(shells < 0)[0]
            LOGGER.error(f'Failed to assign shell id for {wrong.shape[0]} vertices: '
                         f'ids = {wrong}, verts = {self.vertices[wrong]}')
            assert shells.min() == 0, 'Failed to assign shell id for some vertices'

        self.attributes[key] = shells.astype(int)
        mtimer.end()
        LOGGER.info(f'{self.tag()} Computed shells for {self.vertices.shape[0]} vertices! took {mtimer} :: {shells.min()}, {shells.max()}')
        return self.attributes[key]

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # fetch faces from the trimesh (cpp) object
    def _fetch_faces(self):

        self.faces = np.array(self.tmesh.get_faces(), dtype=np.uint32) \
            .reshape(-1, 3)

        if not self.periodic:
            return

        self._periodic_faces = np.array(self.tmesh.periodic_faces(), dtype=np.uint32) \
            .reshape(-1, 3)
        self._trimmed_faces = np.array(self.tmesh.trimmed_faces(), dtype=np.uint32) \
            .reshape(-1, 3)

        self._duplicated_vidxs = np.array(self.tmesh.duplicate_ids(), dtype=np.uint32)
        self._duplicated_verts = np.array(self.tmesh.duplicated_vertices(), dtype=np.float32) \
            .reshape(-1, self.vertices.shape[1])

    # --------------------------------------------------------------------------
    # convert the mesh to vtk and pandas format for output
    # --------------------------------------------------------------------------
    def as_vtkpolydata(self,
                       include_attributes=KNOWN_ATTRIBUTES,
                       additional_attributes={}):

        LOGGER.info(f'{self.tag()} Converting to vtk polydata '
                    f'(include={include_attributes}, '
                    f'additional={list(additional_attributes.keys())})')

        include_periodic = self.periodic and (self._duplicated_vidxs is not None)

        polydata = vtk.vtkPolyData()
        points = vtk.vtkPoints()
        cells = vtk.vtkCellArray()

        # ----------------------------------------------------------------------
        if include_periodic:
            verts = np.concatenate((self.vertices, self._duplicated_verts), axis=0)
            faces = np.concatenate((self.faces, self._trimmed_faces), axis=0)
        else:
            verts = self.vertices
            faces = self.faces

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
            if not include_periodic:
                for f in faces:
                    cell = vtk.vtkTriangle()
                    for i in range(3):
                        cell.GetPointIds().SetId(i, f[i])
                    cells.InsertNextCell(cell)

            else:
                xw = 0.5 * (self.bbox[1, 0] - self.bbox[0, 0])
                yw = 0.5 * (self.bbox[1, 1] - self.bbox[0, 1])

                for f in faces:
                    is_periodic = False
                    cell = vtk.vtkTriangle()
                    for i in range(3):
                        cell.GetPointIds().SetId(i, f[i])

                        p = polydata.GetPoint(f[i])
                        q = polydata.GetPoint(f[(i + 1) % 3])
                        is_periodic |= abs(p[0] - q[0]) > xw or \
                                       abs(p[1] - q[1]) > yw

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
            if isinstance(_data, (int, float, str)):
                _data = np.array([_data])

            elif not isinstance(_data, np.ndarray):
                _data = np.array(_data)

            try:
                # duplicate if needed, if this is point data
                if include_periodic and _data.shape[0] == self.vertices.shape[0]:
                    _data = np.concatenate((_data, _data[self._duplicated_vidxs]), axis=0)

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
                    raise Exception('Failed to match size of data')

            except Exception as e:
                LOGGER.warning(f'Could not add "{_name}" with size {_data.shape}. '
                               f'(nverts = {points.GetNumberOfPoints()}, '
                               f'faces = {cells.GetNumberOfCells()}). Error = {e}')

        # ----------------------------------------------------------------------
        for att in include_attributes:
            if att in self.attributes.keys():
                _add_to_polydata(att, self.attributes[att])
            else:
                LOGGER.warning(f'Requested attribute "{att}" is not computed')

        for k, v in additional_attributes.items():
            _add_to_polydata(k, v)

        # ----------------------------------------------------------------------
        LOGGER.info(f'{self.tag()} Converted to vtk polydata')
        return polydata

    # --------------------------------------------------------------------------
    def as_pddataframe(self,
                       include_attributes=KNOWN_ATTRIBUTES,
                       additional_attributes={}):

        import pandas as pd
        LOGGER.info(f'{self.tag()} Converting to Pandas dataframe '
                    f'(include={include_attributes}, '
                    f'additional={list(additional_attributes.keys())})')

        # this function does not list duplicate vertices
        # as_vtkpolydata does!
        # include_periodic = False
        df = pd.DataFrame()

        # ----------------------------------------------------------------------
        def _add_to_dataframe(_name, _data):

            # everything must be a numpy array
            if not isinstance(_data, np.ndarray):
                _data = np.array([_data])

            # now add the data
            if _data.ndim == 1:
                df[_name] = _data

            elif _data.ndim == 2 and _data.shape[1] == self.vertices.shape[1]:
                dims = ['x', 'y', 'z']
                for d in range(self.vertices.shape[1]):
                    df[f'{_name}_{dims[d]}'] = _data[:, d]

            else:
                LOGGER.warning(f'Could not add "{_name}" with size {_data.shape}')

        # ----------------------------------------------------------------------
        _add_to_dataframe('pos', self.vertices)
        for att in include_attributes:
            if att in self.attributes.keys():
                _add_to_dataframe(att, self.attributes[att])
            else:
                LOGGER.warning(f'Requested attribute "{att}" is not computed')

        for k, v in additional_attributes.items():
            _add_to_dataframe(k, v)

        # ----------------------------------------------------------------------
        LOGGER.info(f'{self.tag()} Converted to Pandas dataframe')
        return df

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    def write_vtp(self, filename, additional_attrs={}):

        polydata = self.as_vtkpolydata(TriMesh.KNOWN_ATTRIBUTES, additional_attrs)

        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(filename)
        writer.SetInputData(polydata)
        writer.Write()
        LOGGER.info(f'{self.tag()} Written to ({filename})')

    def write_pd(self, filename, additional_attrs={}):

        df = self.as_pddataframe(TriMesh.KNOWN_ATTRIBUTES, additional_attrs)

        df.to_csv(filename)
        LOGGER.info(f'{self.tag()} Written to ({filename})')

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
