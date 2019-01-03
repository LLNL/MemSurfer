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
# This file provides some simple utilities used by MemSurfer
# ------------------------------------------------------------------------------

import numpy as np
import timeit

# ------------------------------------------------------------------------------
# timing utils
# ------------------------------------------------------------------------------
class Timer(object):

    def __init__(self):
        self.start()

    def start(self):
        self.stime = timeit.default_timer()

    def end(self, print_time=False):

        self.etime = timeit.default_timer()
        if print_time:
            print self

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
# logging utils
# ------------------------------------------------------------------------------
import inspect, logging
LOGGER = logging.getLogger(__name__)

# Configuration globals
#LFORMAT = '%(asctime)s %(name)s:%(funcName)s:%(lineno)s - %(levelname)s - %(message)s'
#LFORMAT = '%(name)s:%(funcName)s: - %(message)s'
LFORMAT = '%(asctime)s [%(levelname)s] %(name)s:%(funcName)s - %(message)s'
ACCEPTED_INPUT = set(["yes", "y"])

# Create an interface to initialize a logger
def create_logger(loglevel, logstdout, logfile, logpath, logfname):
    '''
        logstdout -- enable logging to stdout!
        logfile   -- enable logging to file!
        logfname  -- logfile name!
        logpath   -- logpath name!
        loglevel  -- log level!
    '''
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
        LOGGER.debug('file {} Enabled'.format(logfname))

    return LOGGER

    # Print the level of logging.
    LOGGER.debug('Enabled')
    LOGGER.info('Enabled')
    LOGGER.warning('Enabled')
    LOGGER.error('Enabled')
    LOGGER.critical('Enabled')
    return LOGGER

# ------------------------------------------------------------------------------
# i/o
# ------------------------------------------------------------------------------
def read_off(filename):

    verts = []
    faces = []
    with open(filename,'r') as file:

        s = file.readline()
        s = file.readline().split()

        nv = int(s[0])
        nf = int(s[1])

        s = file.readlines()
        s = [l.split() for l in s]

        verts = [[float(l[0]), float(l[1]), float(l[2])] for l in s[:nv]]
        faces = [[int(l[1]), int(l[2]), int(l[3])] for l in s[nv:nv+nf]]
        file.close()

    return verts, faces

def write_off(filename, verts, faces):

    nv = len(verts)
    nf = len(faces)

    with open(filename,'wb') as file:
        file.write('OFF\n')
        file.write("{} {} {}\n".format(nv, nf, 0))
        for p in verts:
            if(len(p) == 2):
                file.write("{0:0.6f} {1:0.6f} 0.0\n".format(p[0],p[1]))
            else:
                file.write("{0:0.6f} {1:0.6f} {2:0.6f}\n".format(p[0],p[1], p[2]))

        for f in faces:
            number = len(f)
            row = "{0}".format(number)
            for elem in f:
                row += " {0} ".format(elem)
            row += "\n"
            file.write(row)

def read_ply(filename):
    pass

# write a surface as a ply file
def write_ply(filename, verts, faces):

    nv = len(verts)
    nf = len(faces)

    header = '''ply
format ascii 1.0
element vertex {0}
property float x
property float y
property float z
element face {1}
property list uchar int vertex_indices
end_header\n'''.format(nv, nf)

    with open(filename,'wb') as file:
        file.writelines(header)
        for p in verts:
            if(len(p) == 2):
                file.write("{0:0.6f} {1:0.6f} 0.0\n".format(p[0],p[1]))
            else:
                file.write("{0:0.6f} {1:0.6f} {2:0.6f}\n".format(p[0],p[1], p[2]))

        for f in faces:
            number = len(f)
            row = "{0}".format(number)
            for elem in f:
                row += " {0} ".format(elem)
            row += "\n"
            file.write(row)

# ------------------------------------------------------------------------------
# vtk i/o
# ------------------------------------------------------------------------------
import vtk
from vtk.util import numpy_support

def write2vtkpolydata(filename, verts, properties):

    LOGGER.info('Creating vtkdata and writing to [{}]'.format(filename))

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(filename)

    polydata = vtk.vtkPolyData()

    points = vtk.vtkPoints()
    cells = vtk.vtkCellArray()

    # --------------------------------------------------------------------------
    for v in verts:
        if len(v) == 3:     points.InsertNextPoint(v[0], v[1], v[2])
        else:               points.InsertNextPoint(v[0], v[1], 0)
    polydata.SetPoints(points)

    # --------------------------------------------------------------------------
    if 'faces' in properties.keys():
        faces = properties['faces']
        for f in faces:
            cell = vtk.vtkTriangle()
            for i in xrange(3):
                cell.GetPointIds().SetId(i, f[i])
            cells.InsertNextCell(cell)
        polydata.SetPolys(cells)

    # --------------------------------------------------------------------------
    for key in properties.keys():
        if key == 'faces':
            continue

        data = properties[key]
        if not isinstance(data, np.ndarray):
            data = np.array([data])

        VTK_data = numpy_support.numpy_to_vtk(num_array=data)
        VTK_data.SetName(key)

        if data.shape[0] == 1:
            polydata.GetFieldData().AddArray(VTK_data)

        elif data.shape[0] == points.GetNumberOfPoints():
            polydata.GetPointData().AddArray(VTK_data)

        elif data.shape[0] == cells.GetNumberOfCells():
            polydata.GetCellData().AddArray(VTK_data)

    # --------------------------------------------------------------------------
    if 'faces' not in properties.keys():
        for i in xrange(len(verts)):
            cell = vtk.vtkVertex()
            cell.GetPointIds().SetId(0, i)
            cells.InsertNextCell(cell)
        polydata.SetPolys(cells)


    writer.SetInputData(polydata)
    writer.Write()
    LOGGER.info('File [{}] successfully written.'.format(filename))

# ------------------------------------------------------------------------------
