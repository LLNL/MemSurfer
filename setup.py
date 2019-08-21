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

import os, sys
import glob
import numpy
from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize

# ------------------------------------------------------------------------------
# Read the correct paths for the depenencies from the environment
# ------------------------------------------------------------------------------
def read_paths(default):

    try:
        VTK_ROOT = os.path.expandvars(os.environ['VTK_ROOT'])
    except KeyError:
        VTK_ROOT = default

    try:
        CGAL_ROOT = os.path.expandvars(os.environ['CGAL_ROOT'])
    except KeyError:
        CGAL_ROOT = default

    try:
        EIGEN_ROOT = os.path.expandvars(os.environ['EIGEN_ROOT'])
    except KeyError:
        EIGEN_ROOT = default

    try:
        BOOST_ROOT = os.path.expandvars(os.environ['BOOST_ROOT'])
    except KeyError:
        BOOST_ROOT = default

    return VTK_ROOT, CGAL_ROOT, EIGEN_ROOT, BOOST_ROOT

def choose_path(root, libname):

    options = ['lib', 'lib64']
    extns = ['so', 'dylib']

    for o in options:
      for e in extns:
        if os.path.isfile(os.path.join(root, o, '{}.{}'.format(libname, e))):
          return os.path.join(root, o)

    raise Exception('choose_path({},{}) failed!'.format(root, libname))

# ------------------------------------------------------------------------------
# Read the correct paths for the code, including the dependencies
# ------------------------------------------------------------------------------
SRC_PATH = os.path.split(os.path.abspath(sys.argv[0]))[0]
VTK_ROOT, CGAL_ROOT, EIGEN_ROOT, BOOST_ROOT = read_paths(os.path.join(SRC_PATH, 'external'))

print('--- vtk', VTK_ROOT)
print('--- cgal', CGAL_ROOT)
print('--- boost', BOOST_ROOT)
print('--- eigen', EIGEN_ROOT)

# cgal gets installed in lib or lib64
CGAL_lpath = choose_path(CGAL_ROOT, 'libCGAL')
CGAL_ipath = os.path.join(CGAL_ROOT, 'include')

VTK_lpath = os.path.join(VTK_ROOT, 'lib')
VTK_ipath = os.path.join(VTK_ROOT, 'include', 'vtk-8.1')

EIG_ipath = os.path.join(EIGEN_ROOT, 'include', 'eigen3')
BOOST_ipath = os.path.join(BOOST_ROOT, 'include')

# ------------------------------------------------------------------------------
# Collect the names of libraries to link to and code to compile
# ------------------------------------------------------------------------------

# cpp code
headers = os.path.join(SRC_PATH, 'libsurfer', 'include')
sources = glob.glob(os.path.join(SRC_PATH, 'libsurfer', 'src', '*.cpp'))
sources.append(os.path.join(SRC_PATH, 'memsurfer', 'pymemsurfer.i'))

# external libs
libs = glob.glob(os.path.join(VTK_lpath, 'libvtk*'))
libs = [l for l in libs if os.path.islink(l)]
libs = [os.path.basename(l) for l in libs]
libs = [l[3:l.rfind('.')] for l in libs]
libs.append('CGAL')

# define the extension module
ext_mod =  Extension('_pymemsurfer',
                     sources = sources,
                     include_dirs=[headers, numpy.get_include(),
                                   BOOST_ipath, EIG_ipath, CGAL_ipath, VTK_ipath],
                     libraries=libs,
                     library_dirs=[VTK_lpath, CGAL_lpath],
                     language = 'c++',
                     swig_opts=['-c++', '-I'+headers],
                     define_macros=[('VTK_AVAILABLE',1), ('CGAL_AVAILABLE',1)],
                     extra_compile_args=['-std=c++11',
                                         '-Wno-inconsistent-missing-override',
                                         '-Wno-deprecated-declarations',
                                         '-Wno-unknown-pragmas',
                                         '-Wno-misleading-indentation'],
                     extra_link_args=['-std=c++11']
                    )


# ------------------------------------------------------------------------------
# install pypoisson as an Extension module
# included the setup functionality of https://github.com/mmolero/pypoisson
# ------------------------------------------------------------------------------

pp_path = os.path.join(os.getcwd(), 'pypoisson/src/PoissonRecon_v6_13/src/')
pp_files = [os.path.join(pp_path,x) for x in os.listdir(pp_path) if x.endswith('.cpp')]
pp_sources = ['pypoisson/src/pypoisson.pyx'] + pp_files
ext_pp = Extension('pypoisson', pp_sources,
                   language='c++',
                   include_dirs = [numpy.get_include()],
                   extra_compile_args = ['-w','-fopenmp'],
                   extra_link_args=['-fopenmp'],
                   define_macros = [('__linux__', '1')] # todo: need to add for mac
                  )
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# set up!
setup(name='memsurfer',
      version='0.1.0',
      description='Python tool to compute bilayer membranes.',
      author='Harsh Bhatia',
      author_email='hbhatia@llnl.gov',
      setup_requires=['cython'],
      packages=find_packages(),
      package_data={ 'memsurfer': ['_pymemsurfer.so', 'pypoisson.so'] },
      ext_modules=cythonize([ext_mod, ext_pp])
     )
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
