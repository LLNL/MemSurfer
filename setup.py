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
import os
import sys
import glob
import numpy
import socket
import getpass
import subprocess

from distutils import sysconfig
from pkg_resources import parse_version
from setuptools import find_packages, setup, Extension
from setuptools.command.build_ext import build_ext
from setuptools.command.build_py import build_py
from Cython.Build import cythonize

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
__version_info__ = ('1', '1', '0')
__version__ = '.'.join(__version_info__)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
def check_platform():
    if sys.platform.startswith('linux'):
        return 'linux'
    elif sys.platform == 'darwin':
        return 'darwin'
    elif sys.platform.startswith('win'):
        return 'win'
    raise Exception(f'Unknown platform ({sys.platform})!')


def shlib_extn():
    if sys.platform.startswith('linux'):
        return 'so'
    elif sys.platform == 'darwin':
        return 'dylib'
    elif sys.platform.startswith('win'):
        return 'dll'
    raise Exception(f'Unknown platform ({sys.platform})!')


def find_shlib_path(path, libname):
    extn = shlib_extn()
    options = ['lib', 'lib64']
    for o in options:
        if os.path.isfile(os.path.join(path, o, f'{libname}.{extn}')):
            return os.path.join(path, o)
    raise Exception(f'Find_shlib_path({path},{libname}) failed!')


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
def check_dependencies(packages):

    for k, v in packages.items():
        output = subprocess.run(v['cmd'],
                                capture_output=True,
                                universal_newlines=True,
                                shell=True,
                                check=True)

        if len(output.stderr) != 0:
            raise Exception(output.stderr)

        rval = [r for r in output.stdout.splitlines() if len(r) > 0]
        vstring = rval[0].split()[-1]

        if parse_version(vstring) < parse_version(v['min_version']):
            raise Exception(f"Insufficient version of ({k}): "
                            f"{vstring} < {v['min_version']}!")


def fetch_paths(default):
    # look for the paths of these packages
    packages = ['cgal', 'eigen', 'boost']

    # now, create a dictionary of all paths
    paths = {}
    for i in range(len(packages)):
        k = packages[i]
        kr = f'{k.upper()}_ROOT'

        p = os.path.expandvars(os.environ.get(kr, default))
        assert os.path.isdir(p), f'Did not find ({kr}={p})'

        print(f'  > {kr} = ({p})')
        paths[k] = {'root': p, 'include': os.path.join(p, 'include')}

    # now, update the include and lib paths
    paths['eigen'].update(include=os.path.join(paths['eigen']['include'], 'eigen3'))
    paths['cgal'].update({'lib': find_shlib_path(paths['cgal']['root'], 'libCGAL')})

    # make sure these paths exist!
    for pkg, pkgpath in paths.items():
        for k, v in pkgpath.items():
            if not os.path.exists(v):
                raise Exception(f'Did not find expected path ({v}) for ({pkg})')
    return paths


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
class CustomBuildExt(build_ext):

    def build_extensions(self):
        compiler_name = self.compiler.compiler[0]
        if not self.compiler._is_gcc(compiler_name):
            raise Exception(f'Need a GCC compiler. Found ({compiler_name})')
        super().build_extensions()

    def get_ext_filename(self, ext_name):
        filename = super().get_ext_filename(ext_name)
        suffix = sysconfig.get_config_var('EXT_SUFFIX')
        ext = os.path.splitext(filename)[1]
        return filename.replace(suffix, "") + ext


class CustomBuildPy(build_py):
    # https://stackoverflow.com/questions/12491328/python-distutils-not-include-the-swig-generated-module
    # https://stackoverflow.com/questions/29477298/setup-py-run-build-ext-before-anything-else/48942866#48942866
    def run(self):
        self.run_command("build_ext")
        return super().run()


# ------------------------------------------------------------------------------
# main function
# ------------------------------------------------------------------------------
if __name__ == '__main__':

    if sys.version_info[0] < 3 or sys.version_info[1] < 7:
        raise Exception('MemSurfer must be installed using Python 3.7 and above!')

    # figure out the path to the source directory of memsurfer
    PATH_MEM = os.path.dirname(os.path.abspath(__file__))
    print(f'> Installing MemSurfer for ({getpass.getuser()}) '
          f'on ({socket.gethostname()}) (platform={sys.platform})')
    pltform = check_platform()

    # --------------------------------------------------------------------------
    dependencies = {'swig': {'cmd': 'swig -version', 'min_version': '3.0.2'}}
    check_dependencies(dependencies)

    print(f'  > MemSurfer = ({PATH_MEM})')
    PATHS = fetch_paths(os.path.join(PATH_MEM, 'external'))

    # external libs
    extn = shlib_extn()
    LIBS_EXT = glob.glob(os.path.join(PATHS['cgal']['lib'], f'libCGAL*.{extn}'))
    LIBS_EXT = [os.path.basename(l) for l in LIBS_EXT]
    LIBS_EXT = [os.path.splitext(l)[0] for l in LIBS_EXT]
    LIBS_EXT = [l[3:] for l in LIBS_EXT]

    # --------------------------------------------------------------------------
    # install pypoisson as an Extension module
    # included the setup functionality of https://github.com/mmolero/pypoisson
    # --------------------------------------------------------------------------
    PATH_PP = os.path.join(PATH_MEM, 'pypoisson')
    PATH_PPR = os.path.join(PATH_PP,  'PoissonRecon_v6_13/src')

    # get all cpp files
    SRC_PP = [x for x in os.listdir(PATH_PPR) if x.endswith('.cpp')]

    # need to remove these files that expose main function
    exe_files = ['PoissonRecon.cpp', 'SurfaceTrimmer.cpp']
    SRC_PP = [x for x in SRC_PP if x not in exe_files]
    SRC_PP = [os.path.join(PATH_PPR, x) for x in SRC_PP]

    # add the pyx
    SRC_PP = [os.path.join(PATH_PP, 'pypoisson.pyx')] + SRC_PP

    # here, we assume that the compiler supports c++11
    # ideally, we would like to identify the compiler during configuration
    if pltform in ['linux', 'darwin']:
        macros = [('__linux__', 1)]
    else:
        macros = [('__linux__', 0)]

    EXT_PP = Extension('pypoisson', language='c++',
                       define_macros=macros,
                       extra_compile_args=['-w', '-fopenmp'],
                       extra_link_args=['-fopenmp'],
                       sources=SRC_PP,
                       include_dirs=[numpy.get_include()],
                       )

    # --------------------------------------------------------------------------
    # Now, build an extension module for MemSurfer's cpp code
    # --------------------------------------------------------------------------
    # cpp code
    PATH_PM = os.path.join(PATH_MEM, 'memsurfer')
    SRC_MEM = ['memsurfer.i', 'src/PointSet.cpp', 'src/TriMesh.cpp',
               'src/TriMesh_cgal.cpp', 'src/TriMesh_kde.cpp']
    SRC_MEM = [os.path.join(PATH_PM, f) for f in SRC_MEM]
    PATH_PM = os.path.join(PATH_PM, 'src')

    # define the extension module
    # this name should be _[name defined in memsurfer.i]
    # prepending with "memsurfer" puts this in the memsurfer namespace
    EXT_MEM = Extension('memsurfer._memsurfer_cmod', language='c++',
                        swig_opts=['-c++', '-I' + PATH_PM],
                        define_macros=[('CGAL_AVAILABLE', 1)],
                        extra_compile_args=['-std=c++11',
                                            '-Wno-parentheses',
                                            '-Wno-uninitialized',
                                            '-Wno-inconsistent-missing-override',
                                            '-Wno-deprecated-declarations',
                                            '-Wno-unknown-pragmas',
                                            '-Wno-misleading-indentation',
                                            '-Wno-unknown-warning-option'],
                        extra_link_args=['-std=c++11'],
                        sources=SRC_MEM,
                        include_dirs=[PATH_PM, numpy.get_include(),
                                      PATHS['boost']['include'],
                                      PATHS['eigen']['include'],
                                      PATHS['cgal']['include'],
                                      '/opt/local/include'
                                      ],
                        libraries=LIBS_EXT,
                        library_dirs=[PATHS['cgal']['lib']],
                        )

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # set up!
    setup(name='memsurfer',
          version=__version__,
          author='Harsh Bhatia',
          author_email='hbhatia@llnl.gov',
          url='https://github.com/LLNL/MemSurfer',
          description='Python tool to compute bilayer membranes.',
          license='GPL 3',

          install_requires=['numpy>=1.20', 'vtk>=9.1'],
          packages=find_packages(),
          ext_modules=cythonize([EXT_PP, EXT_MEM]),
          cmdclass={'build_py': CustomBuildPy, 'build_ext': CustomBuildExt}
          )

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
