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
import getpass
import glob
import numpy
import os
import socket
import subprocess
import sys

from pkg_resources import parse_version
from distutils import sysconfig
from setuptools import find_packages, setup, Extension
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
__version_info__ = ('1', '0', '0')
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


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
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


def fetch_paths(default):
    # look for the paths of these packages
    packages = ['vtk', 'cgal', 'eigen', 'boost']

    # now, create a dictionary of all paths
    paths = {}
    for i in range(len(packages)):
        k = packages[i]
        try:
            p = os.path.expandvars(os.environ[f'{k.upper()}_ROOT'])
        except KeyError:
            p = default
        print(f'  > {k.upper()} = ({p})')
        paths[k] = {'root': p, 'include': os.path.join(p, 'include')}

    # now, update the include paths
    paths['vtk'].update(include=os.path.join(paths['vtk']['include'], 'vtk-8.1'))
    paths['eigen'].update(include=os.path.join(paths['eigen']['include'], 'eigen3'))

    # now, add the lib paths where needed
    paths['vtk'].update({'lib': os.path.join(paths['vtk']['root'], 'lib')})
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


class CustomInstall(install):

    def run(self):
        self.run_command('build_ext')
        self.do_egg_install()

        # ----------------------------------------------------------------------
        # currently, the _pymemsurfer.so file is getting installed
        # in the egg directory but not in the memsurfer package
        # until I figure out how to fix that, just create a symlink
        if True:
            libname = '_pymemsurfer.so'

            dist_name = self.config_vars['dist_name']
            dist_fullname = self.config_vars['dist_fullname']
            py_version = self.config_vars['py_version_short']
            platform = os.path.basename(self.build_lib)
            platform = platform.replace('lib.', '').replace(f'-{py_version}', '')

            egg_name = f'{dist_fullname}-py{py_version}-{platform}.egg'
            egg_path = os.path.join(self.install_libbase, egg_name)

            pwd = os.getcwd()
            os.chdir(os.path.join(egg_path, dist_name))

            src = os.path.join('..', libname)
            trg = os.path.join('.', libname)
            print(f'Symlinking ({src}) to ({trg}) within ({os.getcwd()})')
            os.symlink(src, trg)
            os.chdir(pwd)
        # ----------------------------------------------------------------------


# ------------------------------------------------------------------------------
# main function
# ------------------------------------------------------------------------------
if __name__ == '__main__':

    if sys.version_info[0] < 3 or sys.version_info[1] < 7:
        raise Exception('MemSurfer must be installed using Python 3.7 and above!')

    # figure out the path to the source directory of memsurfer
    PATH_MEM = os.path.dirname(os.path.abspath(__file__))
    print(f'> Installing MemSurfer for ({getpass.getuser}) '
          f'on ({socket.gethostname}) (platform={sys.platform})')
    pltform = check_platform()

    # --------------------------------------------------------------------------
    dependencies = {'swig': {'cmd': 'swig -version', 'min_version': '3.0.2'}}
    check_dependencies(dependencies)

    print(f'  > MemSurfer = ({PATH_MEM})')
    PATHS = fetch_paths(os.path.join(PATH_MEM, 'external'))

    # --------------------------------------------------------------------------
    # install pypoisson as an Extension module
    # included the setup functionality of https://github.com/mmolero/pypoisson
    # --------------------------------------------------------------------------
    PATH_PP = os.path.join(PATH_MEM, 'pypoisson/src/PoissonRecon_v6_13/src/')
    SRC_PP = [os.path.join(PATH_PP, x) for x in os.listdir(PATH_PP) if x.endswith('.cpp')]
    SRC_PP = ['pypoisson/src/pypoisson.pyx'] + SRC_PP

    # here, we assume that the compiler supports c++11
    # ideally, we would like to identify the compiler during configuration
    if pltform in ['linux', 'darwin']:
        macros = [('__linux__', '1')]
    else:
        macros = [('__linux__', '0')]

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
    INC_MEM = os.path.join(PATH_MEM, 'memsurfer', 'include')
    SRC_MEM = glob.glob(os.path.join(PATH_MEM, 'memsurfer', 'src', '*.cpp'))
    SRC_MEM.append(os.path.join(PATH_MEM, 'memsurfer', 'pymemsurfer.i'))

    # external libs
    LIBS_VTK = glob.glob(os.path.join(PATHS['vtk']['lib'], 'libvtk*'))
    LIBS_CGAL = glob.glob(os.path.join(PATHS['cgal']['lib'], 'libCGAL*'))
    LIBS_EXT = LIBS_VTK + LIBS_CGAL

    LIBS_EXT = list(set([os.path.realpath(l) for l in LIBS_EXT]))
    LIBS_EXT = [os.path.basename(l) for l in LIBS_EXT]
    LIBS_EXT = [os.path.splitext(l)[0] for l in LIBS_EXT]
    LIBS_EXT = [l[3:] for l in LIBS_EXT]

    # define the extension module
    EXT_MEM = Extension('_pymemsurfer', language='c++',
                        swig_opts=['-c++', '-I' + INC_MEM],
                        define_macros=[('VTK_AVAILABLE', 1), ('CGAL_AVAILABLE', 1)],
                        extra_compile_args=['-std=c++11',
                                            '-Wno-parentheses',
                                            '-Wno-inconsistent-missing-override',
                                            '-Wno-deprecated-declarations',
                                            '-Wno-unknown-pragmas',
                                            '-Wno-misleading-indentation',
                                            '-Wno-unknown-warning-option'],
                        extra_link_args=['-std=c++11'],
                        sources=SRC_MEM,
                        include_dirs=[INC_MEM, numpy.get_include(),
                                      PATHS['boost']['include'],
                                      PATHS['eigen']['include'],
                                      PATHS['cgal']['include'],
                                      PATHS['vtk']['include'],
                                      '/opt/local/include'
                                      ],
                        libraries=LIBS_EXT,
                        library_dirs=[PATHS['vtk']['lib'], PATHS['cgal']['lib']],
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

          packages=find_packages(),
          package_data={'memsurfer': ['_pymemsurfer.so', 'pypoisson.so']},
          ext_modules=[EXT_PP,EXT_MEM],
          cmdclass={'build_ext': CustomBuildExt, 'install': CustomInstall}
          )
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
