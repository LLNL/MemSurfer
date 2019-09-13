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
import os, sys, glob
import getpass, socket, subprocess
import numpy

from pkg_resources import parse_version
from Cython.Build import cythonize
from setuptools import find_packages, setup, Extension
from setuptools.command.install import install
from distutils.command.build import build

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
def is_newer_than_python37():
    return sys.version_info[0] == 3 and sys.version_info[1] >= 7

def check_platform():
    if sys.platform.startswith('linux'):    return
    elif sys.platform == 'darwin':          return
    elif sys.platform.startswith('win'):    return
    raise Exception('Unknown platform ({})!'.format(sys.platform))

def check_dependencies(packages):

    def run_cmd(cmd):
      output = subprocess.run(cmd, capture_output=True, universal_newlines=True,
                                   shell=True, check=True)
      if len(output.stderr) != 0:
        raise Exception (output.stderr)
      return output

    for k, v in packages.items():

        output = run_cmd(v['cmd'])
        rval = [r for r in output.stdout.splitlines() if len(r) > 0]
        vstring = rval[0].split()[-1]
        if parse_version(vstring) < parse_version(v['min_version']):
            msg = 'Incorrect version: ({} {}); need at least ({})!'
            raise Exception (msg.format(k, vstring, parse_version(v['min_version'])))

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
def shlib_extn():
    if sys.platform.startswith('linux'):    return 'so'
    elif sys.platform == 'darwin':          return 'dylib'
    elif sys.platform.startswith('win'):    return 'dll'
    raise Exception('Unknown platform ({})!'.format(sys.platform))

def find_shlib_path(path, libname):
    extn = shlib_extn()
    options = ['lib', 'lib64']
    for o in options:
      if os.path.isfile(os.path.join(path, o, '{}.{}'.format(libname, extn))):
        return os.path.join(path, o)
    raise Exception('Find_shlib_path({},{}) failed!'.format(path, libname))

def fetch_paths(default):

    # look for the paths of these packages
    packages = ['vtk', 'cgal', 'eigen', 'boost']

    # now, create a dictionary of all paths
    paths = {}
    for i in range(4):
        k = packages[i]
        try:
          p = os.path.expandvars(os.environ[k.upper()+'_ROOT'])
        except KeyError:
          p = default
        print ('  > {} = ({})'.format(k.upper(), p))
        paths[k] = {'root': p, 'include': os.path.join(p, 'include')}

    # now, update the include paths
    paths['vtk'].update(include = os.path.join(paths['vtk']['include'], 'vtk-8.1'))
    paths['eigen'].update(include = os.path.join(paths['eigen']['include'], 'eigen3'))

    # now, add the lib paths where needed
    paths['vtk'].update({'lib': os.path.join(paths['vtk']['root'], 'lib')})
    paths['cgal'].update({'lib': find_shlib_path(paths['cgal']['root'], 'libCGAL')})

    # TODO: make sure these paths exist!
    return paths

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
class CustomBuild(build):
    def run(self):
        self.run_command('build_ext')
        build.run(self)


class CustomInstall(install):
    def run(self):
        self.run_command('build_ext')
        self.do_egg_install()

# ------------------------------------------------------------------------------
# main function
# ------------------------------------------------------------------------------
if __name__ == '__main__' :

    if not is_newer_than_python37():
        raise Exception('MemSurfer must be installed using Python 3.7 and above!')

    # figure out the path to the source directory of memsurfer
    PATH_SELF = os.path.abspath(sys.argv[0])
    PATH_MEM = os.path.split(PATH_SELF)[0]

    print ('> Installing MemSurfer for ({}) on ({}) (platform={})'.format(getpass.getuser(), socket.gethostname(), sys.platform))

    check_platform()

    dependencies = {'swig': {'cmd': 'swig -version', 'min_version': '3.0.2'}}
    check_dependencies(dependencies)

    print ('  > MemSurfer = ({})'.format(PATH_MEM))
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
    EXT_PP = Extension('pypoisson', SRC_PP,
                       language='c++',
                       include_dirs = [numpy.get_include()],
                       extra_compile_args = ['-w','-fopenmp'],
                       extra_link_args=['-fopenmp'],
                       define_macros = [('__linux__', '1')] # todo: need to add for mac
                      )

    # --------------------------------------------------------------------------
    # Now, build an extension module for MemSurfer's cpp code
    # --------------------------------------------------------------------------
    # cpp code
    INC_MEM = os.path.join(PATH_MEM, 'memsurfer', 'include')
    SRC_MEM = glob.glob(os.path.join(PATH_MEM, 'memsurfer', 'src', '*.cpp'))
    SRC_MEM.append(os.path.join(PATH_MEM, 'memsurfer', 'pymemsurfer.i'))

    # external libs
    LIBS_EXT = glob.glob(os.path.join(PATHS['vtk']['lib'], 'libvtk*'))
    LIBS_EXT = [l for l in LIBS_EXT if os.path.islink(l)]
    LIBS_EXT = [os.path.basename(l) for l in LIBS_EXT]
    LIBS_EXT = [l[3:l.rfind('.')] for l in LIBS_EXT]
    LIBS_EXT.append('CGAL')

    # define the extension module
    EXT_MEM = Extension('_pymemsurfer',
                         sources = SRC_MEM,
                         include_dirs=[INC_MEM, numpy.get_include(),
                                       PATHS['boost']['include'],
                                       PATHS['eigen']['include'],
                                       PATHS['cgal']['include'],
                                       PATHS['vtk']['include']],
                         libraries=LIBS_EXT,
                         library_dirs=[PATHS['vtk']['lib'], PATHS['cgal']['lib']],
                         language = 'c++',
                         swig_opts=['-c++', '-I'+INC_MEM],
                         define_macros=[('VTK_AVAILABLE',1), ('CGAL_AVAILABLE',1)],
                         extra_compile_args=['-std=c++11',
                                             '-Wno-inconsistent-missing-override',
                                             '-Wno-deprecated-declarations',
                                             '-Wno-unknown-pragmas',
                                             '-Wno-misleading-indentation',
                                             '-Wno-unknown-warning-option'],
                         extra_link_args=['-std=c++11']
                        )

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # set up!
    setup(name='memsurfer',
          version='1.0.0',
          author='Harsh Bhatia',
          author_email='hbhatia@llnl.gov',
          url='https://github.com/LLNL/MemSurfer',
          description='Python tool to compute bilayer membranes.',
          license='GPL 3',

          setup_requires=['cython'],
          packages=find_packages(),
          package_data={ 'memsurfer': ['_pymemsurfer.so', 'pypoisson.so'] },
          ext_modules=cythonize([EXT_PP, EXT_MEM]),
          cmdclass={'build': CustomBuild, 'install': CustomInstall}
         )
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
