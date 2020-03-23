## MemSurfer, Version 1.0
#### Released: Apr 17, 2019

##### Author: Harsh Bhatia (hbhatia@llnl.gov) and Peer-Timo Bremer

MemSurfer is an efficient and versatile tool to compute and analyze membrane surfaces found in a wide
variety of large-scale molecular simulations. MemSurfer works independent of the
type of simulation, directly on the 3D point coordinates, and can handle a variety of membranes as well as
atomic simulations. MemSurfer provides many in-built analysis
tasks, such as computing the membrane curvature, density and normals of lipids,
and area per lipid. More importantly, MemSurfer provides a simple-to-use
Python API that may be easily used/extended to perform other types of analysis.

### Dependencies

The interface to MemSurfer is provided via `python 3.7`, whereas the Core functionality is written in `C++ 11` with the following dependencies.
  - [The Computational Geometry Algorithms Library (CGAL)](https://www.cgal.org/ "CGAL"): `v 4.13.0`
  - [Visualization Toolkit (VTK)](https://www.vtk.org/ "VTK"): `v 8.1.2`
  - [boost](https://www.boost.org/): `v 1.66`
  - [Eigen](http://eigen.tuxfamily.org/index.php): `v 3.3.7`

Installing MemSurfer and its dependencies additionally require the following software.

  - `C++` compiler that supports `C++ 11` and [OpenMP](https://www.openmp.org/) (tested with `GNU gcc 7.3.0`)
  - `Python 3` interpreter (tested with `python 3.7.2`)
  - [Cython](https://cython.org/): `v 0.29.10`
  - [Swig](http://www.swig.org/): `v 3.0.12`
  - [CMake](https://cmake.org/): `v 3.13`


### Installation via Spack

MemSurfer is available via [spack](https://spack.io) -- a package manager for HPC. Please
download spack (see instructions provded by `spack`). Once installed,
please do the following
```
$ spack install memsurfer@develop +osmesa %gcc@7.3.0
```
The `+osmesa` variant enables OSMesa support for `VTK` . Here, the installation has been tested with `gcc@7.3.0`.

### Installation from source

MemSurfer can be downloaded using the following link.
```
$ git clone --recursive git@github.com:LLNL/MemSurfer.git
$ MEM_HOME=`pwd`/MemSurfer
```

***Note:*** If the cloning fails with the following error:
```
Cloning into '<your-path>/MemSurfer/pypoisson'...
git@github.com: Permission denied (publickey).
fatal: Could not read from remote repository.
```
This error means that your ssh keys are not registered with github. Please see [here](https://help.github.com/en/articles/connecting-to-github-with-ssh) and [here](https://help.github.com/en/articles/adding-a-new-ssh-key-to-your-github-account) to register your ssh key and retry.

#### 1. Installing MemSurfer's Dependencies

The following describes the installation of `cgal`, `vtk`, `boost` and `eigen`,
whereas the rest of the dependencies are assumed to be standard and available.


##### 1a. Eigen, Boost, CGAL, and VTK (on mac)
On `macosx`, the simplest way to install these is using [`macports`](macports.org). As of Jan 2019, `macports` installs the correct versions.
```
$ port install vtk+python
$ port install cgal           # also installs boost and eigen
```

##### 1b. Eigen, Boost, CGAL, and VTK (from source)


A helper script `$MEM_HOME/install_deps.sh` is provided to install these dependencies.
In order to use this script, you need to specify a `C++` compiler that supports `OpenMP`.

```
$ export CC=`which gcc`
$ export CXX=`which g++`
$ sh install_deps.sh
```

This script installs all these dependencies in a folder called `$MEM_HOME/external`.
Note that you may not need to install all four of these dependencies. Please
edit the script (lines 5--8) to select which ones to install.

Once installed successfully, please add the following to your shell profile to
access these dependencies.
```
$ export PYTHONPATH=$MEM_HOME/external/lib/python3.7/site-packages:$PYTHONPATH
$ export LD_LIBRARY_PATH=$MEM_HOME/external/lib:$MEM_HOME/external/lib64:$LD_LIBRARY_PATH
      # linux
$ export DYLD_LIBRARY_PATH=$MEM_HOME/external/lib:$MEM_HOME/external/lib64:$DYLD_LIBRARY_PATH
      # mac
```

#### 2. MemSurfer

Once the dependencies have been installed, `MemSurfer` can be installed
using `distutils`. However, you need to explicitly supply the path of the external
dependencies.
```
$ export BOOST_ROOT=<path_to_boost>
      # such that boost headers are contained in $BOOST_ROOT/include/boost
$ export VTK_ROOT=<path_to_vtk>
      # such that vtk headers are contained in $VTK_ROOT/include/vtk-8.1
$ export CGAL_ROOT=<path_to_cgal>
      # such that cgal headers are contained in $CGAL_ROOT/include/CGAL
$ export EIGEN_ROOT=<path_to_eigen>
      # such that eigen headers are contained in $EIGEN_ROOT/include/eigen3
```
Note that, all these paths default to `$MEM_HOME/external`. So you if you installed
the dependencies using the script provided above, you do not need to specify these
paths.

```
$ cd $MEM_HOME
$ CC=`which gcc` CXX=`which g++` LDCXXSHARED="`which g++` -bundle -undefined dynamic_lookup" \
  python setup.py install
```

### Examples

* See the `example` directory.

### Change Log

##### Mar 23, 2020

* Correctly normalize Gaussian KDE and support 3D Gaussian kernel.


### License

MemSurfer is released under GPU-3.0 license. See the `LICENSE` file for details.

*`LLNL-CODE-763493`*
