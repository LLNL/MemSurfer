## MemSurfer, Version 1.0.1
#### Released: July 22, 2021

##### Author: Harsh Bhatia (hbhatia@llnl.gov) and Peer-Timo Bremer

MemSurfer is an efficient and versatile tool to compute and analyze membrane surfaces found in a wide
variety of large-scale molecular simulations. MemSurfer works independent of the
type of simulation, directly on the 3D point coordinates, and can handle a variety of membranes as well as
atomic simulations. MemSurfer provides many in-built analysis
tasks, such as computing the membrane curvature, density and normals of lipids,
and area per lipid. More importantly, MemSurfer provides a simple-to-use
Python API that may be easily used/extended to perform other types of analysis.

### Dependencies

The interface to MemSurfer is provided via `python 3.7`, whereas the core functionality is written in `C++ 11` with the following dependencies.
  1. [The Computational Geometry Algorithms Library (CGAL)](https://www.cgal.org/ "CGAL"): `v 4.13.0`
  2. [Visualization Toolkit (VTK)](https://www.vtk.org/ "VTK"): `v 8.1.2`
  3. [boost](https://www.boost.org/): `v 1.66`
  4. [Eigen](http://eigen.tuxfamily.org/index.php): `v 3.3.9`

Installing MemSurfer and its dependencies additionally require the following software.

  - `C++` compiler that supports `C++ 11` and [OpenMP](https://www.openmp.org/) (tested with `GNU gcc 7.3.0` and `GNU gcc 7.5.0`)
    - **Please use `gcc@7`. See known issues below.**
  - `Python 3` interpreter (tested with `3.7.2` and `3.7.11`)
  - [Cython](https://cython.org/) (tested with `0.29.10` and  `0.29.24`)
  - [Swig](http://www.swig.org/) (tested with `3.0.12` and `4.0.2`)
  - [CMake](https://cmake.org/) (tested with `3.13` and `3.20`)
  - [gmp](https://gmplib.org/) (tested with `6.2.1`)
  - [mpfr](https://www.mpfr.org/) (tested with `4.1.0`)

All the additional packages (unnumbered) are usually available on all machines
through standard package managers, for example, `macports`, `apt-get`, `yum`, etc.
Please make use of these package managers to avail these dependencies.

The numbered dependencies may not be available, or package managers may not
give full control on the version number to be installed. Here, we will discuss
the installation of these major dependencies.


### Installation via Spack

MemSurfer is available via [spack](https://spack.io) -- a package manager for HPC. Please
download spack (see instructions provded by `spack`). Once installed,
please do the following
```
$ spack install memsurfer@1.0.1 ^python@3.7.3 %gcc@7.3.0
```


### Installation from source

MemSurfer can be downloaded using the following link.
```
$ git clone --recursive git@github.com:LLNL/MemSurfer.git
$ MEM_HOME=`pwd`/MemSurfer
```

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
$ export CC_COMPILER=`which gcc`
$ export CXX_COMPILER=`which g++`
$ sh install_deps.sh
```

This script installs all these dependencies in a folder called `$MEM_HOME/external`.
Note that you may not need to install all four of these dependencies. Please
edit the script (lines 6--9) to select which ones to install.

Once installed successfully, please add the following to your shell profile to
access these dependencies.
```
$ export PYTHONPATH=$MEM_HOME/external/lib/python3.7/site-packages:$PYTHONPATH

# LD_LIBRARY_PATH for linux
$ export LD_LIBRARY_PATH=$MEM_HOME/external/lib:$MEM_HOME/external/lib64:$LD_LIBRARY_PATH

# DYLD_LIBRARY_PATH for mac
$ export DYLD_LIBRARY_PATH=$MEM_HOME/external/lib:$MEM_HOME/external/lib64:$DYLD_LIBRARY_PATH
```

#### 2. MemSurfer

Once the dependencies have been installed, `MemSurfer` can be installed
using `distutils`. However, you need to explicitly supply the path of the external
dependencies.
```
$ export BOOST_ROOT=<path_to_boost>   # boost headers are contained in $BOOST_ROOT/include/boost
$ export VTK_ROOT=<path_to_vtk>       # vtk headers are contained in $VTK_ROOT/include/vtk-8.1
$ export CGAL_ROOT=<path_to_cgal>     # cgal headers are contained in $CGAL_ROOT/include/CGAL
$ export EIGEN_ROOT=<path_to_eigen>   # eigen headers are contained in $EIGEN_ROOT/include/eigen3
```
Note that, all these paths default to `$MEM_HOME/external`. So if you installed
the dependencies using the script provided above, you do not need to specify these
paths.

```
$ cd $MEM_HOME
$ CC=`which gcc` CXX=`which g++` LDCXXSHARED="`which g++` -bundle -undefined dynamic_lookup" \
  python setup.py install
```

### Trobubleshooting and Known Issues

***Jul 20, 2021:*** It appears that compilation using `gcc@8` results in a `segfault` or
a `malloc error` (see issue [#10](https://github.com/LLNL/MemSurfer/issues/10)).
The problem appears to be residing in a dependency called `pypoisson`. Until
it is resolved, please use `gcc@7`.

***Jan 23, 2019:*** If the cloning fails with the following error:
```
Cloning into '<your-path>/MemSurfer/pypoisson'...
git@github.com: Permission denied (publickey).
fatal: Could not read from remote repository.
```
This error means that your ssh keys are not registered with github. Please see [here](https://help.github.com/en/articles/connecting-to-github-with-ssh) and [here](https://help.github.com/en/articles/adding-a-new-ssh-key-to-your-github-account) to register your ssh key and retry.


### Examples

* See the `example` directory.

### Change Log

##### Jul 22, 2021

* Upgraded to Eigen 3.3.9
* Updated `pypoisson` submodule to [commit #94534a2](https://github.com/mmolero/pypoisson/commit/94534a28e063b2d0ab8b8239e4ad0034c3613ec8).
* Added a specific guideline to use gcc@7.

##### Mar 23, 2020

* Correctly normalize Gaussian KDE and support 3D Gaussian kernel.
* Ported to Python3.


### License

MemSurfer is released under GNU GPL-3.0 license. See the `LICENSE` file for details.

*`LLNL-CODE-763493`*
