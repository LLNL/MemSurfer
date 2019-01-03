## MemSurfer, Version 0.1
#### Released: Jan 10, 2019

##### Author: Harsh Bhatia (hbhatia@llnl.gov) and Peer-Timo Bremer

MemSurfer is a tool to compute and analyze membrane surfaces found in a wide
variety of large-scale molecular simulations. MemSurfer works independent of the
type of simulation, directly on the 3D point coordinates. As a result, MemSurfer
can handle a variety of membranes, such as tethers and vesicles, as well as
atomic simulations. MemSurfer provides many in-built analysis
tasks, such as computing the membrane curvature, density and normals of lipids,
and area per lipid. More importantly, MemSurfer provides a simple-to-use
Python API that may be easily used/extended to perform other types of analysis.

### Installation

Core functionality of MemSurfer is written in `python` and `C++`, with the following
dependencies.
  - [The Computational Geometry Algorithms Library (CGAL)](https://www.cgal.org/ "CGAL"): `v 4.13`
    - [boost](https://www.boost.org/): `v 1.66` and above
    - [Eigen](http://eigen.tuxfamily.org/index.php): `v 3.3.7` and above
  - [Visualization Toolkit (VTK)](https://www.vtk.org/ "VTK"): `v 8.1.1`
    - [CMake](https://cmake.org/): `v 3.13` and above
  - [pypoisson](https://github.com/mmolero/pypoisson): Thanks to [Miguel Molero](https://github.com/mmolero/pypoisson)
    - [OpenMP](https://www.openmp.org/)
    - [Cython](https://cython.org/)
  - [Swig](http://www.swig.org/)

The following describes the installation of `cgal`, `vtk`, `eigen`, and `pypoisson`,
whereas the rest of the dependencies are assumed to be standard and available.

##### 1a. Eigen, CGAL, and VTK (on mac)
On `macosx`, the simplest way to install these is using [`macports`](macports.org). As of Jan 2019, `macports` installs the correct versions.
```
$ port install vtk+python
$ port install cgal           # also installs boost and eigen
```

##### 1b. Eigen, CGAL, and VTK (from source)
Here are the simplified instructions to install these dependencies from source.
More detailed instructions can be found from the respective websites.

*Note* that `PATH_DEP` refers to the path where you want to install the dependencies.

* `eigen3`
```
$ wget http://bitbucket.org/eigen/eigen/get/3.3.7.tar.bz2
$ tar -xf 3.3.7.tar.bz2
$ mkdir -p $PATH_DEP/include/eigen3
$ cp -r eigen-eigen-323c052e1731/Eigen $PATH_DEP/include/eigen3/
```

* `cgal`
```
$ wget https://github.com/CGAL/cgal/archive/releases/CGAL-4.13.tar.gz
$ tar -xf CGAL-4.13.tar.gz
$ mkdir -p cgal-releases-CGAL-4.13/build
$ cd cgal-releases-CGAL-4.13/build
$ cmake -DCMAKE_INSTALL_PREFIX:STRING=$PATH_DEP \
        -DCMAKE_BUILD_TYPE:STRING=Release \
        -DBUILD_SHARED_LIBS:BOOL=ON \
        -DWITH_CGAL_Qt5:BOOL=OFF\
        ..
$ make -j12
$ make install
```

* `vtk`
```
$ wget https://www.vtk.org/files/release/8.1/VTK-8.1.2.tar.gz
$ tar -xf VTK-8.1.2.tar.gz
$ mkdir -p VTK-8.1.2/build
$ cd VTK-8.1.2/build
$ cmake -DCMAKE_INSTALL_PREFIX:STRING=$PATH_DEP \
        -DCMAKE_BUILD_TYPE:STRING=Release \
        -DCMAKE_CXX_FLAGS:STRING="-Wno-inconsistent-missing-override"  \
        -DBUILD_SHARED_LIBS:BOOL=ON \
        -DVTK_WRAP_PYTHON:BOOL=ON \
        ..
$ make -j12
$ make install
```



##### 2. pypoisson
The `pypoisson` package requires a `C++` compiler that supports `openmp`.
*Note* that on `mac`, you may not have `openmp` installed with `clang` compiler. If you
wish to instead use `gnu` compiler, you can force `distutils` to use a compiler
different than the default for your `python` installation.

```
$ git clone --recursive git://github.com/mmolero/pypoisson.git
$ cd pypoisson

  # on linux, or on mac if your python was installed with a compiler supporting openmp
$ python setup.py build

  # build command on mac, when you want to use your default gnu compiler that supports openmp
$ CC=`which gcc` CXX=`which g++` LDCXXSHARED="`which g++` -bundle -undefined dynamic_lookup" python setup.py build

  # installation
$ python setup.py install --prefix=$PATH_DEP
```

##### 3. Testing the dependencies

Please update your python path to link to the dependencies. *Note* that it is
best to have this added to your `.bashrc` file.
```
$ export PYTHONPATH=$PATH_DEP/lib/python2.7/site-packages:$PYTHONPATH
```
Now, check if you are able to load these modules.
```
$ python
Python 2.7.15 (default, Sep 12 2018, 13:32:25)
[GCC 4.2.1 Compatible Apple LLVM 9.1.0 (clang-902.0.39.2)] on darwin
Type "help", "copyright", "credits" or "license" for more information.

>>> import vtk
>>> print vtk.__file__
<PATH_DEP>/lib/python2.7/site-packages/vtk/__init__.pyc

>>> import pypoisson
>>> print pypoisson.__file__
<PATH_DEP>/lib/python2.7/site-packages/pypoisson.so
```

##### 4. MemSurfer

Once the dependencies have been installed, `MemSurfer` can be installed simply
using `distutils`. However, you need to explicitly supply the path of the external
dependencies.
```
$ export BOOST_ROOT=<path_to_boost>
    # such that boost headers are contained in BOOST_ROOT/include/boost
$ export VTK_ROOT=$PATH_DEP
$ export CGAL_ROOT=$PATH_DEP

$ git clone https://github.com/LLNL/MemSurfer.git
$ cd MemSurfer
$ python setup.py install
```

### Examples

* See the `example` directory. The examples therein *(to be added soon)* require `MDAnalysis`, which can be installed as `pip install mdanalysis`.

### License

MemSurfer is released under GPU-3.0 license. See the `LICENSE` file for details.

*`LLNL-CODE-763493`*
