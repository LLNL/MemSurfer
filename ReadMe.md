## MemSurfer, Version 0.1
#### Released: Jan 10, 2019

##### Author: Harsh Bhatia (hbhatia@llnl.gov) and Peer-Timo Bremer

MemSurfer is an efficient and versatile tool to compute and analyze membrane surfaces found in a wide
variety of large-scale molecular simulations. MemSurfer works independent of the
type of simulation, directly on the 3D point coordinates, and can handle a variety of membranes as well as
atomic simulations. MemSurfer provides many in-built analysis
tasks, such as computing the membrane curvature, density and normals of lipids,
and area per lipid. More importantly, MemSurfer provides a simple-to-use
Python API that may be easily used/extended to perform other types of analysis.

### Dependencies

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


### Installation

MemSurfer can be installed using [`spack`](https://spack.io). `spack` is a package manager designed to simplify installation process by providing consistent builds of all dependencies, as well as supporting different configurations simultaneously.
In the following, we use `spack` to install the complete dependency tree of MemSurfer (even including `python`) to simplify the process for the user without figuring out the correct versions, paths, etc. However, this conservative approach can take a long time (up to half an hour) for the complete installation.
Therefore, while we recommend using `spack` for easy installation of MemSurfer, an expert user may also choose to perform manual installation of only the required components.

### Installation using `spack`
```
SPACK_ROOT=`pwd`
```

#### 1. Download `spack`
```
git clone https://github.com/spack/spack.git
```

#### 2. Activate `spack`
```
source $SPACK_ROOT/spack/share/spack/setup-env.sh
```
**Note:** You will need to activate `spack` every time you start a new session. You may consider putting this line in your profile file (e.g., `.bashrc`), in which case, your terminal will do it for you upon start.

#### 3. Configure `spack`.

  **3a) Compiler**. You can ask `spack` to search through all the available compilers on your machine, and then look at the available compilers.
  ```
  spack compiler add
  spack compiler list
  ```
  MemSurfer requires `gcc 7.1.0` (and above). You may already have the required compiler, in which case, you can move to the next step. If not, you can ask `spack` to install the compiler for you.
  ```
  spack install gcc@7.1.0
  ```

  **3b) MPI**. `spack` needs to know which MPI will you want to use.  Typically, MPI implementations are ported with the OS, so you will already have them installed. In this step, we will simply point `spack` to your MPI installation by editing either of the following two files.
  ```
  # Users can override these settings by editing the following files.
  #
  # Per-spack-instance settings (overrides defaults):
  #   $SPACK_ROOT/etc/spack/packages.yaml
  #
  # Per-user settings (overrides default and site settings):
  #   ~/.spack/packages.yaml
  ```
  In `package.yaml`, you need to specify two details: providers for MPI, and the available versions and corresponding paths. For example,
  ```
  packages:
    all:
      providers:
        mpi: [openmpi, mvapich2]
    mvapich2:
      buildable: False
      version: [2.2]
      paths:
        mvapich2@2.2 %gcc@4.9.3 arch=linux-rhel7-x86_64: <path_to_mvapich2@2.2_gcc@4.9.3>
        mvapich2@2.2 %gcc@7.1.0 arch=linux-rhel7-x86_64: <path_to_mvapich2@2.2_gcc@7.1.0>
  ```
  This configuration specifies that when a package is to be compiled with `mvapich2`, `spack` should not try to build MPI, but instead, use the system-installed version at the given path. Most machines install MPI in system locations, such as `\opt` and `\usr`. You also need to specify the architecture of your machine, which you can find using
  ```
  spack arch
  ```
  Finally, for MemSurfer, we want to use `gcc @7.1.0`.
  You can find more details on `spack` configuration [here](https://spack.readthedocs.io/en/latest/build_settings.html#build-settings).

  **3c) OpenGL**. Similar to above, we need to specify which OpenGL implementation to use for installation. `vtk`, which is a dependency of MemSurfer, supports two implementations: `opengl` and `mesa`. You may choose either. After this update, your `package.yaml` would look like this.
  ```
  packages:
    all:
      providers:
        mpi: [openmpi, mvapich2]
        gl: [opengl, mesa]
    mvapich2:
      buildable: False
      version: [2.2]
      paths:
        mvapich2@2.2 %gcc@4.9.3 arch=linux-rhel7-x86_64: <path_to_mvapich2@2.2_gcc@4.9.3>
        mvapich2@2.2 %gcc@7.1.0 arch=linux-rhel7-x86_64: <path_to_mvapich2@2.2_gcc@7.1.0>
    opengl:
      buildable: False
      paths:
        opengl@4.5.0 arch=linux-rhel7-x86_64: <path_to_opengl>
    mesa:
      buildable: False
      paths:
        mesa@17.3 arch=linux-rhel7-x86_64:  <path_to_mesa>
  ```
  ***Note:*** A current limitation of `spack` is that when looking for opengl libraries, it only looks in `path/lib`, and not in `path/lib64`. To get around this limitation, you can use symbolic links to route `spack` to the actual locations of the files. See [here](https://groups.google.com/forum/#!msg/spack/270j1Ooi9VY/mvMpFZ8-BAAJ).

#### 4. Install MemSurfer
```
spack install memsurfer@1.0 +vtkmesa ^mvapich2 %gcc@7.1.0     (if using mesa, see 3c above)
spack install memsurfer@1.0 ~vtkmesa ^mvapich2 %gcc@7.1.0     (otherwise)
```

#### 5. Load MemSurfer
```
source <your_path>/spack/share/spack/setup-env.sh
spack load --dependencies memsurfer@1.0
```

### Manual Installation

The following describes the installation of `cgal`, `vtk`, and `eigen`,
whereas the rest of the dependencies are assumed to be standard and available.

#### 1a. Eigen, CGAL, and VTK (on mac)
On `macosx`, the simplest way to install these is using [`macports`](macports.org). As of Jan 2019, `macports` installs the correct versions.
```
$ port install vtk+python
$ port install cgal           # also installs boost and eigen
```

#### 1b. Eigen, CGAL, and VTK (from source)
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

#### 2. MemSurfer

Once the dependencies have been installed, `MemSurfer` can be installed simply
using `distutils`. However, you need to explicitly supply the path of the external
dependencies.
```
$ export BOOST_ROOT=<path_to_boost>
    # such that boost headers are contained in BOOST_ROOT/include/boost
$ export VTK_ROOT=$PATH_DEP
$ export CGAL_ROOT=$PATH_DEP
$ export EIGEN_ROOT=$PATH_DEP

$ git clone --recursive https://github.com/LLNL/MemSurfer.git
$ cd MemSurfer
$ python setup.py install
```

Please update your python path to link to the dependencies. *Note* that it is
best to have this added to your `.bashrc` file.
```
$ export PYTHONPATH=$PATH_DEP/lib/python2.7/site-packages:$PYTHONPATH
```

### Examples

* See the `example` directory.

### License

MemSurfer is released under GPU-3.0 license. See the `LICENSE` file for details.

*`LLNL-CODE-763493`*
