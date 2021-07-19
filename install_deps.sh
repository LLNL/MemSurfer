#!/bin/bash

# ------------------------------------------------------------------------------
# parameters to control this script!

INSTALL_VTK=true
INSTALL_CGAL=true
INSTALL_BOOST=true
INSTALL_EIGEN=true

PATH_MemSurfer=`pwd`
PATH_Ext=$PATH_MemSurfer/external

NPROCS=20

# ------------------------------------------------------------------------------
# install script for MemSurfer and its depenencies
# ------------------------------------------------------------------------------
echo ''
echo "> Installing MemSurfer's dependencies for (`whoami`) on (`hostname`). platform = (`uname`)"
echo " > MemSurfer Root Directory: $PATH_MemSurfer"

# ------------------------------------------------------------------------------
# check commands

command -v wget >/dev/null 2>&1 || { echo >&2 "Cannot find command 'wget'. Aborting."; exit 1;  }
command -v tar >/dev/null 2>&1 || { echo >&2 "Cannot find command 'tar'. Aborting."; exit 1;    }
command -v cmake >/dev/null 2>&1 || { echo >&2 "Cannot find command 'cmake'. Aborting."; exit 1;  }
command -v swig >/dev/null 2>&1 || { echo >&2 "Cannot find command 'swig'. Aborting."; exit 1;  }
command -v python3 >/dev/null 2>&1 || { echo >&2 "Cannot find command 'python3'. Aborting."; exit 1;  }
command -v cython >/dev/null 2>&1 || { echo >&2 "Cannot find command 'cython'. Aborting."; exit 1;  }

if [ -n "${CXX}" ] ; then
    command -v "${CXX}" >/dev/null 2>&1 || { echo >&2 "Cannot find command '${CXX}'. Aborting."; exit 1;  }
    echo " > Using `${CXX}  --version | head -n1`":  `type ${CXX}`
else
    echo " > Using `gcc --version | head -n1`":  `type gcc`
fi

# ------------------------------------------------------------------------------
# common utilities
# ------------------------------------------------------------------------------
_fetch() {
    httplink=$1
    filename=$2
    logname=$3
    dirname=${filename%.*}

    rm $logname 2>/dev/null
    if ! ls $filename 1> /dev/null 2>&1 ; then
        echo "     > Downloading ($httplink)"
        wget -o $logname $httplink
    fi
    if ! ls $dirname 1> /dev/null 2>&1 ; then
        echo "     > Extracting ($filename) to ($dirname)"
        mkdir $dirname
        tar -xf $filename -C $dirname --strip-components=1
    fi
}
_build() {
    packagename=$1
    logname=$2

    rm $logname 2>/dev/null
    echo "     > Building ($packagename)"
    make -j$NPROCS > $logname
    echo "     > Installing ($packagename)"
    make install -j$NPROCS > $logname
}
_test() {
    testfile=$1
    if ! ls $testfile 1> /dev/null 2>&1 ; then
        found=0
    else
        found=1
    fi
}
_test2() {
   _test $1
   if [[ "$found" == 1 ]]; then
     found=1
   else
     _test $2
     if [[ "$found" == 1 ]]; then
        found=1
     else
        found=0
     fi
   fi
}

# ------------------------------------------------------------------------------
mkdir -p $PATH_Ext/downloads

# ------------------------------------------------------------------------------
# Eigen 3.3.7
# ------------------------------------------------------------------------------
if [ "$INSTALL_EIGEN" = true ] ; then

    VERSION='3.3.7'
    NAME='Eigen-'$VERSION
    FILE=$VERSION.tar.bz2
    URL=https://gitlab.com/libeigen/eigen/-/archive/$VERSION/$FILE
    TEST=$PATH_Ext/include/eigen3

    _test $TEST
    if [[ "$found" == 1 ]]; then
        echo "   > ($NAME) is already installed at the requested location ($PATH_Ext)."

    else
        cd $PATH_Ext/downloads
        echo "   > Installing ($NAME)"

        # fetch the code base
        _fetch $URL $FILE $PATH_Ext/$NAME.download.log

        # configure step is custom for each package
        BUILD_PATH=$dirname/build
        mkdir -p $BUILD_PATH; cd $BUILD_PATH
        echo "     > Configuring out-of-source build (`pwd`)"

        rm $PATH_Ext/$NAME.cmake.log 2>/dev/null
        cmake -DCMAKE_INSTALL_PREFIX:STRING=$PATH_Ext \
              .. > $PATH_Ext/$NAME.cmake.log

        # eigen needs only make install
        rm $PATH_Ext/$NAME.make.log 2>/dev/null
        echo "     > Installing ($NAME)"
        make install -j$NPROCS > $NAME.make.log
    fi
fi

# ------------------------------------------------------------------------------
# Boost 1.66
# ------------------------------------------------------------------------------
if [ "$INSTALL_BOOST" = true ] ; then

    VERSION='1.66.0'
    FILE='boost_1_66_0.tar.gz'
    NAME='boost-'$VERSION
    URL=https://dl.bintray.com/boostorg/release/$VERSION/source/$FILE
    TEST=$PATH_Ext/lib/libboost_graph*

    _test $TEST
    if [[ "$found" == 1 ]]; then
        echo "   > ($NAME) is already installed at the requested location ($PATH_Ext)."

    else
        cd $PATH_Ext/downloads
        echo "   > Installing ($NAME)"

        # fetch the code base
        _fetch $URL $FILE $PATH_Ext/$NAME.download.log

        # configure step is custom for each package
        BUILD_PATH=$dirname
        cd $BUILD_PATH
        echo "     > Configuring $NAME (`pwd`)"

        rm $PATH_Ext/$NAME.cmake.log 2>/dev/null

        ./bootstrap.sh --prefix=$PATH_Ext \
                       --with-python=`which python3` \
                       --with-libraries=graph \
                       > $PATH_Ext/$NAME.bootstrap.log

        echo "     > Building and Installing $NAME"
        ./b2 install > $PATH_Ext/$NAME.make.log

        # test the installation
        _test $TEST
        if [[ "$found" == 1 ]]; then
            echo "   > ($NAME) successfully installed at ($PATH_Ext)."
        else
            echo "   > ($NAME) installation failed. Please see build logs in ($PATH_Ext) for more information."
        fi
    fi
fi

# ------------------------------------------------------------------------------
# CGAL 4.13
# ------------------------------------------------------------------------------
if [ "$INSTALL_CGAL" = true ] ; then

    VERSION='4.13'
    NAME='CGAL-'$VERSION
    FILE=$NAME.tar.gz
    URL=https://github.com/CGAL/cgal/archive/releases/$FILE

    TEST1=$PATH_Ext/lib/libCGAL.*
    TEST2=$PATH_Ext/lib64/libCGAL.*

    _test2 $TEST1 $TEST2
    if [[ "$found" == 1 ]]; then
        echo "   > ($NAME) is already installed at the requested location ($PATH_Ext)."

    else
        cd $PATH_Ext/downloads
        echo "   > Installing ($NAME)"

        # fetch the code base
        _fetch $URL $FILE $PATH_Ext/$NAME.download.log

        # configure step is custom for each package
        BUILD_PATH=$dirname/build
        mkdir -p $BUILD_PATH; cd $BUILD_PATH
        echo "     > Configuring out-of-source build (`pwd`)"

        rm $PATH_Ext/$NAME.cmake.log 2>/dev/null

        cmake -DCMAKE_INSTALL_PREFIX:STRING=$PATH_Ext \
              -DCMAKE_BUILD_TYPE:STRING=Release \
              -DCMAKE_CXX_FLAGS:STRING="-Wno-dev -Wno-unknown-warning-option" \
              -DBUILD_SHARED_LIBS:BOOL=ON \
              -DWITH_CGAL_Qt5:BOOL=OFF \
              -DWITH_GMP:BOOL=OFF \
              -DWITH_MPFR:BOOL=OFF \
              .. > $PATH_Ext/$NAME.cmake.log

        # common functionality to build and install
        _build $NAME $PATH_Ext/$NAME.make.log

        # test the installation
        _test2 $TEST1 $TEST2
        if [[ "$found" == 1 ]]; then
            echo "   > ($NAME) successfully installed at ($PATH_Ext)."
        else
            echo "   > ($NAME) installation failed. Please see build logs in ($PATH_Ext) for more information."
        fi
    fi
fi

# ------------------------------------------------------------------------------
# VTK 8.1.2
# ------------------------------------------------------------------------------
if [ "$INSTALL_VTK" = true ] ; then

    # install vtk 8.1.2
    VERSION='8.1'
    NAME='VTK-'$VERSION'.2'
    FILE=$NAME.tar.gz
    URL=https://www.vtk.org/files/release/$VERSION/$FILE
    TEST=$PATH_Ext/lib/libvtkCommonCore-$VERSION*

    _test $TEST
    if [[ "$found" == 1 ]]; then
        echo "   > ($NAME) is already installed at the requested location ($PATH_Ext)."

    else
        cd $PATH_Ext/downloads
        echo "   > Installing ($NAME)"

        # fetch the code base
        _fetch $URL $FILE $PATH_Ext/$NAME.download.log

        # configure step is custom for each package
        BUILD_PATH=$dirname/build
        mkdir -p $BUILD_PATH; cd $BUILD_PATH
        echo "     > Configuring out-of-source build (`pwd`)"

        rm $PATH_Ext/$NAME.cmake.log 2>/dev/null

        cmake -DCMAKE_INSTALL_PREFIX:STRING=$PATH_Ext \
              -DCMAKE_BUILD_TYPE:STRING=Release \
              -DBUILD_SHARED_LIBS:BOOL=ON \
              -DCMAKE_CXX_FLAGS:STRING="-Wno-inconsistent-missing-override -Wno-deprecated-declarations -Wno-dev -Wno-unknown-warning-option" \
              -DPYTHON_EXECUTABLE=`which python3` \
              -DVTK_PYTHON_VERSION=3 \
              -DVTK_WRAP_PYTHON:BOOL=ON \
              -DVTK_Group_Rendering:BOOL=OFF \
              -DVTK_Group_StandAlone:BOOL=OFF \
              -DModule_vtkCommonDataModel:BOOL=ON \
              -DModule_vtkFiltersGeneral:BOOL=ON \
              -DModule_vtkIOXML:BOOL=ON \
              .. > $PATH_Ext/$NAME.cmake.log

        # common functionality to build and install
        _build $NAME $PATH_Ext/$NAME.make.log

        # test the installation
        _test $TEST
        if [[ "$found" == 1 ]]; then
            echo "   > ($NAME) successfully installed at ($PATH_Ext)."
        else
            echo "   > ($NAME) installation failed. Please see build logs in ($PATH_Ext) for more information."
        fi
    fi
fi

# ------------------------------------------------------------------------------
# end of the install script
# ------------------------------------------------------------------------------
