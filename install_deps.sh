#!/bin/bash

# ------------------------------------------------------------------------------
# parameters to control this script!

INSTALL_CGAL=true
INSTALL_BOOST=true
INSTALL_EIGEN=true
NPROCS=20


if [ -n "$BASH" ] ;then
    SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
elif [ -n "$ZSH_NAME" ] ;then
    SCRIPT_DIR="$(realpath `dirname ${0}`)"
fi
PATH_Ext=$SCRIPT_DIR/external

# ------------------------------------------------------------------------------
# install script for MemSurfer and its dependencies
# ------------------------------------------------------------------------------
echo ''
echo "> Installing MemSurfer's dependencies for (`whoami`) on (`hostname`). platform = (`uname`)"
echo "    > Installation path = ("$PATH_Ext")"

unameout=`uname`
if [ "$unameout" == "Linux" ]; then
    shlib_extn=".so"
elif [ "$unameout" == "Darwin" ]; then
    shlib_extn=".dylib"
else
    echo "Script works only for Linux and Darwin"
    exit
fi

# ------------------------------------------------------------------------------
# check commands
# ------------------------------------------------------------------------------
command -v wget >/dev/null 2>&1 || { echo >&2 "Cannot find command 'wget'. Aborting."; exit 1;  }
command -v tar >/dev/null 2>&1 || { echo >&2 "Cannot find command 'tar'. Aborting."; exit 1;    }
command -v cmake >/dev/null 2>&1 || { echo >&2 "Cannot find command 'cmake'. Aborting."; exit 1;  }
command -v swig >/dev/null 2>&1 || { echo >&2 "Cannot find command 'swig'. Aborting."; exit 1;  }
command -v cython >/dev/null 2>&1 || { echo >&2 "Cannot find command 'cython'. Aborting."; exit 1;  }

if [ -z "${CC_COMPILER}" ] ; then
    export CC_COMPILER=`which gcc`
fi
if [ -z "${CXX_COMPILER}" ] ; then
    export CXX_COMPILER=`which g++`
fi
export CC=${CC_COMPILER}
export CXX=${CXX_COMPILER}
export PYTHON=`which python3`


command -v "${CC}" >/dev/null 2>&1 || { echo >&2 "Cannot find C compiler '${CC}'. Aborting."; exit 1;  }
command -v "${CXX}" >/dev/null 2>&1 || { echo >&2 "Cannot find CXX compiler '${CXX}'. Aborting."; exit 1;  }
command -v "${PYTHON}" >/dev/null 2>&1 || { echo >&2 "Cannot find Python '${PYTHON}'. Aborting."; exit 1;  }

CVERSION=`${CC} --version | head -n1`
CXXVERSION=`${CXX} --version | head -n1`
PYVERSION=`${PYTHON} --version`
echo "    > Using gcc = (${CC}) [${CVERSION}]"
echo "    > Using g++ = (${CXX}) [${CXXVERSION}]"
echo "    > Using python = (${PYTHON}) [${PYVERSION}]"

# ------------------------------------------------------------------------------
# common utilities
# ------------------------------------------------------------------------------
_fetch() {
    modname=$1
    httplink=$2
    logname=$3

    filename=`basename $httplink`
    dirname=${filename%.*}

    rm $logname 2>/dev/null
    if ! ls $filename 1> /dev/null 2>&1 ; then
        echo "    ($modname) Downloading ($httplink)"
        wget -o $logname $httplink

        rm -rf $dirname 2>/dev/null
    fi
    if ! ls $dirname 1> /dev/null 2>&1 ; then
        mkdir $dirname
        echo "    ($modname) Extracting ($dirname)"
        tar -xf $filename -C $dirname --strip-components=1
    fi
}

# ------------------------------------------------------------------------------
mkdir -p $PATH_Ext/downloads
pushd $PATH_Ext/downloads > /dev/null

# ------------------------------------------------------------------------------
# Eigen 3.3.9
# ------------------------------------------------------------------------------
if [ "$INSTALL_EIGEN" = true ] ; then

    VERSION='3.3.9'
    NAME='eigen-3.3.9'
    FILE='eigen-3.3.9.tar.bz2'
    URL="https://gitlab.com/libeigen/eigen/-/archive/$VERSION/$FILE"
    TEST_DIR="$PATH_Ext/include/eigen3"

    if [ -d $TEST_DIR ]; then
      echo "    ($NAME) Already installed. Found ($TEST_DIR)."
    else
      rm $PATH_Ext/$NAME.*.log 2>/dev/null

      ts=$(date "+%Y%m%d-%H%M%S")
      _fetch $NAME $URL $PATH_Ext/$NAME.download-$ts.log

      BUILD_PATH=$dirname/build-$ts
      mkdir -p $BUILD_PATH
      pushd $BUILD_PATH > /dev/null
      echo "    ($NAME) Configuring (`pwd`)"

      cmake -DCMAKE_C_COMPILER=${CC_COMPILER} \
            -DCMAKE_CXX_COMPILER=${CXX_COMPILER} \
            -DCMAKE_INSTALL_PREFIX:STRING=$PATH_Ext \
            -DCMAKE_BUILD_TYPE:STRING=Release \
            .. > $PATH_Ext/$NAME.cmake-$ts.log 2>&1

      # eigen needs only make install
      echo "    ($NAME) Building and installing"
      make install -j$NPROCS >$PATH_Ext/$NAME.make-$ts.log 2>&1

      # test the installation
      if [ -d $TEST_DIR ]; then
        echo "    ($NAME) Successfully installed at ($PATH_Ext)."
      else
        echo "    ($NAME) Installation failed. Please see build logs in ($PATH_Ext) for more information."
      fi
      popd > /dev/null
    fi
fi


# ------------------------------------------------------------------------------
# Boost 1.66
# ------------------------------------------------------------------------------
if [ "$INSTALL_BOOST" = true ] ; then

    VERSION='1.66.0'
    FILE='boost_1_66_0.tar.gz'
    NAME='boost-1.66.0'
    URL="https://boostorg.jfrog.io/artifactory/main/release/$VERSION/source/$FILE"
    TEST_FILE="$PATH_Ext/lib/libboost_graph.a"

    if [ -f $TEST_FILE ]; then
      echo "    ($NAME) Already installed. Found ($TEST_FILE)."
    else
      rm $PATH_Ext/$NAME.*.log 2>/dev/null

      ts=$(date "+%Y%m%d-%H%M%S")
      _fetch $NAME $URL $PATH_Ext/$NAME.download-$ts.log

      BUILD_PATH=$dirname
      mkdir -p $BUILD_PATH
      pushd $BUILD_PATH > /dev/null
      echo "    ($NAME) Configuring (`pwd`)"

      ./bootstrap.sh --prefix=$PATH_Ext \
                     --with-toolset=gcc \
                     --with-python=$PYTHON \
                     --with-libraries=atomic,thread,graph,chrono,date_time \
                     > $PATH_Ext/$NAME.bootstrap-$ts.log 2>&1

      echo "    ($NAME) Building and installing"
      ./b2 install > $PATH_Ext/$NAME.make-$ts.log 2>&1

      # test the installation
      if [ -f $TEST_FILE ]; then
        echo "    ($NAME) Successfully installed at ($PATH_Ext)."
      else
        echo "    ($NAME) Installation failed. Please see build logs in ($PATH_Ext) for more information."
      fi
      popd > /dev/null
    fi
fi


# ------------------------------------------------------------------------------
# CGAL 4.13
# ------------------------------------------------------------------------------
if [ "$INSTALL_CGAL" = true ] ; then

    VERSION='4.13'
    NAME='CGAL-4.13'
    FILE='CGAL-4.13.tar.gz'
    URL="https://github.com/CGAL/cgal/archive/releases/$FILE"
    TEST_FILE1="$PATH_Ext/lib/libCGAL$shlib_extn"
    TEST_FILE2="$PATH_Ext/lib64/libCGAL$shlib_extn"

    if [ -f $TEST_FILE1 ]; then
      echo "    ($NAME) Already installed. Found ($TEST_FILE1)."
    elif [ -f $TEST_FILE2 ]; then
      echo "    ($NAME) Already installed. Found ($TEST_FILE2)."
    else
      rm $PATH_Ext/$NAME.*.log 2>/dev/null

      ts=$(date "+%Y%m%d-%H%M%S")
      _fetch $NAME $URL $PATH_Ext/$NAME.download-$ts.log

      BUILD_PATH=$dirname/build-$ts
      mkdir -p $BUILD_PATH
      pushd $BUILD_PATH > /dev/null
      echo "    ($NAME) Configuring (`pwd`)"

      cmake -DCMAKE_C_COMPILER=${CC_COMPILER} \
            -DCMAKE_CXX_COMPILER=${CXX_COMPILER} \
            -DCMAKE_INSTALL_PREFIX:STRING=$PATH_Ext \
            -DCMAKE_BUILD_TYPE:STRING=Release \
            -DCMAKE_CXX_FLAGS:STRING="-Wno-dev -Wno-unknown-warning-option" \
            -DBUILD_SHARED_LIBS:BOOL=ON \
            -DWITH_CGAL_Qt5:BOOL=OFF \
            -DWITH_GMP:BOOL=OFF \
            -DWITH_MPFR:BOOL=OFF \
            .. > $PATH_Ext/$NAME.cmake-$ts.log 2>&1

      # build and install
      echo "    ($NAME) Building and Installing"
      make -j$NPROCS > $PATH_Ext/$NAME.make-$ts.log 2>&1
      make install -j$NPROCS > $PATH_Ext/$NAME.make-$ts.log 2>&1

      # test the installation
      if [ -f ${TEST_FILE1} -o -f ${TEST_FILE2} ]; then
        echo "    ($NAME) Successfully installed at ($PATH_Ext)."
      else
        echo "    ($NAME) Installation failed. Please see build logs in ($PATH_Ext) for more information."
      fi
      popd > /dev/null
    fi
fi

# ------------------------------------------------------------------------------
# end of the install script
# ------------------------------------------------------------------------------
