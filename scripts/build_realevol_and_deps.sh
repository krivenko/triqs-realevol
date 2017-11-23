#!/bin/sh

export CC=gcc                               # C compiler to build TRIQS/realevol
export CXX=g++                              # C++ compiler to build TRIQS/realevol
export FC=gfortran                          # FORTRAN compiler to build ARPACK-NG
BUILD_TYPE=Debug                            # CMake build type (Release, Debug, RelWithDebInfo)
INSTALL_DIR=$(pwd)/installed                # Installation directory for realevol and prereqs
REALEVOL_SRC_DIR=$(pwd)/realevol.git        # Directory with realevol sources
BUILD_DIR=$(pwd)/build                      # Directory to build realevol and its prerequisites in
PYTHON_INTERPRETER=/usr/bin/python2.7       # Path to Python interpreter
PYTHON_LIBRARY=/usr/lib64/libpython2.7.so   # Path to Python shared library

set -x
set -e

#
# Create/clean installtion directory
#

mkdir -p ${INSTALL_DIR}
[ -e ${INSTALL_DIR} ] && rm -rf ${INSTALL_DIR}/*

#
# Create build directory
#
mkdir -p $BUILD_DIR
pushd $BUILD_DIR

#
# Install TRIQS library
#

if [ ! -d triqs.git ]; then
  git clone https://github.com/TRIQS/triqs.git triqs.git
fi
pushd triqs.git && git checkout 1.4.1 && popd
mkdir -p triqs.build
pushd triqs.build
rm -fv CMakeCache.txt
rm -frv CMakeFiles
cmake ../triqs.git \
  -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
  -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
  -DPYTHON_INTERPRETER=${PYTHON_INTERPRETER} \
  -DPYTHON_LIBRARY=${PYTHON_LIBRARY} \
  -DBuild_Documentation=OFF \
  -DBuild_Tests=ON
make -j4 && make test && make install
popd

#
# Install ARPACK-NG
#

if [ -d arpack-ng.git ]; then
  pushd arpack-ng.git && git checkout master && git pull && popd
else
  git clone -b master https://github.com/opencollab/arpack-ng.git arpack-ng.git
fi
pushd arpack-ng.git
rm -fv CMakeCache.txt
rm -frv CMakeFiles
cmake . \
  -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
  -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
  -DENABLE_STATIC=OFF
make -j4 && make check && make test && make install
popd

#
# Install triqs_arpack
#

if [ ! -d triqs_arpack.git ]; then
  git clone https://github.com/krivenko/triqs_arpack.git triqs_arpack.git
fi
pushd triqs_arpack.git && git checkout 0.5 && popd
pushd triqs_arpack.git
rm -fv CMakeCache.txt
rm -frv CMakeFiles
cmake . \
  -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
  -DTRIQS_PATH=${INSTALL_DIR} \
  -DARPACK_DIR=${INSTALL_DIR}
make -j4 && make test && make install
popd

#
# Build realevol
#

mkdir -p realevol.build
pushd realevol.build
rm -fv CMakeCache.txt
rm -frv CMakeFiles
cmake ${REALEVOL_SRC_DIR} \
  -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
  -DTRIQS_PATH=${INSTALL_DIR} \
  -DARPACK_DIR=${INSTALL_DIR} \
  -DTests=ON
make -j4 && make test
popd

popd
