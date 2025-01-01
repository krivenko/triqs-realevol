# syntax=docker/dockerfile:1
FROM flatironinstitute/triqs:3.1.0
LABEL maintainer="Igor Krivenko"
LABEL description="Real time evolution solver based on TRIQS"
LABEL version="0.10.0"

USER root
RUN useradd -m -s /bin/bash -u 999 build && echo "build:build" | chpasswd
RUN usermod -aG sudo build
RUN apt-get update && \
    apt-get install -y --no-install-recommends make g++-10 apt-utils

ENV OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

# Download a C++20 compatible version of Boost
ARG BOOST_URL=https://archives.boost.io/release/1.78.0/source/boost_1_78_0.tar.bz2
RUN curl -O -L $BOOST_URL &&             \
    tar -xf boost_1_78_0.tar.bz2 &&       \
    mv boost_1_78_0 /home/build/boost && \
    rm boost_1_78_0.tar.bz2

# Install ARPACK-NG
USER build
WORKDIR /home/build
RUN git clone https://github.com/opencollab/arpack-ng.git arpack-ng.git
RUN mkdir arpack-ng.build
WORKDIR /home/build/arpack-ng.build
RUN cmake ../arpack-ng.git                              \
        -DCMAKE_INSTALL_PREFIX=/usr                     \
        -DCMAKE_BUILD_TYPE=Release                      \
        -DBUILD_SHARED_LIBS=ON                          \
        -DICB=ON                                        \
        -DMPI=ON
RUN make -j6 && ctest --output-on-failure
USER root
RUN make install

# Install realevol
USER build
WORKDIR /home
COPY --chown=build . /home/build/realevol.git
RUN mkdir /home/build/realevol.build
WORKDIR /home/build/realevol.build
ENV CC=gcc-10 CXX=g++-10 REPO=/build/repo
RUN cmake ../realevol.git                               \
        -DCMAKE_INSTALL_PREFIX=/usr                     \
        -DCMAKE_BUILD_TYPE=Release                      \
        -DBoost_INCLUDE_DIR=$HOME/boost                 \
        -DBUILD_SHARED_LIBS=ON                          \
        -DBuild_Tests=ON                                \
        -DBuild_Benchmarks=OFF                          \
        -DBUILD_DEBIAN_PACKAGE=ON
RUN make -j6 VERBOSE=1 && ctest --output-on-failure && cpack
USER root
RUN make install && mkdir -p $REPO && mv *.deb $REPO

# Cleanup build files
USER root
RUN rm -rf /home/build/boost /home/build/arpack-ng.* /home/build/realevol.*
