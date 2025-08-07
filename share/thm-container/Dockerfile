FROM docker.io/ubuntu:24.04

# FROM docker.io/ecpe4s/e4s-gpu:22.02
# FROM nvcr.io/nvidia/nvhpc:23.3-devel-cuda11.8-ubuntu22.04
# FROM nvcr.io/nvidia/nvhpc:25.5-devel-cuda12.9-ubuntu24.04
# FROM docker.io/nvidia/cuda:12.9.0-devel-ubuntu24.04
# FROM nvcr.io/nvidia/nvhpc:25.3-devel-cuda_multi-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV SINGTAG=KRAKEN
ENV NVIDIA_VISIBLE_DEVICES=all
ENV NVIDIA_DRIVER_CAPABILITIES=compute,utility

# Install necessary build tools and dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    curl \
    gfortran \
    git \
    wget \
    python3 \
    python3-dev \
    python3-pip \
    python3-venv \
    gfortran \
    libomp-dev \    
    zlib1g-dev \
    m4 \
    automake \
    libtool \
    vim \
    rsync \
    eog \
    exuberant-ctags \
    binutils \
    pkg-config \
    libboost-all-dev \
    libblas-dev liblapack-dev \
    openmpi-bin openmpi-doc libopenmpi-dev \
    libnetcdf-dev netcdf-bin \
    && rm -rf /var/lib/apt/lists/*

# Create directory for installations
WORKDIR /opt/src

#
# INSTALL PETSc
#
RUN wget https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-3.23.5.tar.gz && \
    tar -xzf petsc-3.23.5.tar.gz 

RUN cd petsc-3.23.5 &&\
    ./configure --prefix=/usr/local/petsc \
        --PETSC_ARCH=arch-linux-c-opt \
        --with-debugging=0 \
        --with-mpi \
        --download-hypre \
        --with-blaslapack \
         && \
    make all && \
    make install && \
    cd .. 

ENV PETSC_DIR="/usr/local/petsc" \
    PETSC_ARCH="arch-linux-c-opt" \
    LD_LIBRARY_PATH="/usr/local/petsc/lib:${LD_LIBRARY_PATH}" 

#
# INSTALL Trilinos
#
RUN git clone https://github.com/trilinos/Trilinos.git 
RUN cd Trilinos && git checkout trilinos-release-14-4-0

#     cd Trilinos && git checkout trilinos-release-16-1-0  # latest - was not compiling
RUN apt-get update && apt-get install -y libmatio-dev xorg-dev

RUN cd Trilinos && mkdir build && cd build && cmake  \
    -DTrilinos_ENABLE_Amesos=ON                      \
    -DTrilinos_ENABLE_Epetra=ON                      \
    -DTrilinos_ENABLE_EpetraExt=ON                   \
    -DTrilinos_ENABLE_Ifpack=ON                      \
    -DTrilinos_ENABLE_AztecOO=ON                     \
    -DTrilinos_ENABLE_Sacado=ON                      \
    -DTrilinos_ENABLE_SEACAS=ON                      \
    -DTrilinos_ENABLE_Teuchos=ON                     \
    -DTrilinos_ENABLE_MueLu=ON                       \
    -DTrilinos_ENABLE_ML=ON                          \
    -DTrilinos_ENABLE_NOX=ON                         \
    -DTrilinos_ENABLE_ROL=ON                         \
    -DTrilinos_ENABLE_Tpetra=ON                      \
    -DTrilinos_ENABLE_COMPLEX=ON                     \
    -DTrilinos_ENABLE_FLOAT=ON                       \
    -DTrilinos_ENABLE_Zoltan=ON                      \
    -DTrilinos_VERBOSE_CONFIGURE=OFF                 \
    -DTPL_ENABLE_MPI=ON                              \
    -DBUILD_SHARED_LIBS=ON                           \
    -DCMAKE_VERBOSE_MAKEFILE=OFF                     \
    -DCMAKE_BUILD_TYPE=RELEASE                       \
    -DCMAKE_INSTALL_PREFIX:PATH=/opt/trilinos        \
    -DTrilinos_ENABLE_EXPLICIT_INSTANTIATION=ON      \  
    ../         && \
    make -j10 && make install

# 
# INSTALL P4est
#
RUN wget https://p4est.github.io/release/p4est-2.8.7.tar.gz && \
    tar zxvf p4est-*gz

RUN apt-get update && apt-get install -y libjansson-dev 

RUN cd p4est-*/ && \
    ./configure --enable-mpi --enable-shared \
        --disable-vtk-binary --without-blas \
        CFLAGS="-O2" CPPFLAGS="-DSC_LOG_PRIORITY=SC_LP_ESSENTIAL" \
        --prefix=/opt/p4est && \ 
    make -j10 && make install

#
# INSTALL Deal II
#
RUN wget https://dealii.org/downloads/dealii-9.7.0.tar.gz && tar xvf dealii-*gz

RUN apt-get update && apt-get install -y libgmsh-dev libmumps-dev

# Ref: https://dealii.org/developer/users/cmake_dealii.html
RUN cd dealii-*/ && mkdir build && cd build && \
    cmake -DWITH_MPI=ON \
          -DWITH_GMSH=ON \
          -DWITH_LAPACK=ON \
          -DWITH_MUMPS=ON \
          -DWITH_P4EST=ON \
          -DWITH_PETSC=ON \
          -DWITH_TRILINOS=ON \
          -DWITH_ZLIB=ON \
          -DP4EST_DIR="/opt/p4est" \
    .. && \
    make -j10 && make install


# Set default command
CMD ["/bin/bash"]

#
# To configure a python environment (if necessary for the container construction)

# # Install netcdf
# ENV VENV=/opt/venv
# RUN python3 -m venv $VENV
# ENV PATH="$VENV/bin:$PATH"
# RUN pip3 install netcdf4

