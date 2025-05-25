FROM docker.io/nvidia/cuda:12.6.2-devel-ubuntu24.04

ENV DEBIAN_FRONTEND=noninteractive \
    LIBMESH_DIR=/usr/local \
    OMPI_MCA_btl=^openib \
    PATH=$PATH:/opt/petsc/arch-linux-c-opt/bin:/usr/local/cuda/bin:/opt/paraview/ParaView-5.13.0-RC1-MPI-Linux-Python3.10-x86_64/bin \
    PETSC_ARCH=arch-linux-c-opt \
    PETSC_DIR=/opt/petsc \
    SLEPC_DIR=/opt/petsc/arch-linux-c-opt \
    __SINGTAG=CHIMAS-ENV

RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        curl \
        dbus-x11 \
        eog \
        exuberant-ctags \
        fonts-freefont-otf \
        freeglut3-dev \
        g++ \
        gcc \
        gdb \
        gfortran \
        git \
        ksh \
        less \
        libblas-dev \
        libboost-all-dev \
        libcairo2-dev \
        libcgal-dev \
        libfltk1.3-dev \
        libfreetype-dev \
        libglpk-dev \
        libgmp-dev \
        liblapack-dev \
        libmetis-dev \
        libmpfr-dev \
        libnuma-dev \
        libocct-*-dev \
        libopenmpi-dev \
        libpthread-stubs0-dev \
        libtbb-dev \
        libxcb-xinerama0 \
        locales \
        m4 \
        make \
        octave \
        python3 \
        python3-pip \
        python3-venv \
        qt5-qmake \
        qtbase5-dev \
        rapidjson-dev \
        rsync \
        shared-mime-info \
        sudo \
        texlive-full \
        ttf-mscorefonts-installer \
        unzip \
        valgrind \
        valgrind-mpi \
        vim \
        wget \
        xclip \
        xterm \
        zip && \
    rm -rf /var/lib/apt/lists/*

RUN locale-gen en_US en_US.UTF-8 pt_BR.UTF-8

ENV LC_ALL=en_US.UTF-8

COPY download/* Makefile.* /opt

RUN make -f /opt/Makefile.petsc N=10

RUN make -f /opt/Makefile.vtk N=10

RUN make -f /opt/Makefile.gmsh N=10

RUN make -f /opt/Makefile.paraview N=10

RUN make -f /opt/Makefile.libmesh N=10

RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        cabextract \
        fontforge && \
    rm -rf /var/lib/apt/lists/*

COPY ttf-vista-fonts-installer.sh /tmp

RUN /tmp/ttf-vista-fonts-installer.sh -q -0 -

RUN rm -f /opt/*gz /opt/Makefile.*

RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        libxcb-cursor0 \
        libxcb-icccm4 \
        libxcb-image0 \
        libxcb-keysyms1 \
        libxcb-randr0 \
        libxcb-render-util0 \
        libxcb-xfixes0 \
        libxcb-xinerama0 \
        libxkbcommon-x11-0 \
        qtwayland5 && \
    rm -rf /var/lib/apt/lists/*

RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        zathura && \
    rm -rf /var/lib/apt/lists/*


