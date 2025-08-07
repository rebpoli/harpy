NPROC = 10

Stage0 += baseimage(image='docker.io/nvidia/cuda:12.6.2-devel-ubuntu24.04')
# Stage0 += baseimage(image='docker.io/nvidia/cuda:12.2.2-devel-ubuntu22.04')
# Stage0 += baseimage(image='docker.io/nvidia/cuda:11.1.1-cudnn8-devel-ubuntu20.04')

## Ambiente
Stage0 += environment(
    variables={
        '__SINGTAG'         :'CHIMAS-ENV', # tag para identificar o container

        'PATH'              :'$PATH:/opt/petsc/arch-linux-c-opt/bin:/usr/local/cuda/bin:/opt/paraview/ParaView-5.13.0-RC1-MPI-Linux-Python3.10-x86_64/bin',

        'PETSC_DIR'         : '/opt/petsc',
        'PETSC_ARCH'        : 'arch-linux-c-opt',
        'SLEPC_DIR'         :'/opt/petsc/arch-linux-c-opt',

        'LIBMESH_DIR'       :'/usr/local',

        ## Evita mensagens do OpenMPI: https://github.com/open-mpi/ompi/issues/10693
        'OMPI_MCA_btl'      :'^openib',
        'DEBIAN_FRONTEND'   : 'noninteractive'

    }, _export=False)


## Pacotes
Stage0 += packages(
    ospackages=[
        'cmake',
        'gcc',
        'g++',
        'gfortran',
        'valgrind',
        'valgrind-mpi',
        'make',
        'wget',
        'ksh',
        'm4',
        'git',

        'octave',

        'build-essential',
        'python3',
        'python3-pip' ,
        'python3-venv' ,

        'libopenmpi-dev',
        'libboost-all-dev',

        'libcairo2-dev',
        'freeglut3-dev',
        'libblas-dev',
        'liblapack-dev',
        'libtbb-dev',
        'libmpfr-dev',
        'libglpk-dev',
        'libnuma-dev',
        'libpthread-stubs0-dev',
        'libgmp-dev',
        'libfreetype-dev',
        'libfltk1.3-dev',
        'libocct-*-dev',

        'rapidjson-dev',
        'vim',
        'xterm',
        'zip', 'less', 'unzip', 'gdb', 'rsync' , 'sudo',

        # Para funcionar o EOG (eye of gnome)
        'eog', 'shared-mime-info', 'dbus-x11',

        'locales',   # para rodar o locale-gen

        # Para funcionar fontes no octave
        'qtbase5-dev', 'qt5-qmake',
        'fonts-freefont-otf',
        'libcgal-dev',

        # Extras
        'exuberant-ctags',
        'curl', 
        'xclip',
        'libxcb-xinerama0',
        'libmetis-dev',

        # Fontes microsoft
        'ttf-mscorefonts-installer',

        # latex
        'texlive-full'
    ])

## ajusta linguagem para acertar acentos no bash
Stage0 += shell(commands=[f"locale-gen en_US en_US.UTF-8 pt_BR.UTF-8"], _args=False, chdir=False)
Stage0 += environment( variables={ 'LC_ALL':'en_US.UTF-8' }, _export=False)  # 'LANG':'en_US.UTF-8'

# Instala pacotes que requerem compilacao, download etc
Stage0 += copy(src=f"download/* Makefile.*", dest=f"/opt")
for i in ["petsc", "vtk", "gmsh", "paraview", "libmesh"] :
    Stage0 += shell(commands=[f"make -f /opt/Makefile.{i} N={NPROC}"], _args=False, chdir=False)

# Install Fonte Cambria math (generate graphs with the same MS for used in MSWord equations)
Stage0 += packages( ospackages=[ 'cabextract', 'fontforge' ] )
Stage0 += copy(src="ttf-vista-fonts-installer.sh", dest="/tmp")
Stage0 += shell(commands=[f"/tmp/ttf-vista-fonts-installer.sh -q -0 -"], _args=False, chdir=False)

# Cleanup  - it seems that podman does not like the caches after cleanups. So keep it last
Stage0 += shell(commands=[f"rm -f /opt/*gz /opt/Makefile.*"], _args=False, chdir=False)

# So that python matplotlib works
Stage0 += packages( ospackages=[ 'qtwayland5', 'libxkbcommon-x11-0', 'libxcb-cursor0', 'libxcb-icccm4',
                                 'libxcb-image0', 'libxcb-keysyms1', 'libxcb-randr0', 'libxcb-render-util0',
                                 'libxcb-xinerama0', 'libxcb-xfixes0' ] )

# For tex edition in VIM
Stage0 += packages( ospackages=[ 'zathura' ] )
#  For autodiff
Stage0 += packages( ospackages=[ 'libeigen3-dev' ] )   

Stage0 += packages( ospackages=[ 'xauth' ] )   

Stage0 += environment( variables={ 'SINGTAG' :'HARPY' } )
