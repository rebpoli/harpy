
# Para construir no unix de casa
# SINGULARITY_CACHEDIR := /mnt/data/tmp/singularity_${USER}
# BASEDIR = /mnt/data

# Para construir na UT
# SINGULARITY_CACHEDIR := /tmp/singularity_${USER}
# BASEDIR = ${HOME}/build

# Cuidar com o storage. Configurar aqui:
#    ~/.config/containers/storage.conf
# ou aqui
#    /etc/containers/storage.conf
#
# [storage]
#   driver = "overlay"
#   runroot = "/run/user/1000"
#   graphroot = "/mnt/data/tmp"
#   [storage.options]
#     mount_program = "/bin/fuse-overlayfs"


# Petrobras
SINGULARITY_CACHEDIR := /raid/bfq9/tmp/singularity_${USER}
BASEDIR = /raid/bfq9

# WSL
# SINGULARITY_CACHEDIR := /home/renato/temp/singularity_cache
# BASEDIR = /home/renato/build


SINGULARITY_TMPDIR := $(SINGULARITY_CACHEDIR)
export SINGULARITY_CACHEDIR
export SINGULARITY_TMPDIR
# TMPDIR := /mnt/data/tmp/
TMPDIR := /tmp
export TMPDIR

MASK := lm23_ub24
INSTALLDIR := $(BASEDIR)/build/$(MASK)

TAG := $(MASK)_$(USER)
PYFILE := $(MASK).py
DOCKERFILE := $(INSTALLDIR)/Dockerfile_$(MASK)
TARFILE := $(INSTALLDIR)/$(MASK).tar
SIFFILE := $(INSTALLDIR)/$(MASK).sif

all: build
download: .force
	rm -rf download
	mkdir -p download
	cp Makefile.download download/Makefile
	make -j5 -C download


build: $(TARFILE)
	@echo "**"
	@echo "** Rodando singularity (Docker tar => sif)"
	@echo "**"
	singularity build --force $(SIFFILE) docker-archive:$(TARFILE)

$(TARFILE): $(DOCKERFILE)
	@echo "**"
	@echo "** Rodando Docker (Dockerfile => Docker tar)"
	@echo "**"
	podman build -f $(DOCKERFILE) -t $(TAG) .
	rm -f $(TARFILE)
	podman save $(TAG) -o $(TARFILE)

$(DOCKERFILE): $(PYFILE) .force
	@echo "**"
	@echo "** Rodando hpccm (py => Dockerfile)"
	@echo "**"
	mkdir -p $(INSTALLDIR)
	hpccm --recipe $(PYFILE) --format docker > $(DOCKERFILE) 

clean:
	rm -rf $(DOCKERFILE) $(TARFILE)

install:
	mv $(SIFFILE) $(INSTALLDIR)/../../sif

sync_sif_to_sdu:
	scp /raid/bfq9/build/lm23/lm23.sif lncc:/petrobr/poromec/sif

.force:
