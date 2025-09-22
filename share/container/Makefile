
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
# SINGULARITY_CACHEDIR := /raid/bfq9/tmp/singularity_${USER}
# BASEDIR = /raid/bfq9

# WSL
# SINGULARITY_CACHEDIR := /home/renato/temp/singularity_cache
# BASEDIR = /home/renato/build

#Alien
SINGULARITY_CACHEDIR := /tmp/sing_cache
BASEDIR = build/

SINGULARITY_TMPDIR := $(SINGULARITY_CACHEDIR)
export SINGULARITY_CACHEDIR
export SINGULARITY_TMPDIR
# TMPDIR := /mnt/data/tmp/
TMPDIR := /tmp
export TMPDIR

MASK := harpy
BUILDDIR := $(BASEDIR)/tmp/build/$(MASK)
INSTALLDIR := $(BASEDIR)/sif/$(MASK)

TAG := $(MASK)_$(USER)
PYFILE := $(MASK).py
DOCKERFILE := $(BUILDDIR)/Dockerfile_$(MASK)
TARFILE := $(BUILDDIR)/$(MASK).tar
SIFFILE := $(INSTALLDIR)/$(MASK).sif

all: build
download: .force
	rm -rf download
	mkdir -p download
	cp Makefile.download download/Makefile
	$(MAKE) -C download

build: $(TARFILE)
	@echo "**"
	@echo "** Rodando singularity (Docker tar => sif)"
	@echo "**"
	mkdir -p $(INSTALLDIR)
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
	mkdir -p $(BUILDDIR)
	hpccm --recipe $(PYFILE) --format docker > $(DOCKERFILE) 

clean:
	rm -rf $(DOCKERFILE) $(TARFILE)

install:
	mv $(SIFFILE) $(INSTALLDIR)
	scp $(SIFFILE) atenadev01:/dfs_geral_ep/res/santos/unbs/gger/er/er01/USR/bfq9/sif/harpy.sif

# push:
#     podman login ghcr.io -u rebpoli
#     podman tag localhost/harpy:latest ghcr.io/rebpoli/harpy:latest
#     podman push ghcr.io/rebpoli/harpy:latest

# pull:
#     podman pull ghcr.io/rebpoli/harpy:latest

push:
	singularity push $(SIFFILE) docker://ghcr.io/rebpoli/harpy:latest

.force:
