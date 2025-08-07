#SINGULARITY_CACHEDIR := /raid/bfq9/tmp/singularity_${USER}
SINGULARITY_CACHEDIR := /tmp/singularity_${USER}
SINGULARITY_TMPDIR := $(SINGULARITY_CACHEDIR)
TMPDIR := /tmp
export SINGULARITY_CACHEDIR
export SINGULARITY_TMPDIR
export TMPDIR


BUILD_DIR = build
INSTALL_DIR = /opt/container/
TAG = kraken
DOCKERFILE = Dockerfile
TARFILE = $(BUILD_DIR)/$(TAG).tar
SIFFILE = $(TAG).sif

all: .force build install

build: .force
	rm -rf $(BUILD_DIR)
	mkdir -p $(BUILD_DIR) $(INSTALL_DIR)
	podman build -f $(DOCKERFILE) -t $(TAG) .
	rm -f $(TARFILE)
	podman save $(TAG) -o $(TARFILE)
	singularity build --force $(BUILD_DIR)/$(SIFFILE) docker-archive:$(TARFILE)

install: .force
	sudo mv $(BUILD_DIR)/$(SIFFILE) $(INSTALL_DIR)
	sudo chmod a+rx $(INSTALL_DIR)/$(SIFFILE)

.force:
