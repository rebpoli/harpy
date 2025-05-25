#
#
#
BDIR := /opt/cmake-3.27.7
cmake:
	cd /opt && wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 0 https://github.com/Kitware/CMake/releases/download/v3.27.7/cmake-3.27.7.tar.gz && tar -xzf cmake-3.27.7.tar.gz 
	cd $(BDIR) && ./configure && make -j$(N) all && make install 

