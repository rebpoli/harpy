
### DEFAULT HEADER
toprel := $(shell i=. ; for c in 0 1 2 3 4 5 6 7 8; do [ ! -f $$i/.chimas_top ] && i=$$i/..; done; echo $$i)
topabs := $(realpath $(toprel))
###

CHDEBUG=10

all: build/opt build/devel build/dbg

build/%: .force cmake-%
	$(MAKE) -C $@

cmake: cmake-opt cmake-devel cmake-dbg
cmake-%: 
	mkdir -p build/$*
	echo $(topabs)
	cd build/$* && export METHOD=$* && cmake $(topabs)

tags: .force
	ctags -R --links=no include src /opt/libmesh/include /opt/libmesh/src

clean:
	rm -rf build

# Count number of lines of code
cloc:
	find -type f -name "*.txt" > exclude.txt
	find -type f -name "*.e" >> exclude.txt
	find -type f -name "*.csv" >> exclude.txt
	find -type f -name "*.json" >> exclude.txt
	find -type f -name "*.svg" >> exclude.txt
	find build/* >> exclude.txt
	cloc --exclude-list-file exclude.txt .
	rm exclude.txt

.force:
