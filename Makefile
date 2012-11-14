
CMAKE_FLAGS = \
	-DBOOST_ROOT=$(HOME)/opt/local \
	-DDEALII_PREFIX=$(HOME)/dev/deal.II \
	-DGETFEM_PREFIX=$(HOME)/opt/getfem \
	-DBUILD_TESTS=1

default: build

build:
	mkdir -p build && cd build && cmake $(CMAKE_FLAGS) .. && make

clean:
	rm -rf ./build

.PHONY: build

