
CMAKE_FLAGS = \
	-DGETFEM_PREFIX=$(HOME)/opt/getfem \
	-DDEALII_PREFIX=$(HOME)/dev/deal.II \
	-DBUILD_TESTS=1

default: build

build:
	mkdir -p build && cd build && cmake $(CMAKE_FLAGS) .. && make

clean:
	rm -rf ./build

.PHONY: build

