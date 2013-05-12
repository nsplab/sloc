CMAKE_BUILD_TYPE = Debug
BUILD_TESTS = 0

CMAKE_FLAGS = \
	-DBOOST_ROOT=$(HOME)/opt/local \
	-DDEALII_PREFIX=$(HOME)/opt/local \
	-DGETFEM_PREFIX=$(HOME)/opt/local \
	-DCMAKE_BUILD_TYPE="$(CMAKE_BUILD_TYPE)" \
	-DBUILD_TESTS="${BUILD_TESTS}"

default: build

build:
	mkdir -p build && cd build && cmake $(CMAKE_FLAGS) .. && make

clean:
	rm -rf ./build

.PHONY: build

