# CMakeLists.txt
#
# Things to build:
# (1) libsloc.a library (for common code)
# (2) programs in "bin" directory
# (3) programs in "tests" directory
#

project(sloc)
cmake_minimum_required(VERSION 2.8)

# Set path for CMake dependencies
if(NOT CMAKE_MODULE_PATH)
    set(CMAKE_MODULE_PATH) # init to empty list
endif()
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/config")

# Build out a list of common libraries (initialize to empty list)
set(COMMON_LIBRARIES)

# Decide whether to build the tests or not (default is 'false')
if(NOT DEFINED BUILD_TESTS)
    set(BUILD_TESTS FALSE)
endif()

# Set a default build type for single-configuration
# CMake generators if no build type is set.
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

# Select flags
set(CMAKE_CXX_FLAGS "-Wall -Wno-unused-variable")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -ltbb")
set(CMAKE_CXX_FLAGS_DEBUG  "-O0 -g -ltbb")

# -----------------------------------------------------------------------------

# Configure MPI
find_package(MPI REQUIRED)
if(MPI_CXX_FOUND)
    #set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
    set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS})
    set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS})
    include_directories(${MPI_CXX_INCLUDE_PATH})
    list(APPEND COMMON_LIBRARIES ${MPI_CXX_LIBRARIES})

    message(STATUS "MPI settings:")
    message(STATUS "  MPI_CXX_COMPILER = ${MPI_CXX_COMPILER}")
    message(STATUS "  MPI_CXX_COMPILE_FLAGS = ${MPI_CXX_COMPILE_FLAGS}")
    message(STATUS "  MPI_CXX_INCLUDE_PATH = ${MPI_CXX_INCLUDE_PATH}")
    message(STATUS "  MPI_CXX_LIBRARIES = ${MPI_CXX_LIBRARIES}")
    message(STATUS "  CMAKE_CXX_COMPILE_FLAGS = ${CMAKE_CXX_COMPILE_FLAGS}")
    message(STATUS "  CMAKE_CXX_LINK_FLAGS = ${CMAKE_CXX_LINK_FLAGS}")
endif()

# Configure Boost
set(BOOST_USE_STATIC_LIBS ON)
set(BOOST_USE_MULTITHREADED ON)
set(REQUIRED_BOOST_VERSION 1.49)
find_package(Boost ${REQUIRED_BOOST_VERSION}
             COMPONENTS
             exception
             system
             thread
             serialization
             filesystem
             program_options
             REQUIRED
)
include_directories(${Boost_INCLUDE_DIR})
list(APPEND COMMON_LIBRARIES ${Boost_LIBRARIES})
message(STATUS "Boost settings:")
message(STATUS "  BOOST_ROOT = ${BOOST_ROOT}")

# Configure GetFEM
find_package(GetFEM REQUIRED)
add_definitions(${GETFEM_CFLAGS})
add_definitions("-Wno-return-type-c-linkage")
add_definitions("-Wno-deprecated-declarations")
list(APPEND COMMON_LIBRARIES ${GETFEM_LIBS})

# Configure dealII
find_package(DealII)
if(DEALII_FOUND)
    add_definitions(${DEALII_CFLAGS})
    list(APPEND COMMON_LIBRARIES ${DEALII_LDFLAGS})
endif()

# Configure VCGLib
set(VCGLIB_DIR "$ENV{HOME}/dev/vcglib")
include_directories(${VCGLIB_DIR})

# Configure pkg-config (used for configuring other libraries: see "pkg-config --list-all")
find_package(PkgConfig)

# Configure MuParser
pkg_check_modules(MUPARSER muparser)
include_directories(${MUPARSER_INCLUDE_DIRS})
list(APPEND COMMON_LIBRARIES ${MUPARSER_LDFLAGS})
message(STATUS "Found these muParser settings:")
message(STATUS "  MUPARSER_INCLUDE_DIRS = ${MUPARSER_INCLUDE_DIRS}")
message(STATUS "  MUPARSER_LDFLAGS = ${MUPARSER_LDFLAGS}")

# See also http://www.cmake.org/Wiki/CMake_Useful_Variables
message(STATUS "Other project settings:")
message(STATUS "  PROJECT_SOURCE_DIR = ${PROJECT_SOURCE_DIR}")
message(STATUS "  CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")

# -----------------------------------------------------------------------------

# Load build targets and list of libsloc.a sources
include("targets.txt")

# Build libsloc.a
include_directories("${PROJECT_SOURCE_DIR}/include")
add_library(sloc STATIC ${LIBSLOC_SOURCES}) 

# Build stuff in "bin"
foreach(target ${BIN_TARGETS})
    add_executable(${target} "bin/${target}.cc")
    target_link_libraries(${target} sloc ${COMMON_LIBRARIES})
endforeach()

# Build stuff in "tests"
if(BUILD_TESTS)
    foreach(target ${TEST_TARGETS})
        add_executable(${target} "tests/${target}.cc")
        target_link_libraries(${target} sloc ${COMMON_LIBRARIES})
    endforeach()
endif()

# EOF
