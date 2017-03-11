# States that CMake required version must be greater than 2.6
cmake_minimum_required(VERSION 2.6) 
project("pbpic")

enable_language (Fortran)

# convert relative paths to absolute paths 
cmake_policy(SET CMP0015 NEW)

# MPI libraries
find_package(MPI REQUIRED)
add_definitions(${MPI_Fortran_COMPILE_FLAGS})
include_directories(${MPI_Fortran_INCLUDE_PATH})
link_directories(${MPI_Fortran_LIBRARIES})

set(HDF5_USE_STATIC_LIBRARIES ON)
find_package(HDF5 COMPONENTS Fortran REQUIRED)
#find_path(HDF5_INCLUDE include/hdf5.mod HINTS ${HDF5_ROOT})
add_definitions(${HDF5_DEFINITIONS})
if(DEFINED ENV{HDF5_INCL})
    include_directories($ENV{HDF5_INCL})
    MESSAGE($ENV{HDF5_INCL})
elseif(DEFINED ENV{TACC_HDF5_INC})
    include_directories($ENV{TACC_HDF5_INC})
    MESSAGE($ENV{TACC_HDF5_INC})
endif()
include_directories(${HDF5_INCLUDE_DIRS})
link_directories(${HDF5_Fortran_LIBRARIES})
MESSAGE(${HDF5_INCLUDE_DIRS})
MESSAGE(${HDF5_Fortran_LIBRARIES})

#find_package(FFTW)
if(DEFINED ENV{FFTW_INCLUDES})
    include_directories($ENV{FFTW_INCLUDES})
    MESSAGE($ENV{FFTW_INCLUDES})
endif()
if(DEFINED ENV{FFTW_LIBS})
    link_directories($ENV{FFTW_LIBS})
    MESSAGE($ENV{FFTW_LIBS})
endif()

# include_directories(${CMAKE_SOURCE_DIR}/src/third_party/FLAP/static/mod)
# link_directories(${CMAKE_SOURCE_DIR}/src/third_party/FLAP/static)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  # set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
  # set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -fimplicit-none \
    -Wline-truncation -fwhole-file -std=f2008 -fall-intrinsics -fPIC -cpp")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fimplicit-none -Wall \
    -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing \
    -Wimplicit-interface -Wunused-parameter -fwhole-file -fcheck=all \
    -std=f2008 -pedantic -fbacktrace -fPIC -cpp")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  #set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  # set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
  set (CMAKE_Fortran_FLAGS_RELEASE "-xHost -O3 -no-prec-div -fPIC -cpp")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -traceback -check all \
    -check bounds -check uninit -debug all -check arg_temp_created \
    -fpe0 -ftrapuv -fPIC -cpp")
elseif (Fortran_COMPILER_NAME MATCHES "ftn.*")
  set (CMAKE_Fortran_FLAGS_RELEASE "-xHost -O3 -no-prec-div -fPIC -cpp")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fPIC -cpp")
elseif (Fortran_COMPILER_NAME MATCHES "g77")
  # g77
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

# Define some directories
set(SRC ${CMAKE_SOURCE_DIR}/src)
set(LIB ${CMAKE_SOURCE_DIR}/lib)
set(BIN ${CMAKE_SOURCE_DIR}/bin)

# Have the .mod files placed in the lib folder
set(CMAKE_Fortran_MODULE_DIRECTORY ${LIB})

# Add subdirectories
add_subdirectory(src/modules ${LIB})
add_subdirectory(src/programs ${BIN})
#add_subdirectory(src/distribution ${BIN})

# build a CPack driven installer package
include (InstallRequiredSystemLibraries)
# set (CPACK_RESOURCE_FILE_LICENSE  
#      "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
# set (CPACK_PACKAGE_VERSION_MAJOR "${Tutorial_VERSION_MAJOR}")
# set (CPACK_PACKAGE_VERSION_MINOR "${Tutorial_VERSION_MINOR}")
set (CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
include (CPack)

# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/distclean.cmake
)
