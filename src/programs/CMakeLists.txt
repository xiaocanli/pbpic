#include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# # For a large number of source files you can create it in a simpler way
# # using file() function:
# file(GLOB analysis_SOURCES *.f90)

# Create an executable files.
add_executable(pbpic.exec pic2d.f90)
target_link_libraries(pbpic.exec ${MPI_Fortran_LIBRARIES} ${HDF5_Fortran_LIBRARIES}
    fftw3 fftw3_mpi)

set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
install (TARGETS pbpic.exec DESTINATION ${CMAKE_INSTALL_PREFIX})
