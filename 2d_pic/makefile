#
## Location of HDF5 binaries (with include/ and lib/ underneath)
HDF5 = /home/xli/hdf5

## Location of fftw3 binaries (with include/ and lib/ underneath)
FFTW3 = /home/xli/local/fftw-3.3.4

# Location of External Libraries
LIBZ    = /home/xli/local/zlib-1.2.8/lib/libz.a
LIBSZ   = /home/xli/local/szip-2.1/lib/libsz.a
#
# # Compiler
FC = mpif90
#
# ------ No machine-specific paths/variables after this  -----
#
FORTRANLIB=-I$(HDF5)/include $(HDF5)/lib/libhdf5_fortran.a
#
FSOURCE = pic2d
OBJECTS = pic2d.o
#
FLAGS = -c
FORTRANLIBS=
LIBSHDF   = $(FORTRANLIB) $(HDF5)/lib/libhdf5.a
FFTWHEADER = -I$(FFTW3)/include
FFTWLIB = -L$(FFTW3)/lib
LIB = $(LIBZ) $(LIBSZ) $(FFTWHEADER)  $(FFTWLIB) -lfftw3_mpi -lfftw3 -ldl
#
all: $(FSOURCE)
#
pic2d: pic2d.f90 
	$(FC) -o $@ pic2d.f90 $(LIBSHDF) $(LIB)

clean:
	rm -f $(FSOURCE) $(OBJECTS) *.h5

.f90.o: $(FSOURCE)
		$(FC) $(FLAGS) $?
.SUFFIXES:.o.c.f90
