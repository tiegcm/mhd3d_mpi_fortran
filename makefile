#  WARNING: THIS MAKEFILE ONLY WORKS WITH GNU MAKE

########## PROGRAM name is program root name

PROGNAME      = ./mhd3d2

##########  Find out Machine Type ####################

ifdef debug
   DEBUG_OPT = -g
   no_opt    =  1
else
   DEBUG_OPT = 
endif

# HDF5 paths to include directory, libhdf5_fortran.a, libhdf5.a and other external hdf5 libraries used
# HDF5 = /opt/hdf5/1.8.12-intel13.0
#
F90FLAGS    =  -fp-model precise -O3 -r8\
               -Wl,--start \
               -L$(HDF5_DIR)/lib -lhdf5 -lhdf5_hl -I$(HDF5_DIR)/include -lhdf5_fortran -lhdf5hl_fortran \
                 $(HDF5_DIR)/lib/libhdf5.a $(HDF5_DIR)/lib/libhdf5_hl.a \
                 $(HDF5_DIR)/lib/libhdf5_fortran.a $(HDF5_DIR)/lib/libhdf5hl_fortran.a -lz \
               -Wl,--end 

CF          = mpif90
CPP         = cpp

FLAGS       = $(F90FLAGS) $(F90_OPT) $(DEBUG_OPT)

$(PROGNAME).x: $(PROGNAME).f90
	$(CF) -o $(PROGNAME).x $(PROGNAME).f90 $(FLAGS)
	rm -f *.o *.mod
$(PROGNAME): 
	echo 'Error: You need the .x extension'
