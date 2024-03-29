#
# Makefile for Mike Fisher's conjugate gradient optimizer for TM5 4DVAR
#

# --- help --------------------------------------

help:
	@echo " "
	@echo " usage:"
	@echo " "
	@echo "   make congrad.exe # compile executable"
	@echo " "
	@echo "   make clean       # remove objects"
	@echo "   make clean-all   # remove objects and executable"
	@echo " "

# --- compiler and flags ------------------------

FC = ifort

FFLAGS_DEF = -m64 -fp-model strict -fomit-frame-pointer
#FFLAGS_DEF = -m64 -xSSE4.2 -fp-model strict -fomit-frame-pointer
FFLAGS_WP = -r8
FFLAGS_OPT = -O3 -align -parallel
#FFLAGS_PAR = -openmp
FFLAGS_CHECK = -O2 -check all -g -traceback -fpe0 -ftrapuv

NETCDF_HOME = $(NETCDF4_DIR)

FFLAGS = $(FFLAGS_DEF) $(FFLAGS_WP) $(FFLAGS_OPT)
#FFLAGS = $(FFLAGS_DEF) $(FFLAGS_WP) $(FFLAGS_CHECK)

NETCDF_INCS = $(shell nf-config --fflags)
NETCDF_LIBS = $(shell nf-config --flibs)

# --- libs etc ---------------------------------------

INCS = $(NETCDF_INCS)
LIBS = $(NETCDF_LIBS)

# --- implicit rules ---------------------------------

# how to form object files from Fortran source:
%.o: %.F90
	$(FC) $< -c -o $@ $(FFLAGS) $(INCS)
	@echo ' '

lapack.o         	: lapack.F90
m_option_parser.o	: m_option_parser.F90
conjuGrad.o      	: conjuGrad.F90 lapack.o m_option_parser.o
conjuGrad_main.o 	: conjuGrad_main.F90 conjuGrad.o

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# targets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

congrad.exe: conjuGrad_main.o
	$(FC) -o $@ $(FFLAGS) *.o $(LIBS)

clean:
	rm -f *.{mod,o}

clean-all:
	rm -f *.{mod,o} congrad.exe
