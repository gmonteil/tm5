#
# Makefile for TM5
#

# ------------------------------------------------
# compiler and flags
# ------------------------------------------------

include Makefile_flags


# ------------------------------------------------
# implicit rules
# ------------------------------------------------

# how to form object files from F90 source:
%.o: %.F90
	$(FC) -c -o $@ $< $(FFLAGS)
	@echo ' '

# how to form object files from f90 source:
%.o: %.f90
	$(FC) -c -o $@ $< $(FFLAGS)
	@echo ' '


# ------------------------------------------------
# special rules
# ------------------------------------------------

# compilation of mdf fails on some machine if optim flags are enabled ..
mdf.o: mdf.F90
	$(FC) -c -o $@ $< $(FFLAGS_BASIC)
	@echo ' '


# ------------------------------------------------
# dependencies
# ------------------------------------------------

include Makefile_deps


# ------------------------------------------------
# end
# ------------------------------------------------

