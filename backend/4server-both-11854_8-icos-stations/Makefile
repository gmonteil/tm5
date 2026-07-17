#======================================================================
#
#   Name:               Makefile
#
#   Description:        Forward and Inversion Codes for FIT-IC
#                       The Inversion Lab Thomas Kaminski Consulting
#                       Thomas.Kaminski(at)Inversion-Lab.com
#
#   Date:            	April and May 2026
#
#======================================================================

PYPLOT="/home/thomas/sw/miniforge3/envs/fit-ic/bin/python /home/thomas/tm5/postrun/inversion_post.py"

FC = gfortran
FFLAGS = -Wall

ifeq ($(findstring gfortran,$(FC)),gfortran)
  FFOPTI = -O2
  ifdef DBG
    FFLAGS = -ggdb -Wall -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow
    FFLAGS += -pedantic
    FFLAGS += -finit-real=snan -finit-integer=-1000000
    FFOPTI = -O0
  endif
  ifdef PROFILE
    FFLAGS += -ggdb -pg
    FFOPTI = -O2
  endif
endif

FFLAGS += $(FFOPTI)

NETCDF_PREFIX = /usr
FFINC  = -I$(NETCDF_PREFIX)/include -I.
LDNETCDF = $(NETCDF_PREFIX)/lib64/libnetcdff.a -lnetcdf -ldl -lm
FFLAGS	+= $(FFOPTI) $(FFINC) 

#-- use MODE = obs1d to run forward/inversion with 1D observation vector
MODE = 
ifeq ($(MODE),obs1d)
  FSRCS = forward.f90 utilities_obs1D.f90
else
  FSRCS	= forward.f90 utilities.f90
endif
FOBJS    = ${FSRCS:.f90=.o}

ifeq ($(MODE),obs1d)
ISRCS	= eispack/eispack.f90 inversion.f90 utilities_obs1D.f90
else
  ISRCS	= eispack/eispack.f90 inversion.f90 utilities.f90
endif
IOBJS    = ${ISRCS:.f90=.o}

TIMEIT = /usr/bin/time -f "CPU time: %U seconds, Memory: %M Kbytes"

# forward run
frun: forward 
	cp foj.nc output/fc.nc
	$(TIMEIT) ./$^

forward: $(FOBJS) 
	$(FC) $(FFOPTI) -o $@ $^ $(LDNETCDF)

# inversion
run: inversion 
	cp foj.nc output/fcpost.nc
	cp foj.nc output/fepost.nc
	$(TIMEIT) ./$^

inversion: $(IOBJS) 
	$(FC) $(FFOPTI) -o $@ $^ $(LDNETCDF)

# emission plots
PFLAGZ = --extent 0 18 42 58 --outname zoom.png
PFLAG = --cbmax 0.14

plotpost:
	$(PYPLOT) emis_visu --tosqm $(PFLAG) $(PFLAGZ) --outdir post fepost.nc 

plotprior:
	$(PYPLOT) emis_visu --tosqm $(PFLAG) $(PFLAGZ) --outdir prior fe.nc 

plotpriorold:
	bash; source ~/sw/miniforge3/bin/activate; conda activate fit-ic; python /home/thomas/tm5/postrun/inversion_post.py emis_visu --tosqm $(PFLAG) $(PFLAGZ) --outdir prior fe.nc 

# concentration plot
SITE=CBW
ID=test
tsprpo: 
	(xmgrace -settype xy output/cpr.dat -pexec 's0 legend "prior"' \
	        -settype xy output/cpost.dat -pexec 's1 legend "post"' \
	        -settype xy output/cobs.dat -pexec 's2 legend "obs"' \
		-pexec 's2 line type 0' -pexec 's2 symbol 8' -pexec 's2 symbol size 0.5' \
		-pexec 'title "$(SITE) January 2021"' -pexec 'world 0, 0, 32, 30' -pexec 'legend 0.85, 0.9' \
		-pexec 'xaxis  tick major 2' -pexec 'xaxis  tick minor ticks 1' \
		-pexec 'xaxis  ticklabel stop type spec' -pexec 'xaxis  ticklabel stop 31.000000' \
		-pexec 'xaxis  ticklabel start type spec' -pexec 'xaxis  ticklabel start 1.000000' \
		-pexec 'yaxis  tick major 2' -pexec 'yaxis  tick minor ticks 1' \
		-pexec 'xaxis  label "day"' -pexec 'yaxis  label "ppb"' \
	        -hdevice PNG -printfile $(SITE)-$(ID)-prpo.png -saveall $(SITE)-$(ID)-prpo.agr -hardcopy ) 

# suffic rule
%.o: %.f90
	$(FC) $(FFLAGS) $(FFXTRA) -c -o $@ $<

# cleaning up
clean:
	rm -f *~ .*~ *.o *.mod $(IOBJS) inversion $(FOBJS) forward
