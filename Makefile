envcontainer: intel-ubuntu.simg
	# install libraries needed by TM5: HDF4, HDF5, netCDF4, udunits, and their own dependencies
	sudo singularity build images/tm5env-intel-ubuntu20.simg images/hdf4-netcdf4-intel-ubuntu20.def

intel: images/tm5env.sif
	# Construct a new clean container based on the old tm5env.sif 
	singularity build images/intel-ubuntu.simg images/intel.def