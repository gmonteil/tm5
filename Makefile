
tm5:
	# Actual TM5 container
	singularity build images/tm5-intel-ubuntu20.simg images/tm5-intel-ubuntu20.def

envcontainer: images/intel-ubuntu.simg
	# install libraries needed by TM5: HDF4, HDF5, netCDF4, udunits, and their own dependencies
	# might need sudo ...

	singularity build images/tm5env-intel-ubuntu20.simg images/hdf4-netcdf4-intel-ubuntu20.def

intel: images/tm5env.sif
	# Construct a new clean container based on the old tm5env.sif 
	# might need sudo ...
	singularity build images/intel-ubuntu.simg images/intel.def

dev:
	singularity run \
		-B /data/projects/TM5:/input \
		-B /data:/data \
		-B /data/projects/TM5/meteo:/meteo \
		-B .:/tm5 \
		-B output:/output \
		images/tm5-intel-ubuntu20.simg
