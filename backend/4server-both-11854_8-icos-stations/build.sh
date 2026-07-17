
gfortran -o forward utilities_obs1D.f90 forward.f90 -I/usr/lib64/gfortran/modules -I/usr/include -lnetcdff
gfortran -o inversion eispack/eispack.f90 inversion.f90 utilities_obs1D.f90 -I/usr/lib64/gfortran/modules -I/usr/include -lnetcdff
