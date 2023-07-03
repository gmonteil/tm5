# runtm5

## Development notes

### Region definitions

There is a system of templates in pycasso that generates (at least) the following files:

- dims_grid.F90

Since I use only a global zoom, I just replaced that whole mess by the generated fortran file, in *src/proj/regions/glb6x4*. The procedure to create another region definition would be to create another subproject in *src/proj/regions* and include it instead.