# runtm5

## Project layout

The wrapper is composed of several layers: 

1. The actual **TM5** source code (in fortran)
2. The **pycasso** scripts (in python2.7), which control the compilation and execution of TM5
3. The **pyshell** python package (python2.7) which performs the inversions as well as lots of pre-postprocessing tasks, and controls pycasso.
4. The **runtm5** (python3.9+) python package controls pyshell
5. the **tm5** (python3.9+) python package is meant to control TM5 without using the pyshell/pycasso at all

TM5/pycasso/pyshell are ran inside a singularity container (with image and definition files located under the *images*) folder.

Functions from pycasso/pyshell are progressively being implemented in runflex, the aim is that it should be able to control TM5 directly.

## Building TM5

Building TM5 can be done using the `tmpy.py` script from the `tm5` python package:

```
tmpy.py build config.yaml
```

## Running TM5 simulations with runtm5

### Forward CO2 simulation:
```
tm5 --dev pyshell forward --rc forward.yaml
```
(the call to *pyshell* will disappear once it is fully deprecated)

### Background extraction

#### with tmflex

```
tm5 --dev pyshell forward --rc background.yaml
```

The difference with the simple forward run is the presence of a `background` section in the `output` section of the yaml file:
```
output :
    ...
    background :
        tracers : ${run.tracers}
        lat_range : [30, 74]
        lon_range : [-18, 36]
```

#### Native TM5 implementation (mask.apply)

### Inversions

```
tm5 --dev pyshell optim --rc optim.yaml
```
* the key `optim.emissions` must be set to 1

### Coarsening of meteorology

## Projects 

### Pyshell-bypassing

#### Compilation

#### Forward run

## Development notes

### Region definitions

There is a system of templates in pycasso that generates (at least) the following files:

- dims_grid.F90

Since I use only a global zoom, I just replaced that whole mess by the generated fortran file, in *src/proj/regions/glb6x4*. The procedure to create another region definition would be to create another subproject in *src/proj/regions* and include it instead.