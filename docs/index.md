# runtm5

## Project layout

The wrapper is composed of several layers: 
1. The actual **TM5** source code (in fortran)
2. The **pycasso** scripts (in python2.7), which control the compilation and execution of TM5
3. The **pyshell** python package (python2.7) which performs the inversions as well as lots of pre-postprocessing tasks, and controls pycasso.
4. The **runtm5** (python3.9+) python package controls pyshell

TM5/pycasso/pyshell are ran inside a singularity container (with image and definition files located under the *images*) folder.

Functions from pycasso/pyshell are progressively being implemented in runflex, the aim is that it should be able to control TM5 directly.

## Running TM5 simulations with runtm5

### Forward CO2 simulation:
```
tm5 --dev pyshell forward --rc forward.yaml
```
(the call to *pyshell* will disappear once it is fully deprecated)

### Background extraction

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