# Standard TM5 workflows

To facilitate the implementation and testing of FIT-IC, some test cases were defined, which should be run once in a while to ensure that the code runs as expected

## Meteo coarsening

```bash
export TM5_HOST=laptop
conda activate tm5
python coarsen_meteo.py coarsen_meteo.yaml
```


## Forward CO$_2$

```bash
export TM5_HOST=laptop
conda activate tm5
python forward.py forward.yaml
```