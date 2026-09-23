#!/usr/bin/env python

from pathlib import Path
from pandas import Timestamp
import sys
from argparse import Namespace
from omegaconf import OmegaConf
from tm5.emissions import prepare_emissions


def get_emfile_name(conf, cache_dir: Path) -> Path:
    period = f'{Timestamp(conf.start):%Y%m%d}-{Timestamp(conf.end):%Y%m%d}'
    scenario = str(conf.name).replace(' ', '_')
    return Path(cache_dir) / f'{scenario}_{period}.nc'


def gen_emfile(conf: dict, cache_dir: Path) -> Path:

    # Hack to make stuff in "postrun" importable ...
    repo_root = str(Path(__file__).resolve().parents[1])
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)
    from postrun.fitic_input4inversion import subcmd_monthly_emissions_for_inversion

    outname = get_emfile_name(conf, cache_dir)
    if outname.exists():
        return outname

    region_names = list(conf.regions.keys())
    finest = [r for r in region_names if conf.regions[r]['dlon'] == 1 and conf.regions[r]['dlat'] == 1]
    coarser = [r for r in region_names if r not in finest]

    def region_block(names):
        return {
            name: {
                'lons': [conf.regions[name]['west'], conf.regions[name]['east'], conf.regions[name]['dlon']],
                'lats': [conf.regions[name]['south'], conf.regions[name]['north'], conf.regions[name]['dlat']],
            }
            for name in names
        }

    def categories_block(scope):
        return {
            catname: {'path': (cat[scope] if scope in cat else cat['global']).file,
                      'field': (cat[scope] if scope in cat else cat['global']).field}
            for catname, cat in conf.categories.items()
        }

    scratch_dir = Path(cache_dir) / 'scratch' / outname.stem
    scratch_dir.mkdir(parents=True, exist_ok=True)

    prefix = f'{scratch_dir}/ch4emis'

    if coarser:
        prepare_emissions(OmegaConf.create({
            'run': {'start': str(conf.start), 'end': str(conf.end), 'regions': coarser},
            'regions': region_block(coarser),
            'emissions': {'CH4': {'prefix': prefix, 'categories': categories_block('global')}},
        }))

    if finest:
        prepare_emissions(OmegaConf.create({
            'run': {'start': str(conf.start), 'end': str(conf.end), 'regions': finest},
            'regions': region_block(finest),
            'emissions': {'CH4': {'prefix': prefix, 'categories': categories_block('regional')}},
        }))

    args = Namespace(
        tm5emisdir=str(scratch_dir),
        time_range=[Timestamp(conf.start), Timestamp(conf.end)],
        regions=region_names,
        outdir=None,
        outname=str(outname),
    )
    subcmd_monthly_emissions_for_inversion(args)

    return outname
