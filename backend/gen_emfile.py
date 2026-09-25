#!/usr/bin/env python

from pathlib import Path
from pandas import Timestamp
import sys
from argparse import Namespace
from omegaconf import OmegaConf
from loguru import logger
from tm5.emissions import prepare_emissions


def get_emfile_name(conf, cache_dir: Path) -> Path:
    period = f'{Timestamp(conf.start):%Y%m%d}-{Timestamp(conf.end):%Y%m%d}'
    scenario = str(conf.name).replace(' ', '_')
    return Path(cache_dir) / f'{scenario}_{period}.nc'


def gen_emfile(conf: dict, cache_dir: Path) -> Path:

    msg = f"...emissions request for conf ******************************\n" \
        f"{conf}\n******************************"
    logger.debug(msg)
    # Hack to make stuff in "postrun" importable ...
    repo_root = str(Path(__file__).resolve().parents[1])
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)
    from postrun.fitic_input4inversion import subcmd_monthly_emissions_for_inversion

    outname = get_emfile_name(conf, cache_dir)
    msg = f"...determined outname ***{str(outname)}*** exists={outname.exists()} " \
        f"(cache_dir={str(cache_dir)})"
    logger.debug(msg)
    
    if outname.exists():
        return outname

    region_names = list(conf.regions.keys())
    finest = [r for r in region_names if conf.regions[r]['dlon'] == 1 and conf.regions[r]['dlat'] == 1]
    coarser = [r for r in region_names if r not in finest]

    msg = f"coarser ==>{coarser}<==  finest ==>{finest}<=="
    logger.debug(msg)

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
        emis_dict = {
            'run': {'start': str(conf.start), 'end': str(conf.end), 'regions': coarser},
            'regions': region_block(coarser),
            'emissions': {'CH4': {'prefix': prefix, 'categories': categories_block('global')}},
        }

    if finest:
        emis_dict = {
            'run': {'start': str(conf.start), 'end': str(conf.end), 'regions': finest},
            'regions': region_block(finest),
            'emissions': {'CH4': {'prefix': prefix, 'categories': categories_block('regional')}},
        }

    msg = f"...calling prepare_emissions with emis_dict **********\n" \
        f"{emis_dict}\n**********"
    logger.debug(msg)
        
    emis_dconf = OmegaConf.create(emis_dict)
    prepare_emissions(emis_dconf)
    #
    #-- generate spatially 1D emissions for FIT-IC Fortran system.
    #
    args = Namespace(
        tm5emisdir=str(scratch_dir),
        time_range=[Timestamp(conf.start), Timestamp(conf.end)],
        regions=region_names,
        outdir=None,
        outname=str(outname),
    )
    subcmd_monthly_emissions_for_inversion(args)

    return outname


def gen_emfile_new(conf: dict, cache_dir: Path) -> Path:

    msg = f"...emissions request for conf ******************************\n" \
        f"{conf}\n******************************"
    logger.debug(msg)
    # Hack to make stuff in "postrun" importable ...
    repo_root = str(Path(__file__).resolve().parents[1])
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)
    from postrun.fitic_input4inversion import subcmd_monthly_emissions_for_inversion

    outname = get_emfile_name(conf, cache_dir)
    msg = f"...determined outname ***{str(outname)}*** exists={outname.exists()} " \
        f"(cache_dir={str(cache_dir)})"
    logger.debug(msg)
    
    if outname.exists():
        return outname

    region_names = list(conf.regions.keys())

    def region_block(names):
        return {
            name: {
                'lons': [conf.regions[name]['west'], conf.regions[name]['east'], conf.regions[name]['dlon']],
                'lats': [conf.regions[name]['south'], conf.regions[name]['north'], conf.regions[name]['dlat']],
            }
            for name in names
        }

    scratch_dir = Path(cache_dir) / 'scratch' / outname.stem
    scratch_dir.mkdir(parents=True, exist_ok=True)

    #
    #-- convert emissions dconf to make it suitable for prepare_emissions
    #
    cat_dict = {}
    for catname, cat in conf.categories.items():
        cat_dict[catname] = {}
        cat_dict[catname]['path'] = conf[catname]['global']['filename']
        cat_dict[catname]['feld'] = conf[catname]['global']['field']
        if 'regional' in conf[catname]:
            cat_dict[catname]['overwrite'] = {}
            cat_dict[catname]['overwrite']['path'] = conf[catname]['regional']['filename']
            cat_dict[catname]['overwrite']['field'] = conf[catname]['regional']['field']
    prefix = f'{scratch_dir}/ch4emis'
    emis_dict = {
        'run': {'start': str(conf.start), 'end': str(conf.end), 'regions': region_names},
        'regions': region_block(region_names),
        'emissions': {'CH4': {'prefix': prefix, 'categories': cat_dict } },
    }
    msg = f"...calling prepare_emissions with emis_dict **********\n" \
        f"{emis_dict}\n**********"
    logger.debug(msg)
        
    emis_dconf = OmegaConf.create(emis_dict)
    msg = f"...prepare_emissions being called with emis_dconf ******************************\n" \
    f"{emis_dconf}\n******************************"
    logger.debug(msg)
    prepare_emissions(emis_dconf)
    #
    #-- generate spatially 1D emissions for FIT-IC Fortran system.
    #
    args = Namespace(
        tm5emisdir=str(scratch_dir),
        time_range=[Timestamp(conf.start), Timestamp(conf.end)],
        regions=region_names,
        outdir=None,
        outname=str(outname),
    )
    subcmd_monthly_emissions_for_inversion(args)

    return outname
