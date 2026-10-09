#!/usr/bin/env python

from pathlib import Path
from pandas import Timestamp
import sys
import os
from argparse import Namespace
from omegaconf import OmegaConf
from loguru import logger


def get_emfile_name(conf : dict, conf_checksum: str, cache_dir: Path) -> Path:
    logger.debug(f"*****\n{conf}\n*****")
    period = f'{Timestamp(conf.start):%Y%m%d}-{Timestamp(conf.end):%Y%m%d}'
    scenario = str(conf.name).replace(' ', '_')
    return Path(cache_dir) / f'{scenario}_{period}_{conf_checksum}.nc'


def gen_emfile(conf: dict, conf_checksum : str, cache_dir: Path) -> Path:
    from tm5.emissions import prepare_emissions
    from postrun.fitic_input4inversion import subcmd_monthly_emissions_for_inversion

    outname = get_emfile_name(conf, conf_checksum, cache_dir)
    msg = f"...determined outname ***{str(outname)}*** exists={outname.exists()} " \
        f"(cache_dir={str(cache_dir)})"
    logger.debug(msg)
    #
    #-- emissions file should be fully determined by checksum of it's configuration
    #
    if outname.exists():
        return outname
    else:
        emisfile_list = list(cache_dir.glob('*.nc'))
        msg = f"for configuration checksum ==>{conf_checksum}<== looking up equal scenario in " \
            f"*****{emisfile_list}*****"
        logger.debug(msg)
        
        for ifile,emisfile in enumerate(emisfile_list):
            msg = f"ifile={ifile}: ***{str(emisfile)}***"
            logger.debug(msg)
            msg = f"emisfile.stem ==>{emisfile.stem}<=="
            logger.debug(msg)
            tokens = emisfile.stem.split('_')
            msg = f"tokens ==>{tokens}<=="
            logger.debug(msg)
            #-- checksum is last token in filename
            cur_chksum = tokens[-1]
            msg = f"...found existing checksum ==>{cur_chksum}<=="
            logger.debug(msg)

            if cur_chksum==conf_checksum:
                msg = f"...detected equal emission configuration for existing file ***{str(emisfile)}***"
                logger.debug(msg)
                dir_fd = os.open(str(cache_dir), os.O_RDONLY)
                srcfile = emisfile.name
                dstfile = outname.name
                os.symlink(src=srcfile, dst=dstfile, dir_fd=dir_fd)
                return outname
            else:
                continue
    #
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
        # msg = f"@{catname} -->\n{cat}\n<--"
        # logger.debug(msg)
        cat_dict[catname] = {}
        cat_dict[catname]['path'] = cat['global']['file']
        cat_dict[catname]['field'] = cat['global']['field']
        if 'regional' in cat:
            cat_dict[catname]['overwrite'] = {}
            cat_dict[catname]['overwrite']['path'] = cat['regional']['file']
            cat_dict[catname]['overwrite']['field'] = cat['regional']['field']
    prefix = f'{scratch_dir}/ch4emis'
    emis_dict = {
        'run': {'start': str(conf.start), 'end': str(conf.end), 'regions': region_names},
        'regions': region_block(region_names),
        'emissions': {'CH4': {'prefix': prefix, 'categories': cat_dict } },
    }
#    msg = f"...calling prepare_emissions with emis_dict **********\n" \
#        f"{emis_dict}\n**********"
#    logger.debug(msg)
        
    emis_dconf = OmegaConf.create(emis_dict)
#    msg = f"...prepare_emissions being called with emis_dconf ******************************\n" \
#    f"{emis_dconf}\n******************************"
#    logger.debug(msg)
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
    # msg = f"...calling subcmd_monthly_emissions_for_inversion with args ==>{args}<=="
    # logger.debug(msg)
    subcmd_monthly_emissions_for_inversion(args)

    return outname
