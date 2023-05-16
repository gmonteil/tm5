#!/usr/bin/env python
import tempfile

from omegaconf import DictConfig
from pandas import Timestamp, date_range, Timedelta
from typing import Set, List
from pathlib import Path
import tempfile
from loguru import logger
from tqdm import tqdm
from subprocess import Popen, PIPE
from dataclasses import dataclass
from pandas.tseries.frequencies import to_offset


@dataclass(kw_only=True)
class Meteo:
    regions : List[str]
    levels : str
    path : Path | str
    archive : str

    # Unused (so far) but valid arguments :
    coarsened : bool = False
    output : bool = False
    output_path : str = None

    def gen_file_list(self, prefix: str, suffix: str, fields: List[str], start: Timestamp, end: Timestamp) -> Set[str]:
        files = []
        for day in date_range(start, end, inclusive='left'):
            for field in fields :
                files.append(day.strftime(prefix) + field + day.strftime(suffix))
        return set(files)

    def get(self, prefix: str, suffix: str, fields: List[str], start: Timestamp, end: Timestamp) -> bool:
        files = self.gen_file_list(prefix, suffix, fields, start, end)
        success = self.check_retrieve(files, msg = prefix + "{fields}" + suffix)
        if not success :
            [logger.error(f'Failed to download {f}') for f in files if not Path(self.path).joinpath(f).exists()]
        return success

    def check_retrieve(self, files: Set[str], attempts : int = 0, max_attempts : int = 3, msg : str = None) -> bool:
        while any(not Path(self.path).joinpath(f).exists() for f in files) and attempts < max_attempts:
            if attempts > 0 :
                missing_files = [f for f in files if not Path(self.path).joinpath(f).exists()]
            self.retrieve(files, msg)
            attempts += 1
        return True

    def retrieve(self, files: Set[str], msg: str = None) -> None:
        with tempfile.NamedTemporaryFile(delete=False) as fid :
            fid.writelines((file + '\n').encode() for file in files)

        cmd = f'rclone copy -P {self.archive} --include-from {fid.name} {self.path}'
        logger.info(cmd)

        # adapted from https://gist.github.com/wholtz/14b3a3479fe9c70eefb2418f091d2103
        with tqdm(total=100, unit='%', desc=msg) as pbar :
            with Popen(cmd.split(), stdout=PIPE, bufsize=1, universal_newlines=True) as proc:
                for line in proc.stdout:
                    line = line.strip()
                    if line.startswith('Transferred:') and line.endswith('%'):
                        percent = float(line.split(',')[1].split('%')[0])
                        pbar.n = percent
                        pbar.refresh()

        # Remove the temporary file, just for cleanliness:
        Path(fid.name).unlink(missing_ok=True)

    def setup_files(self, start: Timestamp, end: Timestamp) -> None:
        """
        Download meteo data
        """
        if self.coarsened:
            regions = self.regions
            levels = self.levels
        else :
            regions = ['glb100x100']
            levels = 'ml137'
        for region in regions:
            for tt in date_range(Timestamp(start).strftime('%Y%m01'), Timestamp(end) + Timedelta(days=1), freq='MS', inclusive='left'):
                self.get(prefix=f'ec/ea/h06h18tr3/{levels}/{region}/{tt:%Y/%m}/', suffix='_%Y%m%d_00p03.nc', fields=['convec', 'mfuv', 'mfw', 'q', 'sp', 't', 'tsp'], start = tt, end = tt + to_offset('MS'))
                self.get(prefix=f'ec/ea/h06h18tr1/sfc/{region}/{tt:%Y/%m}/', suffix='_%Y%m%d_00p01.nc', fields=['blh', 'ewss', 'nsss', 'sshf', 'slhf', 'u10m', 'v10m'], start=tt, end=tt + to_offset('MS'))
                self.get(prefix=f'ec/ea/an0tr1/sfc/{region}/{tt:%Y/%m}/', start = tt, end = tt + to_offset('MS'), suffix='_%Y%m%d_00p01.nc', fields=['sr']) #, 'albedo', 'veg'])
                self.get(prefix=f'ec/ea/an0tr1/sfc/{region}/{tt:%Y/%m}/', start = tt, end = tt + to_offset('MS'), suffix='_%Y%m.nc', fields = ['srols'])
                self.get(prefix=f'ec/ea/an0tr1/sfc/glb100x100/', suffix='.nc', fields=['oro', 'lsm'], start = tt, end = tt + to_offset('MS'))

