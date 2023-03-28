#!/usr/bin/env python

from pathlib import Path
from typing import Union, List, Tuple
from loguru import logger
import re
from types import SimpleNamespace


def readline(fid) -> str:
    """
    Read a file line by line, handling the continued lines, as in: 
    key : x \
        y \
        z
    """
    def readnextline(line: str, fid):
        line = line.strip()
        if line.endswith('\\'):
            line = line.split('\\')[0] + ' ' + fid.readline().strip()
            return readnextline(line, fid)
        return line
    
    for l in fid:
        yield readnextline(l, fid)


class RcFile(dict):
    def __init__(self, filename: Union[str, Path]):
        self.root = Path(filename).parent.absolute()
        self.unresolved_lines = []
        self.unresolved_includes = []
        self.read_file(Path(filename).name)
        
    def read_file(self, filename):
        logger.info(f'Reading file {str(self.root / filename)}')
        
        includes = []
        with open(self.root / filename, 'r') as fid:
            lines = list(readline(fid))
            
            for line in lines:
                # logger.debug(line)
                if line.startswith('#include'):
                    includes.append(line.split()[1]) 
                else:
                    line = line.split('!')[0]
                    if line:
                        try :
                            k, v = line.split(':', maxsplit=1)
                            k = self.resolve_string(k)
                            self[k] = v
                            logger.warning(k)
                            if k == 'ndyn':
                                import pdb; pdb.set_trace()
                        except KeyError :
                            self.unresolved_lines.append(line)
        self.unresolved_includes.extend(includes)
        self.resolve()        

    def resolve(self) -> None:
        includes = list(self.unresolved_includes)
        unresolved = []
        for fname in includes:
            try :
                fname = self.resolve_string(fname)
                self.read_file(fname)
            except KeyError:
                self.unresolved.append(fname)
        self.unresolved_includes = unresolved                

    def resolve_string(self, txt: str) -> str:
        subst = re.findall("\$\{[\d\w\.]*\}", txt)
        for chars in subst:
            txt = txt.replace(chars, self.get(chars[2:-1]))
        return txt.strip()

    def get(self, key: str) -> str:
        return self.resolve_string(self[self.resolve_string(key)])

    def write(self, dest: Path, resolve : bool = True) -> Path:
        if resolve:
            self.resolve()
            
        with open(dest, 'w') as fid:
            for include in self.unresolved_includes:
                fid.write(f'#include {str(self.root / include)}')
            for k, v in sorted(self.items()):
                fid.write(f'{k:>30s} : {str(v)}\n')