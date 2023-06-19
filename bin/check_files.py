import sys
from pathlib import Path

for file in Path(sys.argv[1]).glob(f'ec/ea/*/*/{sys.argv[2]}/*/*/*.status'):
    with open(file, 'r') as fid:
        status = fid.readline().strip()

    if status == 'in-progress':
        print(file.parent / Path(file.stem))
        (file.parent / Path(file.stem)).unlink()
    file.unlink()
