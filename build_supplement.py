#!/usr/bin/env python3
"""Compile the self-contained Supplementary Material."""
from __future__ import annotations
import shutil, subprocess
from pathlib import Path

def main() -> None:
    root = Path(__file__).resolve().parent
    pdflatex = shutil.which('pdflatex')
    if pdflatex is None:
        raise SystemExit('pdflatex was not found. Install a TeX distribution and add it to PATH.')
    name = 'supplementary_material_gamma'
    for current in range(1, 4):
        print(f'Compiling {name}: pass {current}/3', flush=True)
        subprocess.run([pdflatex, '-interaction=nonstopmode', '-halt-on-error', name + '.tex'], cwd=root, check=True)
    print('Supplementary material compiled successfully.')

if __name__ == '__main__':
    main()
