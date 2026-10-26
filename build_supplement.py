#!/usr/bin/env python3
"""Compile and check the standalone mathematical supplement."""
from __future__ import annotations
import hashlib
import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path


def main() -> None:
    root = Path(__file__).resolve().parent
    executable = shutil.which('pdflatex')
    if executable is None:
        raise SystemExit('pdflatex was not found. Install a TeX distribution and add it to PATH.')
    source = root/'supplementary_material_gamma.tex'
    stem = source.stem
    with tempfile.TemporaryDirectory(prefix='gamma_supplement_') as directory:
        work = Path(directory)
        shutil.copy2(source, work/source.name)
        for iteration in range(1, 4):
            print(f'Compiling the supplement: pass {iteration}/3', flush=True)
            run = subprocess.run([executable, '-interaction=nonstopmode', '-halt-on-error', source.name],
                                 cwd=work, capture_output=True, text=True, errors='replace')
            if run.returncode:
                raise SystemExit('LaTeX compilation failed:\n'+run.stdout[-6000:]+run.stderr[-1000:])
        log = (work/(stem+'.log')).read_text(errors='replace')
        problems = [name for name, pattern in [
            ('undefined references or citations', r'(?:Reference|Citation)[\s\S]{0,180}?undefined|There were undefined references'),
            ('duplicate labels', r'multiply defined'),
            ('overfull boxes', r'Overfull \\[hv]box'),
            ('unresolved rerun request', r'Rerun to get cross-references right'),
        ] if re.search(pattern, log)]
        if problems:
            raise SystemExit('LaTeX validation failed: '+', '.join(problems))
        match = re.search(r'Output written on .*?\((\d+) pages?', log, re.S)
        if not match:
            raise SystemExit('Could not identify the generated PDF in the LaTeX log.')
        pages = int(match.group(1))
        shutil.copy2(work/(stem+'.pdf'), root/(stem+'.pdf'))
    report = {'status': 'PASS', 'pages': pages, 'passes': 3,
              'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
              'undefined_references_or_citations': 0, 'duplicate_labels': 0, 'overfull_boxes': 0,
              'external_source_dependencies': []}
    destination = root/'results'/'verification'/'supplement_build.json'
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(f'Supplement compiled and checked: {pages} pages.')


if __name__ == '__main__':
    main()
