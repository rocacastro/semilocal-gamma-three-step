#!/usr/bin/env python3
"""Check or generate SHA-256 digests of the distributed package files.

Use --write only to generate a manifest for a deliberately prepared package.
Normal verification never changes files. Regenerated execution metadata and
PDFs can have different byte-level digests; see REPRODUCIBILITY.md.
"""
from __future__ import annotations
import argparse
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = 'SHA256SUMS.txt'
SKIP_DIRECTORIES = {'.git', '.venv', 'venv', '__pycache__', '.pytest_cache', '.mypy_cache', '.idea', '.vscode'}
SKIP_SUFFIXES = {'.pyc', '.pyo', '.aux', '.toc', '.out', '.fls', '.fdb_latexmk', '.synctex.gz'}


def files() -> list[Path]:
    return sorted(p for p in ROOT.rglob('*') if p.is_file()
                  and p.name not in {MANIFEST, '.DS_Store', 'Thumbs.db', 'supplementary_material_gamma.log'}
                  and not any(part in SKIP_DIRECTORIES for part in p.relative_to(ROOT).parts)
                  and p.suffix not in SKIP_SUFFIXES)


def digest(path: Path) -> str:
    checksum = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            checksum.update(block)
    return checksum.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true', help='Generate a manifest from the present package files.')
    args = parser.parse_args()
    manifest = ROOT/MANIFEST
    if args.write:
        lines = [f'{digest(p)}  {p.relative_to(ROOT).as_posix()}' for p in files()]
        manifest.write_text('\n'.join(lines)+'\n', encoding='utf-8')
        print(f'Wrote {len(lines)} SHA-256 entries.')
        return
    if not manifest.is_file():
        raise SystemExit(f'Missing {MANIFEST}.')
    listed = {}
    for number, line in enumerate(manifest.read_text().splitlines(), 1):
        checksum, separator, relative = line.partition('  ')
        if not separator or len(checksum) != 64 or any(c not in '0123456789abcdef' for c in checksum):
            raise SystemExit(f'Invalid manifest entry at line {number}.')
        candidate = ROOT/relative
        if not candidate.resolve().is_relative_to(ROOT.resolve()) or relative in listed:
            raise SystemExit(f'Unsafe or duplicate path in the manifest: {relative}.')
        listed[relative] = checksum
    actual = {p.relative_to(ROOT).as_posix(): p for p in files()}
    failures = [f'Missing: {name}' for name in listed.keys()-actual.keys()]
    failures += [f'Unlisted: {name}' for name in actual.keys()-listed.keys()]
    failures += [f'Changed: {name}' for name in listed.keys() & actual.keys() if digest(actual[name]) != listed[name]]
    if failures:
        raise SystemExit('\n'.join(failures))
    print(f'PASS: all {len(listed)} file digests match.')


if __name__ == '__main__':
    main()
