#!/usr/bin/env python3
"""Run exact, symbolic and multiprecision supplementary checks.

  python code/run_verifications.py --rational-only
  python code/run_verifications.py
  python code/run_verifications.py --experiments

The last command also recomputes eleven Nyström experiments at 1000 digits.
Outputs are written below results/. Exact checks and numerical checks are
identified separately in results/verification/summary.json.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path


def main() -> None:
    if not __debug__:
        raise SystemExit('Do not use -O: the mathematical checks require Python assertions.')
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--rational-only', action='store_true', help='Run only checks using fractions.Fraction.')
    group.add_argument('--experiments', action='store_true', help='Also reproduce all eleven S5 experiments and their summaries.')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    logs = root/'results'/'verification'
    logs.mkdir(parents=True, exist_ok=True)
    jobs = [
        ('criterion', 'exact rational', ['verify_gamma_criterion.py', '--out', 'results']),
        ('comparison_2015', 'exact rational', ['compare_cordero2015.py', '--out', 'results']),
        ('obstruction', 'exact rational', ['verify_exact_obstruction.py']),
        ('lotfi_orders', 'exact rational series', ['verify_lotfi_orders.py']),
    ]
    if not args.rational_only:
        jobs.extend([
            ('algebra', 'symbolic exact', ['verify_algebra.py']),
            ('symbolic_majorant', 'symbolic exact', ['symbolic_majorant_verification.py']),
            ('published_results', 'exact checks and labelled multiprecision reproductions',
             ['verify_published_results.py', '--dps', '2000', '--out', 'results/literature_checks']),
            ('maroju2022', 'symbolic exact', ['verify_maroju2022.py']),
            ('maroju_bounds', 'exact bounds and labelled multiprecision roots',
             ['verify_maroju_bounds.py', '--out', 'results/maroju_bounds.json']),
            ('lotfi_numerical', 'multiprecision consistency', ['check_lotfi_numerically.py']),
        ])
    if args.experiments:
        for name, nodes in [('all', ['30']), ('chandrasekhar', ['20', '40', '60']), ('bratu', ['20', '40', '60'])]:
            jobs.append(('experiments_'+name, 'multiprecision floating point',
                         ['reproduce_experiments.py', '--problem', name, '--N', *nodes, '--out', 'results/experiments']))
        jobs.append(('experiment_summary', 'multiprecision profiles',
                     ['reproduce_experiments.py', '--summarize', '--out', 'results/experiments']))
    report = []
    for label, arithmetic, task in jobs:
        script = root/'code'/task[0]
        command = [sys.executable, '-B', str(Path('code')/task[0]), *task[1:]]
        print('Running:', label, flush=True)
        start = time.perf_counter()
        with (logs/(label+'.log')).open('w', encoding='utf-8') as log:
            completed = subprocess.run(command, cwd=root, stdout=log, stderr=subprocess.STDOUT, check=False)
        report.append({'check': label, 'arithmetic': arithmetic, 'command': ['python', *command[1:]],
                       'script_sha256': hashlib.sha256(script.read_bytes()).hexdigest(),
                       'returncode': completed.returncode, 'elapsed_seconds': time.perf_counter()-start})
        (logs/'summary.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
        if completed.returncode:
            raise SystemExit(f'Check failed: {label}. See results/verification/{label}.log.')
    if args.experiments:
        subprocess.run([sys.executable, '-B', 'code/verify_experiment_data.py'], cwd=root, check=True)
    print('All requested checks completed successfully.', flush=True)


if __name__ == '__main__':
    main()
