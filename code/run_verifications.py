#!/usr/bin/env python3
"""Runs the supplementary verifications and, optionally, the experiments.

From the package directory:
  python code/run_verifications.py --rational-only
  python code/run_verifications.py
  python code/run_verifications.py --experiments

The last command reruns eleven experiments at 1000 digits.
The supplied results already include these eleven runs.
"""
from __future__ import annotations
import argparse, json, subprocess, sys, time
from pathlib import Path

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--rational-only', action='store_true', help='Only the criterion and the rational comparison with 2015')
    ap.add_argument('--experiments', action='store_true', help='Also repeat the eleven experiments in S5')
    args = ap.parse_args()
    if args.rational_only and args.experiments:
        ap.error('Do not combine --rational-only with --experiments')
    code = Path(__file__).resolve().parent
    root = code.parent
    results = root / 'results'
    logs = results / 'combined_run'
    logs.mkdir(parents=True, exist_ok=True)
    jobs = [
        ('criterion', ['verify_gamma_criterion.py', '--out', str(results)]),
        ('comparison_2015', ['compare_cordero2015.py', '--out', str(results)]),
    ]
    if not args.rational_only:
        jobs.extend([
            ('obstruction', ['verify_exact_obstruction.py']),
            ('algebra', ['verify_algebra.py']),
            ('symbolic_majorant', ['symbolic_majorant_verification.py']),
            ('audit', ['verify_audit.py', '--dps', '2000', '--out', str(results / 'audit')]),
            ('tjm2022', ['independent_verification_tjm2022.py']),
            ('additional_2022', ['check_additional_2022.py', '--out', str(results / 'additional_2022.json')]),
        ])
    if args.experiments:
        out = str(results / 'experiments')
        for name, ns in [('all', ['30']), ('chandrasekhar', ['20', '40', '60']), ('bratu', ['20', '40', '60'])]:
            jobs.append(('experiments_' + name, ['reproduce_experiments.py', '--problem', name, '--N', *ns, '--out', out]))
        jobs.append(('experiment_summary', ['reproduce_experiments.py', '--summarize', '--out', out]))
    report = []
    for label, task in jobs:
        command = [sys.executable, str(code / task[0]), *task[1:]]
        start = time.perf_counter()
        print('Running:', label, flush=True)
        with (logs / (label + '.log')).open('w', encoding='utf-8') as log:
            result = subprocess.run(command, cwd=root, stdout=log, stderr=subprocess.STDOUT, check=False)
        report.append({'test': label, 'returncode': result.returncode, 'seconds': time.perf_counter() - start})
        (logs / 'summary.json').write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
        if result.returncode:
            raise SystemExit('Failure in ' + label + '. See ' + str(logs / (label + '.log')))
    print('All requested runs completed successfully.', flush=True)

if __name__ == '__main__':
    main()
