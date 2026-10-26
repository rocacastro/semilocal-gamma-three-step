#!/usr/bin/env python3
"""Check the eleven S5 records, hashes, increments, ACOC and table values.

These are consistency and regression checks on multiprecision data.
They are not interval validation or a bound on Nyström discretization error.
"""
from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    directory = ROOT/'results'/'experiments'
    references = json.loads((ROOT/'tests'/'reference_values.json').read_text())
    script_hash = hashlib.sha256((ROOT/'code'/'reproduce_experiments.py').read_bytes()).hexdigest()
    expected = {(name, 30) for name in references['table3']}
    expected |= {(name, n) for name in ('chandrasekhar', 'bratu') for n in (20, 40, 60)}
    records = {}
    status = []
    for path in sorted(directory.glob('*/result.json')):
        data = json.loads(path.read_text())
        key = data['problem'], data['N']
        require(key not in records, f'Duplicate experiment: {key}.')
        records[key] = data
        require(key in expected, f'Unexpected experiment: {key}.')
        require(data['precision_digits'] == 1000 and data['cycles'] == 5,
                f'Protocol mismatch: {key}.')
        require(data['script_sha256'] == script_hash, f'Source hash mismatch: {key}.')
        require(len(data['history']) == 6, f'History length mismatch: {key}.')
        vectors = json.loads(path.with_name('complete_iterates.json').read_text())
        require([v['n'] for v in vectors] == list(range(6)), f'Incomplete iterates: {key}.')
        require(all(len(v['x']) == key[1] for v in vectors), f'Wrong vector dimension: {key}.')
        require(all(len(v.get('previous_y', [])) == key[1] and
                    len(v.get('previous_z', [])) == key[1] for v in vectors[1:]),
                f'Incomplete intermediate stages: {key}.')
        with mp.workdps(1100):
            x = [[mp.mpf(v) for v in row['x']] for row in vectors]
            increments = [max(abs(a-b) for a,b in zip(x[n], x[n-1])) for n in range(1,6)]
            for n, increment in enumerate(increments, 1):
                stored = mp.mpf(data['history'][n]['increment_inf'])
                require(abs(increment/stored-1) < mp.mpf('6e-24'), f'Increment mismatch: {key}, n={n}.')
            for n in (3,4,5):
                acoc = mp.log(increments[n-1]/increments[n-2])/mp.log(increments[n-2]/increments[n-3])
                stored = mp.mpf(data['history'][n]['acoc'])
                require(abs(acoc-stored) < mp.mpf('6e-15'), f'ACOC mismatch: {key}, n={n}.')
                require(data['history'][n]['acoc_status'] == 'usable', f'Unusable reported ACOC: {key}.')
            require(all(mp.mpf(v) < mp.mpf('1e-970') for v in data['quadrature_moment_errors']),
                    f'Quadrature moment residual mismatch: {key}.')
            require(all(mp.mpf(v) < mp.mpf('1e-970') for v in data['linear_equation_residuals']),
                    f'Linear-solve residual mismatch: {key}.')
            if key[1] == 30:
                ref = references['table3'][key[0]]
                residual = mp.mpf(data['history'][4]['residual_inf'])
                require(abs(residual/mp.mpf(ref['residual_4'])-1) < mp.mpf('5.1e-7'),
                        f'Table 3 residual mismatch: {key}.')
                require(abs(mp.mpf(data['history'][5]['acoc'])-mp.mpf(ref['acoc_5'])) <= mp.mpf('5e-5'),
                        f'Table 3 ACOC mismatch: {key}.')
        status.append({'problem': key[0], 'N': key[1], 'protocol_checked': True,
                       'source_hash_matches': True, 'increments_and_acoc_recomputed': True})
    require(set(records) == expected, 'The eleven required experiments are not all present.')
    with (directory/'refinement.csv').open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    require(len(rows) == 8, 'Eight profile comparisons are required.')
    with mp.workdps(100):
        for row in rows:
            name, n = row['problem'], row['N']
            if n in references['table4'][name]:
                target = mp.mpf(references['table4'][name][n])
                field = 'difference_vs_exact_Bratu' if name == 'bratu' else 'difference_vs_N60'
                value = mp.mpf(row[field])
                # Half a last printed decimal place for a three-significant-digit value.
                tolerance = mp.mpf('0.5001') * mp.power(10, mp.floor(mp.log10(target))-2)
                require(abs(value-target) <= tolerance, f'Table 4 rounding mismatch: {name}/{n}.')
        ref = json.loads((directory/'Bratu_reference.json').read_text())
        a = mp.mpf(ref['a'])
        require(0 < a < 2 and abs(2*a*a/mp.cosh(a/2)**2-2) < mp.mpf('1e-95'),
                'Bratu reference does not satisfy the stated equation.')
    with (directory/'profiles_grid401.csv').open(newline='') as stream:
        profiles = list(csv.DictReader(stream))
    require(len(profiles) == 8*401, 'Incomplete 401-point profile data.')
    output = {'status': 'PASS', 'experiment_count': len(records), 'records': status,
              'table3_and_table4_roundings_match': True, 'profile_row_count': len(profiles),
              'reproduction_script_sha256': script_hash,
              'scope': 'Consistency of numerical records, not a discretization-error certificate.'}
    (directory/'RUN_STATUS.json').write_text(json.dumps(output, indent=2)+'\n', encoding='utf-8')
    print('PASS: 11 experiments, 55 increments, 33 ACOC values, table roundings, source hashes and profile data.')


if __name__ == '__main__':
    main()
