#!/usr/bin/env python3
# Copyright (c) 2026 The minisketch developers
# Distributed under the MIT software license, see the accompanying
# file LICENSE or http://www.opensource.org/licenses/mit-license.php.

"""Check doc/spec/claims.md against the Lean project and the Python oracle.

For every claim:
- the `status` is one of the known values;
- each `oracle:` name (first dotted component) is defined in tests/pyminisketch.py,
  unless marked `(planned)`;
- a `proved` or `stated` claim has a `lean:` theorem that exists;
- a `proved` claim's theorem uses no sorryAx and only the standard axioms
  (plus Lean.ofReduceBool when the claim says `native_decide: allowed`).

Lean checks need a built project (`lake build`); run from anywhere:
    python3 verification/check_claims.py [--no-lean]
"""

import argparse
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CLAIMS = os.path.join(ROOT, 'doc', 'spec', 'claims.md')
ORACLE = os.path.join(ROOT, 'tests', 'pyminisketch.py')
LEAN_DIR = os.path.join(ROOT, 'verification')

STATUSES = {'proved', 'stated', 'tested', 'open', 'not-a-theorem'}
STD_AXIOMS = {'propext', 'Classical.choice', 'Quot.sound'}


def parse_claims(path):
    claims = []
    cur = None
    with open(path, encoding='utf-8') as f:
        for line in f:
            m = re.match(r'^### (C\d+): (.*)$', line)
            if m:
                cur = {'id': m.group(1), 'title': m.group(2).strip()}
                claims.append(cur)
                continue
            m = re.match(r'^- ([a-z_]+): (.*)$', line)
            if m and cur is not None:
                cur[m.group(1)] = m.group(2).strip()
    return claims


def oracle_names(path):
    with open(path, encoding='utf-8') as f:
        src = f.read()
    return set(re.findall(r'^(?:def|class) (\w+)', src, re.M)) | set(re.findall(r'^(\w+) = ', src, re.M))


def lean_axioms(names):
    """Map each theorem name to its axiom list (None if it does not exist)."""
    script = 'import Minisketch\n' + ''.join(f'#print axioms {n}\n' for n in names)
    proc = subprocess.run(['lake', 'env', 'lean', '--stdin'], input=script, cwd=LEAN_DIR,
                          capture_output=True, text=True)
    out = proc.stdout + proc.stderr
    result = {}
    for n in names:
        if re.search(rf"unknown (constant|identifier) '?{re.escape(n)}'?", out):
            result[n] = None
            continue
        m = re.search(rf"'{re.escape(n)}' depends on axioms: \[([^\]]*)\]", out)
        if m:
            result[n] = {a.strip() for a in m.group(1).split(',') if a.strip()}
        elif re.search(rf"'{re.escape(n)}' does not depend on any axioms", out):
            result[n] = set()
        else:
            result[n] = None
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--no-lean', action='store_true', help="skip the Lean checks")
    args = parser.parse_args()

    claims = parse_claims(CLAIMS)
    names = oracle_names(ORACLE)
    errors = []
    ids = [c['id'] for c in claims]
    if len(ids) != len(set(ids)):
        errors.append("duplicate claim ids")

    for c in claims:
        status = c.get('status')
        if status not in STATUSES:
            errors.append(f"{c['id']}: bad status {status!r}")
        oracle = c.get('oracle', '-')
        if oracle != '-' and not oracle.endswith('(planned)') and oracle.split('.')[0] not in names:
            errors.append(f"{c['id']}: oracle {oracle!r} not defined in tests/pyminisketch.py")
        if status in ('proved', 'stated') and c.get('lean', '-') == '-':
            errors.append(f"{c['id']}: {status} but no lean theorem")

    need = [c for c in claims if c.get('status') in ('proved', 'stated')]
    if need and not args.no_lean:
        axioms = lean_axioms([c['lean'] for c in need])
        for c in need:
            ax = axioms[c['lean']]
            if ax is None:
                errors.append(f"{c['id']}: Lean theorem {c['lean']} not found")
                continue
            if c['status'] == 'proved':
                allowed = set(STD_AXIOMS)
                if c.get('native_decide') == 'allowed':
                    allowed.add('Lean.ofReduceBool')
                bad = ax - allowed
                if bad:
                    errors.append(f"{c['id']}: proved, but {c['lean']} uses {sorted(bad)}")

    counts = {}
    for c in claims:
        counts[c.get('status')] = counts.get(c.get('status'), 0) + 1
    print(f"{len(claims)} claims: " + ", ".join(f"{v} {k}" for k, v in sorted(counts.items())))
    for e in errors:
        print("ERROR:", e)
    sys.exit(1 if errors else 0)


if __name__ == '__main__':
    main()
