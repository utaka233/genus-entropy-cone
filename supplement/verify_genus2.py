"""Reproduce the exact genus-two entropy-cone certificates.

Run from the extracted archive directory: python verify_genus2.py
Numerical LP/QR only selects witnesses; integers and Fractions verify them.
"""
import itertools
import json
import math
import pickle
import runpy
from fractions import Fraction
from pathlib import Path

import numpy as np

from genus2_compute import groups, enumerate_subgroups, profile, canonical, PERM_COORD
from genus2_certify import exponents, exact_memberships, exact_exposures

root = Path(__file__).resolve().parent
import os
os.chdir(root)
out = root / 'genus2_data'

canonical_profiles = {}
carriers = {}
for group in groups():
    subs = enumerate_subgroups(group)
    carriers[group['name']] = group
    for ambient in subs:
        inside = [h for h in subs if h & ambient == h]
        for quad in itertools.combinations_with_replacement(inside, 4):
            v = profile(quad, ambient.bit_count())
            w = canonical(v)
            canonical_profiles.setdefault(w, v)
assert len(canonical_profiles) == 6154
stored = pickle.load((out / 'profiles_canonical.pkl').open('rb'))
assert set(canonical_profiles) == set(stored)

E = np.load(out / 'extreme_rays_verified.npz')['indices']
extreme = {tuple(map(int, v)) for v in E}
assert len(extreme) == 127
assert {tuple(v[j] for j in p) for v in extreme for p in PERM_COORD} == extreme
assert len({canonical(v) for v in extreme}) == 18
certificates = json.load((out / 'extreme_rays_certificates.json').open())
assert {tuple(c['indices']) for c in certificates} == extreme
for cert in certificates:
    group = carriers[cert['carrier']]
    ambient = cert['ambient_bitmask']
    quad = cert['subgroups_bitmasks']
    assert ambient in group['subgroups']
    assert all(h in group['subgroups'] and h & ambient == h for h in quad)
    assert profile(quad, ambient.bit_count()) == tuple(cert['indices'])

bases = pickle.load((out / 'lp_basis_certificates.pkl').open('rb'))
expected = {v for v in canonical_profiles if v[-1] > 1 and not np.any(exponents(v)[2])}
assert len(bases) == len(expected) == 6123
assert {tuple(v) for v, _, _ in bases} == expected
exact_memberships()
exact_exposures()

minimum_R = {}
for v in canonical_profiles:
    T = v[-1]
    if T == 1:
        continue
    numerator = math.prod(v[S-1] for S in range(1, 16) if S.bit_count() in (1, 3))
    denominator = math.prod(v[S-1] for S in range(1, 16) if S.bit_count() == 2)
    R = Fraction(numerator, denominator)
    assert R >= Fraction(4, 3)
    minimum_R[T] = min(minimum_R.get(T, R), R)
(out / 'sharp_support_exact_table.json').write_text(
    json.dumps({T: str(R) for T, R in sorted(minimum_R.items())}, indent=2))

# Sequential: these scripts write dependent certificate files.
runpy.run_path(str(root / 'genus2_symmetric.py'), run_name='__main__')
runpy.run_path(str(root / 'genus2_rank_cert.py'), run_name='__main__')
runpy.run_path(str(root / 'genus2_slice_cert.py'), run_name='__main__')
facets = json.load((out / 'symmetric_facets.json').open())
assert sum(bool(c.get('rank_proved')) for c in facets) == 5
assert all(not c['interval_bad'] for c in facets)
slice_cert = json.load((out / 'symmetric_slice_certificate.json').open())
assert len(slice_cert['vertices']) == len(slice_cert['facets']) == 8

# Known 35-ray theorem for Shannon + Ingleton; verify the supplied representatives
# and the precise genus comparison with exact integer arithmetic.
from genus2_polyhedral import shannon_ingleton
from genus2_certify import integer_rank
comparison = json.load((out / 'ingleton_genus_comparison.json').open())
inequalities = shannon_ingleton().astype(int)
rank_rays = {tuple(v) for v in comparison['rays35']}
assert len(rank_rays) == 35
for v in rank_rays:
    rv = np.array(v, dtype=int)
    assert np.all(inequalities @ rv >= 0)
    assert integer_rank(inequalities[inequalities @ rv == 0]) == 14
binary = {tuple(int(n).bit_length()-1 for n in v) for v in extreme
          if all(int(n) & (int(n)-1) == 0 for n in v)}
u24 = tuple(min(S.bit_count(), 2) for S in range(1, 16))
u34 = tuple(min(S.bit_count(), 3) for S in range(1, 16))
coarse = {tuple(min(3, S.bit_count()+int(bool(S & (1 << i))))
                for S in range(1, 16)) for i in range(4)}
assert rank_rays - binary == {u24, u34} | coarse
print('PASS: 6154 canonical profiles; 127 extreme rays; 18 S4 orbits;')
print('      five full-cone facets; complete symmetric slice with eight facets.')
