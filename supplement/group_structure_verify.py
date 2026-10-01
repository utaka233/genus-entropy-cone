"""Exact, dependency-free audit of the manuscript's order-24 carrier.

Run: python supplement/group_structure_verify.py
No SmallGroups identifier is assumed or certified by this Python audit.
"""

from itertools import product
from collections import Counter

elements = tuple(product(range(6), range(2), range(2)))
one = (0, 0, 0)


def mul(x, y):
    a, b, c = x
    d, e, f = y
    return ((a + (-1) ** b * d) % 6, (b + e) % 2, (c + f + b * d) % 2)


def power(x, n):
    out = one
    for _ in range(n):
        out = mul(out, x)
    return out


def generated(gens):
    out = {one}
    while True:
        new = out | {mul(x, y) for x in out for y in gens}
        if new == out:
            return out
        out = new


assert all(mul(one, x) == x == mul(x, one) for x in elements)
assert all(any(mul(x, y) == one == mul(y, x) for y in elements) for x in elements)
assert all(mul(mul(x, y), z) == mul(x, mul(y, z))
           for x in elements for y in elements for z in elements)
U, V, z = (0, 1, 0), (1, 0, 0), (0, 0, 1)
assert power(U, 2) == one and power(V, 6) == one and power(z, 2) == one
assert power(mul(U, V), 2) == z
assert all(mul(z, x) == mul(x, z) for x in elements)
assert len(generated([U, V, z])) == 24
N = generated([power(V, 2)])
K = generated([U, power(V, 3), z])
assert len(N) == 3 and len(K) == 8 and N & K == {one}
assert {mul(n, k) for n in N for k in K} == set(elements)
assert all({mul(mul(x, n), power(x, 23)) for n in N} == N for x in elements)
assert power(mul(U, power(V, 3)), 2) == z
assert len(generated([U, power(V, 3)])) == 8
assert mul(mul(U, power(V, 2)), U) == power(V, 4)
for x in (power(V, 3), z):
    assert mul(x, power(V, 2)) == mul(power(V, 2), x)
orders = Counter(next(j for j in range(1, 25) if power(x, j) == one) for x in elements)
assert orders == {1: 1, 2: 9, 3: 2, 4: 6, 6: 6}
print("PASS: order 24, group axioms, presentation, and exact C3 semidirect D8 splitting")
print("Complement involutions U,V^3; their product has square z and order 4")
print("D8 action on C3: U inverts, V^3 and z commute")
print("Element order counts:", dict(sorted(orders.items())))
