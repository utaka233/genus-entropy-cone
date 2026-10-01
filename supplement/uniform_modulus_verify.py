"""Exact arithmetic checks for four-atom entropy and low-genus RH defects.
The published finite-group classification is an external theorem, not rerun.
"""
from fractions import Fraction as Q
from collections import defaultdict
from uniform_outer_verify import ingleton
import math


def factor(n):
    out = defaultdict(int)
    p = 2
    while p * p <= n:
        while n % p == 0:
            out[p] += 1
            n //= p
        p += 1
    if n > 1:
        out[n] += 1
    return out


def entropy_coefficients(probabilities):
    out = defaultdict(Q)
    for q in probabilities:
        for prime, k in factor(q.numerator).items():
            out[prime] -= q * k
        for prime, k in factor(q.denominator).items():
            out[prime] += q * k
    return {p:v for p,v in out.items() if v}


def main():
    atoms = [(0,0,0,0),(1,0,1,0),(1,0,0,1),(1,1,1,1)]
    probabilities = [Q(1,3),Q(1,6),Q(1,6),Q(1,3)]
    h = []
    for mask in range(1,16):
        marginal = defaultdict(Q)
        for a,p in zip(atoms, probabilities):
            marginal[tuple(a[i] for i in range(4) if mask & (1<<i))] += p
        h.append(entropy_coefficients(marginal.values()))
    j = ingleton()[0]
    value = defaultdict(Q)
    for c,coefficients in zip(j,h):
        for prime,q in coefficients.items():
            value[prime] += c*q
    assert {p:q for p,q in value.items() if q} == {2:Q(3),3:Q(-2)}
    assert h[-1] == {2:Q(1,3),3:Q(1)}
    assert sum(abs(c) for c in j) == 10 and j[-1] == 0
    expected = {120:[(2,4,5)],144:[(2,3,8)],192:[],240:[]}
    for order, triples in expected.items():
        defect = Q(6,order)
        found = [(2,b,c) for b in (3,4) for c in range(b,order+1)
                 if 1-Q(1,2)-Q(1,b)-Q(1,c) == defect]
        assert found == triples
    r = math.log(3,2)
    lower2 = (r-1)**2 / (6*(3*r*r+r-6))
    lower3 = math.log(Q(9,8)) / (10*(math.log(3)+math.log(2)/3))
    c = math.log(Q(16,15))/math.log(120)
    lower4 = lower3-c/10
    print('Exact four-atom result: J=3log2-2log3=log(8/9), T=log3+log2/3.')
    print('Exact RH triangle possibilities passed for orders120,144,192,240.')
    print({'E2_lower_approx':lower2,'E3_lower_approx':lower3,
           'E4_lower_approx':lower4,'E4_upper_exact':'7/440'})


if __name__ == '__main__':
    main()
