"""Finite subgroup entropy enumeration for hyperbolic genus <=2.
Analysis script: exact integer joint support indices, numerical polyhedral routines.
"""
import itertools
from pathlib import Path
import numpy as np
from scipy.optimize import linprog

OUT = Path(__file__).resolve().parent / 'genus2_data'

def make_group(name, elems, op):
    lut={e:i for i,e in enumerate(elems)}
    mul=np.array([[lut[op(a,b)] for b in elems] for a in elems],dtype=np.int16)
    n=len(elems)
    ids=[i for i in range(n) if np.all(mul[i]==np.arange(n)) and np.all(mul[:,i]==np.arange(n))]
    assert len(ids)==1
    e=ids[0]
    assert all(np.any(mul[i]==e) and np.any(mul[:,i]==e) for i in range(n))
    for i in range(n):
        assert np.array_equal(mul[mul[i,:,None],np.arange(n)], mul[i,mul])
    return {'name':name,'elems':elems,'mul':mul,'identity':e}

def groups():
    yield make_group('C10',list(range(10)),lambda a,b:(a+b)%10)
    mats=[t for t in itertools.product(range(3),repeat=4) if (t[0]*t[3]-t[1]*t[2])%3]
    def mm(a,b):
        return ((a[0]*b[0]+a[1]*b[2])%3,(a[0]*b[1]+a[1]*b[3])%3,(a[2]*b[0]+a[3]*b[2])%3,(a[2]*b[1]+a[3]*b[3])%3)
    yield make_group('GL2F3',mats,mm)
    ts=list(itertools.product(range(6),range(2),range(2)))
    def pp(x,y):
        a,b,c=x;d,e,f=y
        return ((a+(-1 if b else 1)*d)%6,(b+e)%2,(c+f+b*d)%2)
    yield make_group('2D12',ts,pp)

def bit_elements(bits):
    while bits:
        b=bits&-bits; yield b.bit_length()-1; bits^=b

def generate_subgroup(mul,identity,H,x):
    gen=list(bit_elements(H))+[x]
    seen=1<<identity; todo=[identity]
    while todo:
        a=todo.pop()
        for g in gen:
            y=int(mul[a,g]); b=1<<y
            if not seen&b:
                seen|=b;todo.append(y)
    return seen

def enumerate_subgroups(G):
    mul=G['mul'];identity=G['identity'];n=len(mul)
    init=1<<identity; seen={init};todo=[init]
    while todo:
        H=todo.pop()
        for x in range(n):
            if not (H>>x)&1:
                J=generate_subgroup(mul,identity,H,x)
                if J not in seen: seen.add(J);todo.append(J)
    subs=sorted(seen,key=lambda h:(h.bit_count(),h))
    G['subgroups']=subs
    return subs

PERMS=list(itertools.permutations(range(4)))
PERM_COORD=[]
for p in PERMS:
    PERM_COORD.append(tuple(sum((1<<p[i]) for i in range(4) if S>>i&1)-1 for S in range(1,16)))

def canonical(v):
    return min(tuple(v[j] for j in perm) for perm in PERM_COORD)

def profile(quad,order):
    a,b,c,d=quad
    ints=(a,b,a&b,c,a&c,b&c,a&b&c,d,a&d,b&d,a&b&d,c&d,a&c&d,b&c&d,a&b&c&d)
    return tuple(order//h.bit_count() for h in ints)

