#!/usr/bin/env python3
"""Exact obstruction to U_{2,4} for the relevant genus-three groups.
Uses integer multiplication tables and subgroup bitsets; no numerical LP.
"""
import itertools
import json
from collections import Counter
from pathlib import Path
import time

def group(name, elems, op):
    lookup={v:i for i,v in enumerate(elems)}
    mul=[[lookup[op(a,b)] for b in elems] for a in elems]
    n=len(elems)
    identities=[i for i in range(n) if all(mul[i][j]==j==mul[j][i] for j in range(n))]
    assert len(identities)==1
    e=identities[0]
    assert all(any(mul[a][b]==e==mul[b][a] for b in range(n)) for a in range(n))
    assert all(mul[mul[a][b]][c]==mul[a][mul[b][c]] for a in range(n) for b in range(n) for c in range(n))
    return {'name':name, 'elements':elems, 'mul':mul, 'identity':e}

def permutation_product(p,q):
    return tuple(p[q[i]] for i in range(len(p)))

def relevant_groups():
    # Projective monomial matrices in three variables with fourth roots of unity.
    # Diagonal exponents have first coordinate zero; scalar matrices are factored out.
    perms=list(itertools.permutations(range(3)))
    elems=[(a,b,p) for a in range(4) for b in range(4) for p in perms]
    def monomial(t,u):
        a,b,p=t;c,d,q=u
        z=(0,a,b);w=(0,c,d)
        v=tuple((z[i]+w[p[i]])%4 for i in range(3))
        return ((v[1]-v[0])%4,(v[2]-v[0])%4,permutation_product(q,p))
    # With the action convention above the resulting composition permutation is q o p.
    yield group('C4^2:S3 (96,64)', elems, monomial)
    # C4 central product SL(2,3): quotient pairs (A,t) by (-I,2).
    mats=[t for t in itertools.product(range(3),repeat=4) if (t[0]*t[3]-t[1]*t[2])%3==1]
    def canonical(A,t):
        B=tuple((-x)%3 for x in A)
        return min((A,t%4),(B,(t+2)%4))
    elems=sorted({canonical(A,t) for A in mats for t in range(4)})
    def prod(a,b):
        A,t=a;B,u=b
        C=((A[0]*B[0]+A[1]*B[2])%3,(A[0]*B[1]+A[1]*B[3])%3,
           (A[2]*B[0]+A[3]*B[2])%3,(A[2]*B[1]+A[3]*B[3])%3)
        return canonical(C,t+u)
    yield group('C4 central product SL2(F3) (48,33)',elems,prod)
    # Automorphisms r:(x,y)->(zeta_8*x,y), s:(x,y)->(1/x,i*y/x^4), z:(x,y)->(x,-y)
    elems=list(itertools.product(range(8),range(2),range(2)))
    def hyper(t,u):
        a,b,c=t;d,e,f=u
        return ((a+(-1 if b else 1)*d)%8,(b+e)%2,(c+f+b*d+b*e)%2)
    yield group('Hyperelliptic polygon group (32,9)',elems,hyper)
    perms=list(itertools.permutations(range(4)))
    yield group('C2 x S4 (48,48)',[(p,t) for p in perms for t in range(2)],
                lambda a,b:(permutation_product(a[0],b[0]),(a[1]+b[1])%2))

def bits(mask):
    while mask:
        b=mask&-mask;yield b.bit_length()-1;mask^=b

def generate(mul,e,gens):
    seen=1<<e;todo=[e]
    while todo:
        a=todo.pop()
        for g in gens:
            b=mul[a][g]
            if not ((seen>>b)&1):seen|=1<<b;todo.append(b)
    return seen

def enumerate_subgroups(G):
    mul=G['mul'];e=G['identity'];n=len(mul)
    trivial=1<<e
    discovered={trivial:()};todo=[trivial]
    while todo:
        H=todo.pop();gens=discovered[H]
        for x in range(n):
            if not ((H>>x)&1):
                J=generate(mul,e,gens+(x,))
                if J not in discovered:discovered[J]=gens+(x,);todo.append(J)
    return sorted(discovered,key=lambda H:(H.bit_count(),H))

def verify_model_presentations(G):
    mul=G['mul'];e=G['identity'];n=len(mul)
    def power(a,k):
        b=e
        for _ in range(k):b=mul[b][a]
        return b
    def inverse(a):return next(b for b in range(n) if mul[a][b]==e==mul[b][a])
    def conj(a,b):return mul[mul[a][b]][inverse(a)]
    if n==48 and 'central product' in G['name']:
        I=(1,0,0,1)
        a=G['elements'].index((I,1))
        a2=power(a,2)
        four=[b for b in range(n) if power(b,2)==a2 and G['elements'][b][1]==0]
        three=[d for d in range(n) if d!=e and power(d,3)==e]
        candidates=[]
        for b,c,d in itertools.product(four,four,three):
            if (conj(c,b)==mul[a2][b] and conj(d,b)==mul[mul[a2][b]][c]
                    and conj(d,c)==b):
                generated=generate(mul,e,(a,b,c,d))
                if generated.bit_count()==48:candidates.append((a,b,c,d))
        assert candidates, 'No match for the documented (48,33) presentation'
        a,b,c,d=candidates[0]
        assert power(a,4)==e and all(conj(a,t)==t for t in (b,c,d))
        assert generate(mul,e,(a,b,c)).bit_count()==16
        G['presentation_generators']=[G['elements'][i] for i in (a,b,c,d)]
    if n==32:
        r=G['elements'].index((1,0,0));s=G['elements'].index((0,1,0));z=G['elements'].index((0,0,1))
        assert power(r,8)==e and power(z,2)==e and power(s,2)==z
        assert conj(s,r)==mul[z][inverse(r)]
        r2=power(r,2);rs=mul[r][s]
        assert power(rs,2)==e and conj(rs,r2)==inverse(r2)
        D=generate(mul,e,(r2,rs))
        assert D.bit_count()==8 and not ((D>>z)&1)
        assert generate(mul,e,(r2,rs,z)).bit_count()==16
        G['D8xC2_generators']=[G['elements'][i] for i in (r2,rs,z)]

def obstruction(G):
    subgroups=enumerate_subgroups(G)
    rows=[]
    for ambient in subgroups:
        order=ambient.bit_count()
        if order%16:continue
        vertices=[H for H in subgroups if H.bit_count()==order//4 and H&ambient==H]
        colors={H&J for H,J in itertools.combinations(vertices,2) if (H&J).bit_count()==order//16}
        maximum=int(bool(vertices));checked=0
        for size in (2,3,4):
            for family in itertools.combinations(vertices,size):
                if size==4:checked+=1
                K=family[0]&family[1]
                if K.bit_count()!=order//16:continue
                if all(H&J==K for H,J in itertools.combinations(family,2)):
                    assert size<4, (G['name'],order,family,K)
                    maximum=max(maximum,size)
        rows.append({'ambient_order':order,'index_four_subgroups':len(vertices),
                     'intersection_colors':len(colors),'maximum_compatible_family':maximum,
                     'quadruples_checked':checked})
    return {'name':G['name'],'order':len(G['mul']), 'subgroup_count':len(subgroups),
            'subgroup_orders':dict(sorted(Counter(H.bit_count() for H in subgroups).items())),
            'eligible_ambient_count':len(rows),'by_ambient_order':dict(sorted(Counter(r['ambient_order'] for r in rows).items())),
            'maximum_compatible_family':max(r['maximum_compatible_family'] for r in rows),
            'presentation_generators':G.get('presentation_generators'),
            'D8xC2_generators':G.get('D8xC2_generators'),
            'rows':rows}

def verify_controls():
    # Four pairwise transverse two-dimensional subspaces in F_2^4
    # (the lines of slopes 0,1,alpha,infinity in F_4^2).
    control=group('positive control C2^4',list(range(16)),lambda a,b:a^b)
    alpha=(0,2,3,1)
    lines=[sum(1<<(a+4*b) for a,b in pts) for pts in (
        [(a,0) for a in range(4)],[(a,a) for a in range(4)],
        [(a,alpha[a]) for a in range(4)],[(0,b) for b in range(4)])]
    subgroups=enumerate_subgroups(control)
    assert all(H in subgroups and H.bit_count()==4 for H in lines)
    assert all(H&J==1 for H,J in itertools.combinations(lines,2))
    assert len(enumerate_subgroups(control))==67

def main():
    start=time.time();summary=[]
    verify_controls()
    for G in relevant_groups():
        print('Checking',G['name'],flush=True)
        verify_model_presentations(G)
        s=obstruction(G);summary.append(s)
        print('PASS',s['subgroup_count'],'subgroups;',s['eligible_ambient_count'],'eligible ambient subgroups;',
              'largest compatible family',s['maximum_compatible_family'],flush=True)
    out=Path(__file__).with_name('genus3_obstruction_certificates.json')
    out.write_text(json.dumps(summary,indent=2))
    print('PASS: no genus-three carrier has the required four-subgroup configuration.',flush=True)
    print('elapsed seconds',round(time.time()-start,3),flush=True)
if __name__=='__main__':main()
