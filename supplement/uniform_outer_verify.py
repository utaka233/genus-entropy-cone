"""Fraction-only verification of a uniform E_4 <= 1/60 certificate.
No LP solver, numerical arithmetic or polytope enumeration is used.
"""
from fractions import Fraction as Q
from pathlib import Path
import json
from itertools import combinations, permutations

def lin(terms):
 a=[0]*15
 for m,c in terms:
  if m:a[m-1]+=c
 return a

def mi(i,j,k=0):return lin([(i|k,1),(j|k,1),(k,-1),(i|j|k,-1)])
def add(*terms):return [sum(c*a[k] for c,a in terms) for k in range(15)]
def dot(a,b):return sum(Q(x)*Q(y) for x,y in zip(a,b))
def shannon():
 a=[lin([(15,1),(15^(1<<i),-1)]) for i in range(4)]
 for i,j in combinations(range(4),2):
  r=15^(1<<i)^(1<<j)
  for k in range(16):
   if k&r==k:a.append(mi(1<<i,1<<j,k))
 return a

def ingleton():
 a=[]
 for i,j in combinations(range(4),2):
  k,l=[x for x in range(4) if x not in (i,j)]
  i,j,k,l=[1<<x for x in (i,j,k,l)]
  a.append(add((1,mi(i,j,k)),(1,mi(i,j,l)),(1,mi(k,l)),(-1,mi(i,j))))
 return a

def zy():
 a=[]
 for p in permutations(range(4)):
  i,j,k,l=[1<<x for x in p]
  a.append(add((2,mi(i,j,k)),(1,mi(i,k,j)),(1,mi(j,k,i)),(1,mi(i,j,l)),(1,mi(k,l)),(-1,mi(i,j))))
 return sorted(set(tuple(x) for x in a))

def main():
 c=json.loads(Path(__file__).with_name('uniform_outer_certificates.json').read_text())
 S,J,Z=shannon(),ingleton(),zy()
 assert S==c['shannon_rows'] and J==c['ingleton_rows']
 assert [list(z) for z in Z]==c['zhang_yeung_rows']
 d=[Q(x) for x in c['repair_delta']]
 assert d==[Q(a,40) for a in J[0]]
 assert dot(J[0],J[0])==10 and max(abs(t) for t in d)==Q(1,40)
 assert J[0][-1]==0
 B=S+[[-a for a in J[0]]]
 assert len(c['repair_farkas_certificates'])==34
 for k,cert in enumerate(c['repair_farkas_certificates']):
  assert cert['target_index']==k
  a=(S+J)[k]
  target=[Q(x)-4*dot(a,d)*Q(y) for x,y in zip(a,J[0])]
  weights={int(i):Q(v) for i,v in cert['farkas'].items()}
  assert all(v>=0 for v in weights.values())
  assert all(sum(v*B[i][s] for i,v in weights.items())==target[s] for s in range(15))
 A=S+[list(z) for z in Z]
 weights={int(i):Q(v) for i,v in c['zy_bound_farkas_certificate'].items()}
 assert all(v>=0 for v in weights.values())
 assert all(sum(v*A[i][s] for i,v in weights.items())==Q(J[0][s])+Q(s==14,6) for s in range(15))
 x=[Q(y) for y in c['zy_bound_sharp_witness']]
 assert x[-1]==1 and dot(J[0],x)==-Q(1,6) and all(dot(a,x)>=0 for a in A)
 print('Exact repair proof passed for all 34 Shannon/Ingleton rows.')
 print('Exact Zhang-Yeung bound J >= -T/6 and sharp outer witness passed.')
 print('Using L4 subset C4, normalized entropy slice error <= 1/60.')
if __name__=='__main__':main()
