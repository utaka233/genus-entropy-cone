"""Exact polynomial and rational-interval arithmetic for face-rank checks."""
from fractions import Fraction

def trim(p):
 p=list(map(int,p))
 while len(p)>1 and p[-1]==0:p.pop()
 return p

def add(a,b):
 c=[0]*max(len(a),len(b))
 for i,x in enumerate(a):c[i]+=x
 for i,x in enumerate(b):c[i]+=x
 return trim(c)

def neg(a):return [-x for x in a]

def sub(a,b):return add(a,neg(b))

def mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  if not x:continue
  for j,y in enumerate(b):c[i+j]+=x*y
 return trim(c)

def divexact(a,b):
 a=trim(a);b=trim(b);assert any(b)
 if not any(a):return [0]
 assert len(a)>=len(b)
 c=[0]*(len(a)-len(b)+1)
 while any(a) and len(a)>=len(b):
  n=len(a)-len(b);q,rest=divmod(a[-1],b[-1]);assert rest==0,(a,b)
  c[n]=q
  for i,x in enumerate(b):a[i+n]-=q*x
  a=trim(a)
 assert not any(a),(a,b)
 return trim(c)

def interval(c,lo,hi):
 left=right=Fraction(0)
 for x in reversed(c):
  vals=(left*lo,left*hi,right*lo,right*hi);left=min(vals)+x;right=max(vals)+x
 return left,right

def log_bounds(z,N=80):
 z=Fraction(z);s=2*sum((z**(2*k+1)/Fraction(2*k+1) for k in range(N)),Fraction(0));tail=2*z**(2*N+1)/(Fraction(2*N+1)*(1-z*z));return s,s+tail
