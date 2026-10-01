"""Exact polynomial certificates for S4-invariant facets of genus-2 cone."""
import pickle,itertools,json,math
from functools import reduce
from fractions import Fraction
import numpy as np
from scipy.spatial import ConvexHull
from genus2_compute import PERM_COORD
D=pickle.load(open('genus2_data/profiles_canonical.pkl','rb'))
keys=[x for x in D if x[-1]>1]
P=np.log(np.array(keys,dtype=float));P/=P[:,-1,None]
B=np.array([[1. if S.bit_count()==k else 0. for S in range(1,16)] for k in range(1,4)])
Q=P@B.T;H=ConvexHull(Q)

def exps(n):
 a=b=c=0
 while n%2==0:n//=2;a+=1
 while n%3==0:n//=3;b+=1
 while n%5==0:n//=5;c+=1
 assert n==1
 return a,b,c

def sym_profile(k):
 p=[exps(n) for n in k]
 return [np.array([sum(p[S-1][j] for S in range(1,16) if S.bit_count()==i) for j in range(2)],dtype=np.int64) for i in range(1,4)]+[np.array(p[-1][:2],dtype=np.int64)]

def add(a,b):
 n=max(len(a),len(b));c=np.zeros(n,dtype=np.int64);c[:len(a)]+=a;c[:len(b)]+=b;return c

def mul(a,b):return np.convolve(a,b)
def det(M):
 n=len(M);s=np.array([0],dtype=np.int64)
 for p in itertools.permutations(range(n)):
  prod=np.array([1],dtype=np.int64)
  for i,j in enumerate(p):prod=mul(prod,M[i][j])
  sign=(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
  s=add(s,sign*prod)
 return s

def primitive(v):
 cont=reduce(math.gcd,(abs(int(z)) for x in v for z in x))
 return [x//cont for x in v]

def evaluate(c,r):return sum(float(x)*r**j for j,x in enumerate(c))
def main():
 results=[]
 for e in sorted({tuple(np.round(e,10)) for e in H.equations}):
  e=np.array(e);ids=[int(i) for i in H.vertices if abs(Q[i]@e[:3]+e[3])<1e-7]
  if len(ids)!=3:continue
  rows=[sym_profile(keys[i]) for i in ids if i!=33]
  M=[[np.array([x],dtype=np.int64) for x in [1,3,3,1]]]+rows
  v=primitive([(-1)**j*det([[row[k] for k in range(4) if k!=j] for row in M]) for j in range(4)])
  if sum(evaluate(x,1.584962500721156)*float(y) for x,y in zip(v,-e))<0:v=[-x for x in v]
  vf=np.array([evaluate(x,1.584962500721156) for x in v]);slack=np.c_[Q,np.ones(len(Q))]@vf
  zeros=abs(slack)<1e-7;Z=np.vstack([P[zeros][:,p] for p in PERM_COORD]);rank=np.linalg.matrix_rank(Z,tol=1e-7)
  print('VERTS',ids,'ROW',[x.tolist() for x in v],'rank',rank,'minpositive',slack[~zeros].min(),flush=True)
  patterns={}
  for key in keys:
   eg=[exps(n) for n in key]
   if any(x[2] for x in eg):continue
   vals=sym_profile(key);pol=np.array([0],dtype=np.int64)
   for a,b in zip(v,vals):pol=add(pol,mul(a,b))
   c=tuple(int(x) for x in pol);patterns[c]=patterns.get(c,0)+1
  exactzeros=sum(num for c,num in patterns.items() if all(x==0 for x in c))
  print('polynomialpatterns',len(patterns),'exactzeros',exactzeros,flush=True)
  # Every polynomial positive throughout narrow rational interval containing log2(3).
  lo=Fraction(1584962500721156,10**15);hi=Fraction(1584962500721157,10**15)
  bad=[]
  for c,num in patterns.items():
   if not any(c):continue
   left=right=Fraction(0)
   for coeff in reversed(c):
    cand=(left*lo,left*hi,right*lo,right*hi);left=min(cand)+coeff;right=max(cand)+coeff
   if left<=0:bad.append((c,num,str(left),str(right)))
  print('bad interval',bad[:3],len(bad),flush=True)
  results.append({'vertices':ids,'row':[x.tolist() for x in v],'rank':int(rank),'exactzeros':exactzeros,'min_positive_numeric':float(slack[~zeros].min()),'polynomialpatterns':len(patterns),'interval_bad':bad})
 open('genus2_data/symmetric_facets.json','w').write(json.dumps(results,indent=2))

if __name__ == "__main__":
 main()
