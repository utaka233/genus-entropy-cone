import itertools,json,math
from fractions import Fraction
import numpy as np

def add(a,b):
 n=max(len(a),len(b));c=[0]*n
 for i,x in enumerate(a):c[i]+=int(x)
 for i,x in enumerate(b):c[i]+=int(x)
 return c

def mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=int(x)*int(y)
 return c

def neg(a):return [-int(x) for x in a]
def det(M):
 n=len(M);s=[0]
 for p in itertools.permutations(range(n)):
  prod=[1]
  for i,j in enumerate(p):prod=mul(prod,M[i][j])
  sign=(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
  s=add(s,prod if sign==1 else neg(prod))
 return s

def cramer(M,rhs):
 denominator=det(M);num=[]
 for j in range(len(M)):
  N=[[rhs[i] if k==j else x for k,x in enumerate(row)] for i,row in enumerate(M)];num.append(det(N))
 return denominator,num

def interval(c,lo,hi):
 left=right=Fraction(0)
 for coeff in reversed(c):
  cand=(left*lo,left*hi,right*lo,right*hi);left=min(cand)+coeff;right=max(cand)+coeff
 return left,right

def log_bounds(z,N=50):
 z=Fraction(z);s=2*sum((z**(2*k+1)/Fraction(2*k+1) for k in range(N)),Fraction(0));tail=2*z**(2*N+1)/(Fraction(2*N+1)*(1-z*z));return s,s+tail
l2,u2=log_bounds(Fraction(1,3));l3,u3=log_bounds(Fraction(1,2));lo=Fraction(1584962500721156,10**15);hi=Fraction(1584962500721157,10**15)
assert lo<l3/u2 and u3/l2<hi
rows=[c['row'] for c in json.load(open('genus2_data/symmetric_facets.json'))];rows.append([[0],[0],[-1],[4]])
# 7 fromdeterminants plusupperplane; exacteightprojectedvertices.
V=[[[4],[6],[4],[1]],[[1],[3],[3],[1]],[[3],[9],[8],[2]],[[2,2],[7,5],[8,4],[2,1]],[[10,1],[24,6],[16,4],[4,1]],[[8,2],[23,6],[16,4],[4,1]],[[8,1],[22,4],[16,4],[4,1]],[[5,1],[13,4],[12,4],[3,1]]]
certs=[];feasible=[];unknown=[]
for inds in itertools.combinations(range(8),3):
 M=[[rows[i][j] for j in range(3)] for i in inds];rhs=[neg(rows[i][3]) for i in inds];D,N=cramer(M,rhs)
 if not any(D):certs.append({'indices':inds,'type':'singular_exact'});continue
 a,b=interval(D,lo,hi)
 if a<=0<=b:unknown.append(('det',inds,D));continue
 ds=1 if a>0 else -1;slacks=[];ok=True
 for row in rows:
  s=mul(row[3],D)
  for j in range(3):s=add(s,mul(row[j],N[j]))
  if not any(s):slacks.append(0);continue
  a,b=interval(s,lo,hi)
  if a<=0<=b:unknown.append(('slack',inds,s));ok=False;break
  sign=ds*(1 if a>0 else -1);slacks.append(sign)
  if sign<0:ok=False
 if ok:
  matches=[]
  for k,v in enumerate(V):
   diffs=[add(mul(N[j],v[3]),neg(mul(D,v[j]))) for j in range(3)]
   if all(not any(s) for s in diffs):matches.append(k)
  if len(matches)!=1:unknown.append(('match',inds,matches));continue
  feasible.append(matches[0]);certs.append({'indices':inds,'type':'feasible','vertex':matches[0],'slacks':slacks})
 else:certs.append({'indices':inds,'type':'infeasible','slacks':slacks})
print('unknown',unknown,'feasibleclusters',sorted(set(feasible)),'feasibletriples',len(feasible),'allcases',len(certs),flush=True)
assert not unknown and set(feasible)==set(range(8)) and len(certs)==56
# Boundedness: positive combination offour facet normals equalszero; three independent.
bounded=None
for inds in itertools.combinations(range(8),4):
 M=[[rows[inds[j]][i] for j in range(3)] for i in range(3)];rhs=[neg(rows[inds[3]][i]) for i in range(3)];D,N=cramer(M,rhs)
 if not any(D):continue
 a,b=interval(D,lo,hi)
 if a<=0<=b:continue
 ds=1 if a>0 else -1
 if all((interval(c,lo,hi)[0]>0 if ds==1 else interval(c,lo,hi)[1]<0) for c in N):
  bounded={'indices':inds,'denominator':D,'numerators':N};break
print('boundednesscert',bounded,flush=True);assert bounded
out={'rational_r_interval':[str(lo),str(hi)],'log2bounds':[str(l2),str(u2)],'log3bounds':[str(l3),str(u3)],'vertices':V,'facets':rows,'intersections':certs,'boundedness':bounded,'proof':'exactintegerpolynomials and Fraction interval arithmetic'}
json.dump(out,open('genus2_data/symmetric_slice_certificate.json','w'),indent=2)
