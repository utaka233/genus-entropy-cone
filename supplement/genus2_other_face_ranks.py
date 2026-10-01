import json
from fractions import Fraction
import numpy as np
from polynomial_tools import trim,add,mul,sub,divexact,interval,log_bounds
from genus2_certify import exponents

indices = np.load('genus2_data/extreme_rays_verified.npz')['indices']
assert len(indices) == 127
R = []
for vector in indices:
    a,b,c = exponents(vector)
    assert not np.any(c)
    R.append((a.tolist(), b.tolist()))
# The three invariant-supportingplanes that are not fullcone facets.
rows=[[[0],[0],[-1],[4]],[[0],[-1],[3],[-6]],[[-5,-1],[6,2],[-3,-3],[-4,4]]]
l2,u2=log_bounds(Fraction(1,3));l3,u3=log_bounds(Fraction(1,2));lo=l3/u2;hi=u3/l2
out=[]
for row in rows:
 faces=[]
 for j,(A,B) in enumerate(R):
  v=[[sum(A[S-1] for S in range(1,16) if S.bit_count()==i),sum(B[S-1] for S in range(1,16) if S.bit_count()==i)] for i in range(1,4)]+[[A[-1],B[-1]]]
  s=[0]
  for a,b in zip(row,v):s=add(s,mul(a,b))
  if not any(s):faces.append(j)
  else:
   sl,su=interval(s,lo,hi)
   assert sl>0,(j,s,sl,su)
 M=[[[a,b] for a,b in zip(*R[j])] for j in faces];m=len(M);prev=[1];rank=0;cols=[]
 for j in range(15):
  sel=next((i for i in range(rank,m) if any(M[i][j])),None)
  if sel is None:continue
  M[rank],M[sel]=M[sel],M[rank];p=M[rank][j]
  for i in range(rank+1,m):
   q=M[i][j]
   for k in range(j+1,15):M[i][k]=divexact(sub(mul(M[i][k],p),mul(q,M[rank][k])),prev)
   M[i][j]=[0]
  prev=p;rank+=1;cols.append(j)
  if rank==m:break
 assert all(not any(x) for row2 in M[rank:] for x in row2)
 dl,du=interval(prev,lo,hi);assert dl>0 or du<0
 print('row',row,'zeroextremerays',len(faces),'exactrank',rank,flush=True)
 out.append({'row':row,'zero_rays':faces,'exact_rank':rank,'pivot_columns':cols,'last_pivot_determinant':prev,'rank_upper_zero_rows_verified':True})
assert [entry['exact_rank'] for entry in out]==[11,9,11]
json.dump(out,open('genus2_data/other_symmetric_face_ranks.json','w'),indent=2)
print('PASS: the three remaining symmetric inequalities have full-cone equality-face dimensions 11,9,11.',flush=True)
