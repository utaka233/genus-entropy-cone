import json,pickle,itertools,math
from fractions import Fraction
import numpy as np
from scipy.linalg import qr
from genus2_compute import PERM_COORD
from genus2_symmetric import exps,add,mul

def det_integer(M):
 M=[list(map(int,row)) for row in M];n=len(M);prev=1;sign=1
 for k in range(n-1):
  if M[k][k]==0:
   j=next((j for j in range(k+1,n) if M[j][k]),None)
   if j is None:return 0
   M[k],M[j]=M[j],M[k];sign=-sign
  p=M[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):M[i][j]=(M[i][j]*p-M[i][k]*M[k][j])//prev
  for i in range(k+1,n):M[i][k]=0
  prev=p
 return sign*M[-1][-1]

def det_poly(A,B):
 n=len(A);values=[det_integer(A+t*B) for t in range(n+1)]
 differences=[];cur=values
 while cur:
  differences.append(cur[0]);cur=[b-a for a,b in zip(cur,cur[1:])]
 c=[Fraction(0)]*(n+1);basis=[1]
 for k,v in enumerate(differences):
  for j,b in enumerate(basis):c[j]+=Fraction(v*b,math.factorial(k))
  nxt=[0]*(len(basis)+1)
  for j,b in enumerate(basis):nxt[j]-=k*b;nxt[j+1]+=b
  basis=nxt
 assert all(x.denominator==1 for x in c)
 return list(map(int,c))

def interval(c,lo,hi):
 left=right=Fraction(0)
 for coeff in reversed(c):
  cand=(left*lo,left*hi,right*lo,right*hi);left=min(cand)+coeff;right=max(cand)+coeff
 return left,right

r0=math.log2(3);lo=Fraction(1584962500721156,10**15);hi=Fraction(1584962500721157,10**15)
D=pickle.load(open('genus2_data/profiles_canonical.pkl','rb'));keys=[k for k in D if k[-1]>1 and not any(exps(n)[2] for n in k)]
AB=np.array([[[exps(n)[j] for n in k] for j in range(2)] for k in keys],dtype=int)
certs=json.load(open('genus2_data/symmetric_facets.json'))
for cert in certs:
 if cert['rank']!=14:continue
 row=cert['row'];vf=np.array([sum(x*r0**i for i,x in enumerate(c)) for c in row])
 Ps=(AB[:,0]+r0*AB[:,1]);Q=np.c_[*[Ps[:,[S-1 for S in range(1,16) if S.bit_count()==i]].sum(axis=1) for i in range(1,4)],Ps[:,-1]]
 zeros=np.abs(Q@vf)<1e-7
 ZA=np.vstack([AB[zeros,0][:,p] for p in PERM_COORD]);ZB=np.vstack([AB[zeros,1][:,p] for p in PERM_COORD]);Z=ZA+r0*ZB
 _,_,pi=qr(Z.T,pivoting=True,mode='economic');rowids=pi[:14];R=Z[rowids];_,_,pj=qr(R,pivoting=True,mode='economic');cols=pj[:14]
 A=ZA[rowids][:,cols];B=ZB[rowids][:,cols];c=det_poly(A,B);L,U=interval(c,lo,hi)
 print('row',row,'minorpoly',c,'interval',[float(L),float(U)],'success',L>0 or U<0,flush=True)
 cert['rank_minor_full_A']=ZA[rowids].tolist();cert['rank_minor_full_B']=ZB[rowids].tolist();cert['rank_minor_columns']=cols.tolist();cert['rank_minor_coefficients']=c;cert['rank_minor_interval']=[str(L),str(U)];cert['rank_minor_A']=A.tolist();cert['rank_minor_B']=B.tolist();cert['rank_proved']=bool(L>0 or U<0)
json.dump(certs,open('genus2_data/symmetric_facets.json','w'),indent=2)
