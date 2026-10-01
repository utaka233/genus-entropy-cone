import numpy as np,pickle,time,json
from scipy.optimize import linprog
from genus2_compute import OUT

def exponents(v):
 a=[];b=[];c=[]
 for k in v:
  x=int(k);aa=bb=cc=0
  while x%2==0:x//=2;aa+=1
  while x%3==0:x//=3;bb+=1
  while x%5==0:x//=5;cc+=1
  assert x==1
  a.append(aa);b.append(bb);c.append(cc)
 return np.array(a,dtype=np.int64),np.array(b,dtype=np.int64),np.array(c,dtype=np.int64)


# Exact proof support: all arithmetic below uses Python integers/Fractions.
from fractions import Fraction

def log_interval_integer(x,N=120):
 z=Fraction(x-1,x+1);z2=z*z;term=z;s=Fraction(0)
 for k in range(N):
  s+=term/Fraction(2*k+1);term*=z2
 lower=2*s
 error=2*term/(Fraction(2*N+1)*(1-z2))
 return lower,lower+error

def t_interval():
 l2,u2=log_interval_integer(2);l3,u3=log_interval_integer(3)
 return l3/u2,u3/l2

def integer_rank(M):
 # Bareiss fraction-free row elimination for rectangular integral matrices.
 M=[list(map(int,row)) for row in M];nr=len(M);nc=len(M[0]);r=0;prev=1
 for col in range(nc):
  piv=next((i for i in range(r,nr) if M[i][col]),None)
  if piv is None:continue
  M[r],M[piv]=M[piv],M[r];pivot=M[r][col]
  for i in range(r+1,nr):
   factor=M[i][col]
   for j in range(col+1,nc):
    num=pivot*M[i][j]-factor*M[r][j]
    assert num%prev==0
    M[i][j]=num//prev
   M[i][col]=0
  prev=pivot;r+=1
  if r==nr:break
 return r

def inverse_fraction(M):
 n=len(M);rows=[[Fraction(x) for x in row]+[Fraction(int(i==j)) for j in range(n)] for i,row in enumerate(M)]
 for k in range(n):
  piv=next(i for i in range(k,n) if rows[i][k])
  rows[k],rows[piv]=rows[piv],rows[k]
  div=rows[k][k];rows[k]=[x/div for x in rows[k]]
  for i in range(n):
   if i==k:continue
   f=rows[i][k]
   if f:rows[i]=[x-f*y for x,y in zip(rows[i],rows[k])]
 return [row[n:] for row in rows]

def exact_memberships():
 from scipy.linalg import qr
 start=time.time()
 E=np.load(OUT/'extreme_rays_verified.npz')['indices']
 AA=[];BB=[]
 for v in E:
  a,b,c=exponents(v);assert not np.any(c);AA.append(a);BB.append(b)
 AA=np.array(AA).T;BB=np.array(BB).T
 supp=pickle.load((OUT/'lp_basis_certificates.pkl').open('rb'))
 tl,tu=t_interval();t0=Fraction('1.584962500721156181453738943947816508759814407692481060455752654541098227')
 eps=max(abs(tl-t0),abs(tu-t0))
 assert float(eps)<1e-60
 checked=0;worst=Fraction(1);polyevals=0
 for idx,(v,inds,_) in enumerate(supp):
  a,b,c=exponents(v);r=len(inds);A=AA[:,inds];B=BB[:,inds]
  # Each augmented (r+1)-minor has degree <=r+1. r+2 exact zeros prove itidentically zero.
  augA=np.column_stack([A,a]);augB=np.column_stack([B,b])
  for tval in range(r+2):
   assert integer_rank(augA+tval*augB)<=r,(idx,tval)
   polyevals+=1
  # An invertible row minor selects a unique exact solution at true t.
  _,_,rows=qr((A+B*float(t0)).T,pivoting=True,mode='economic');rows=rows[:r]
  Ar=A[rows];Br=B[rows];ar=a[rows];br=b[rows]
  M=[[Fraction(int(Ar[i,j]))+int(Br[i,j])*t0 for j in range(r)] for i in range(r)]
  R=inverse_fraction(M)
  y=[Fraction(int(ar[i]))+int(br[i])*t0 for i in range(r)]
  x=[sum(R[i][j]*y[j] for j in range(r)) for i in range(r)]
  C=[[sum(R[i][k]*int(Br[k,j]) for k in range(r)) for j in range(r)] for i in range(r)]
  d=[sum(R[i][k]*int(br[k]) for k in range(r))-sum(C[i][j]*x[j] for j in range(r)) for i in range(r)]
  cnorm=max(sum(abs(z) for z in row) for row in C)
  dnorm=max(abs(z) for z in d)
  assert eps*cnorm<1
  error=eps*dnorm/(1-eps*cnorm)
  assert min(x)>error,(idx,float(min(x)),float(error))
  worst=min(worst,min(x)-error);checked+=1
  if checked%1000==0:print('exactmembership',checked,'t',round(time.time()-start,1),flush=True)
 # C10 profiles split into prime2 andprime5 common-variable rays, allpure2 rankonerayspresent.
 can=pickle.load((OUT/'profiles_canonical.pkl').open('rb'));common={tuple(a) for a,b,c in map(exponents,E) if not np.any(b) and max(a)==1}
 cyclic=0
 for v in can:
  a,b,c=exponents(v)
  if np.any(c):
   assert not np.any(b)
   assert not np.any(a) or tuple(a) in common
   assert tuple(c) in common
   cyclic+=1
 result={'exact_memberships':checked,'cyclic_prime_split':cyclic,'polynomial_rank_evaluations':polyevals,'minimum_positive_lower_bound_approx':float(worst),'t_lower':str(tl),'t_upper':str(tu),'seconds':time.time()-start}
 (OUT/'exact_membership_summary.json').write_text(json.dumps(result,indent=2))
 print('EXACT MEMBERSHIP SUCCESS',checked,cyclic,'time',time.time()-start,flush=True)
 return result


def interval_add(x,y):return x[0]+y[0],x[1]+y[1]
def interval_mul(x,y):
 z=[a*b for a in x for b in y];return min(z),max(z)
def exact_exposures():
 start=time.time();E=np.load(OUT/'extreme_rays_verified.npz')['indices'];H=np.log(E.astype(float));P=H/H[:,-1,None]
 AA=[];BB=[]
 for v in E:
  a,b,c=exponents(v);assert not np.any(c);AA.append(a);BB.append(b)
 AA=np.array(AA);BB=np.array(BB);tl,tu=t_interval();tiv=(tl,tu);cert=[];minmargin=1
 for i,p in enumerate(P):
  others=np.delete(P,i,axis=0)
  Aub=np.column_stack([-(others-p),np.ones(len(others))])
  obj=np.zeros(16);obj[-1]=-1
  res=linprog(obj,A_ub=Aub,b_ub=np.zeros(len(others)),bounds=[(-1,1)]*15+[(None,None)],method='highs')
  assert res.success and res.x[-1]>1e-9
  D=10**12;c=np.rint(res.x[:15]*D).astype(np.int64)
  ac=[sum(int(c[k])*int(a[k]) for k in range(15)) for a in AA]
  bc=[sum(int(c[k])*int(b[k]) for k in range(15)) for b in BB]
  for j in range(len(E)):
   if i==j:continue
   # Difference numerator c·Pj-c·Pi, withpositive denominators.
   coef0=ac[j]*int(AA[i,-1])-ac[i]*int(AA[j,-1])
   coef1=ac[j]*int(BB[i,-1])+bc[j]*int(AA[i,-1])-ac[i]*int(BB[j,-1])-bc[i]*int(AA[j,-1])
   coef2=bc[j]*int(BB[i,-1])-bc[i]*int(BB[j,-1])
   val=interval_add((Fraction(coef0),Fraction(coef0)),interval_add(interval_mul((Fraction(coef1),Fraction(coef1)),tiv),interval_mul((Fraction(coef2),Fraction(coef2)),interval_mul(tiv,tiv))))
   assert val[0]>0,(i,j,float(val[0]))
  minmargin=min(minmargin,res.x[-1]);cert.append({'ray':i,'objective_numerators':c.tolist(),'denominator':D,'numerical_margin':res.x[-1]})
 (OUT/'exact_exposure_certificates.json').write_text(json.dumps(cert,indent=2))
 print('EXACT EXPOSURE SUCCESS',len(E),'minmargin',minmargin,'time',time.time()-start,flush=True)
 return cert
