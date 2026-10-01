import itertools
import numpy as np

def term(c,S,k):
 if S: c[S-1]+=k

def cmi(i,j,K):
 c=np.zeros(15);term(c,K|1<<i,1);term(c,K|1<<j,1);term(c,K,-1);term(c,K|1<<i|1<<j,-1);return c

def shannon_ingleton():
 A=[]
 for i in range(4):
  c=np.zeros(15);term(c,15,1);term(c,15^(1<<i),-1);A.append(c)
 for i,j in itertools.combinations(range(4),2):
  rest=15^((1<<i)|(1<<j));K=rest
  while True:
   A.append(cmi(i,j,K))
   if K==0:break
   K=(K-1)&rest
 for i,j in itertools.combinations(range(4),2):
  k,l=[x for x in range(4) if x not in (i,j)]
  A.append(cmi(i,j,1<<k)+cmi(i,j,1<<l)+cmi(k,l,0)-cmi(i,j,0))
 return np.array(A)

