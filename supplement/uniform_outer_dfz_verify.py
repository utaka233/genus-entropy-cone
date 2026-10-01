"""Fraction-only verification of the genus-four bound E4 <= 7/440.
Relies only on known DFZ inequalities (50),(52) and the repair lemma.
"""
from uniform_outer_verify import main as verify_repair, shannon,ingleton,mi,add,dot
from fractions import Fraction as Q
from pathlib import Path
import json

T=[0]*14+[1]
coeffs50=(2,4,2,1,2,0,0,2,3)
coeffs52=(2,3,3,2,2,0,0,2,0)
def dfz(c,p):
 i,j,k,l=[1<<x for x in p]
 fs=[mi(i,j),mi(i,j,k),mi(i,k,j),mi(j,k,i),mi(i,j,l),mi(i,l,j),mi(j,l,i),mi(k,l),mi(k,l,i)]
 return add(*zip([-c[0]]+list(c[1:]),fs))

def main():
 verify_repair()
 S=shannon();J=ingleton()[0]
 perms=[(0,1,2,3),(1,0,2,3),(0,1,3,2),(1,0,3,2)]
 fields=[(14,S[k]) for k in range(4)]+[(42,mi(1,2,12)),(14,mi(4,8)),(10,mi(4,8,1)),(10,mi(4,8,2))]
 fields += [(3,dfz(coeffs50,p)) for p in perms]+[(1,dfz(coeffs52,p)) for p in perms]
 assert add(*fields)==add((88,J),(14,T))
 # Verify sparse stored Farkas identity with paper coefficient provenance.
 c=json.loads(Path(__file__).with_name('uniform_outer_dfz_certificate.json').read_text())
 total=[Q(0)]*15
 for u in c['used_inequalities']:
  provenance=u['provenance']
  if 'shannon_index' in provenance:a=S[provenance['shannon_index']]
  else:
   num=provenance['dfz_equation']
   coeff=coeffs50 if num==50 else coeffs52
   assert num in (50,52) and tuple(provenance['tuple'])==coeff
   a=dfz(coeff,provenance['permutation'])
  assert a==u['row']
  w=Q(u['weight']);assert w>=0
  total=[x+w*y for x,y in zip(total,a)]
 assert total==[Q(j)+Q(7,44)*t for j,t in zip(J,T)]
 assert c['uniform_error_bound']=='7/440'
 print('Exact DFZ (50),(52) identity gives J >= -7T/44.')
 print('Using the repair lemma and L4 subset C4, normalized entropy error <= 7/440.')
if __name__=='__main__':main()
