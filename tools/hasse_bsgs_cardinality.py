#!/usr/bin/env python3
"""Certify group cardinalities of elliptic curves over 31..64 bit prime fields.
BSGS in the Hasse interval. Returns only if mathematically unique given
an explicitly certified point order (and 4|#E due to full rational 2-torsion).
"""
import math,random,time,sys
from sympy import isprime, factorint

O=None
class Curve:
 def __init__(self,p,eps,m):
  self.p=p;self.a2=(-(eps+m))%p;self.a4=(eps*m)%p
 def is_on(self,P):
  if P is None:return True
  x,y=P;p=self.p
  return (y*y-(x*x*x+self.a2*x*x+self.a4*x))%p==0
 def neg(self,P):return None if P is None else (P[0],(-P[1])%self.p)
 def add(self,P,Q):
  p=self.p
  if P is None:return Q
  if Q is None:return P
  x,y=P;z,w=Q
  if x==z:
   if (y+w)%p==0:return None
   if y!=w:raise ValueError('not equal, not opposite?')
   lam=((3*x*x+2*self.a2*x+self.a4)*pow(2*y%p,-1,p))%p
  else:lam=((w-y)*pow((z-x)%p,-1,p))%p
  x3=(lam*lam-self.a2-x-z)%p
  y3=(lam*(x-x3)-y)%p
  return (x3,y3)
 def mul(self,n,P):
  if n<0:return self.mul(-n,self.neg(P))
  R=None
  while n:
   if n&1:R=self.add(R,P)
   n//=2
   if n:P=self.add(P,P)
  return R
 def randompoint(self,rng):
  p=self.p
  for _ in range(500):
   x=rng.randrange(p)
   rhs=(x*x*x+self.a2*x*x+self.a4*x)%p
   if rhs==0:return (x,0)
   if pow(rhs,(p-1)//2,p)==1:
    y=mod_sqrt(rhs,p)
    return (x,y)
  raise RuntimeError('point not found')

def mod_sqrt(n,p):
 if n==0:return 0
 if p%4==3:return pow(n,(p+1)//4,p)
 # Tonelli-Shanks
 q=p-1;s=0
 while q%2==0:s+=1;q//=2
 z=2
 while pow(z,(p-1)//2,p)!=p-1:z+=1
 m=s;c=pow(z,q,p);t=pow(n,q,p);r=pow(n,(q+1)//2,p)
 while t!=1:
  i=1;tt=t*t%p
  while tt!=1:tt=tt*tt%p;i+=1
  b=pow(c,1<<(m-i-1),p)
  r=r*b%p;t=t*b*b%p;c=b*b%p;m=i
 return r

def bsgs_interval(E,Q,low,high):
 W=high-low;M=math.isqrt(W)+1
 tab={};cur=None
 for j in range(M):
  tab.setdefault(cur,j)
  cur=E.add(cur,Q)
 giant=E.neg(E.mul(low,Q))
 minus_m=E.neg(cur)
 for i in range(W//M+2):
  if giant in tab:
   n=low+i*M+tab[giant]
   if n<=high and E.mul(n,Q) is None:return n,len(tab)
  giant=E.add(giant,minus_m)
 raise RuntimeError('failed Hasse search')

def get_order(E,Q,candidate):
 assert E.mul(candidate,Q) is None
 N=candidate
 for ell,mult in factorint(N).items():
  for j in range(mult):
   if E.mul(N//ell,Q) is None:N//=ell
   else:break
 return N

def count_certified(p,eps,m,seed=42,maxtrials=8):
 assert isprime(p)
 E=Curve(p,eps,m)
 t=time.perf_counter()
 H=2*math.isqrt(p)+3
 lo=p+1-H;hi=p+1+H
 rng=random.Random(seed)
 orders=[];Q=None
 for i in range(maxtrials):
  Q=E.randompoint(rng)
  found,n= bsgs_interval(E,Q,lo,hi)
  r=get_order(E,Q,found)
  orders.append(r)
  # Curve has three nontrivial Fp rational 2-torsion points; group order multiple 4.
  # Other possible #E are multiples of lcm(4, point orders)
  L=math.lcm(4,*orders)
  start=((lo+L-1)//L)*L
  options=list(range(start,hi+1,L))
  if len(options)==1:
   N=options[0]
   assert all(E.mul(N,E.randompoint(rng)) is None for _ in range(12))
   return dict(N=N,r_cert=r,cardinality_candidates=options,certification_modulus=L,
               trials=i+1,secs=round(time.perf_counter()-t,3),Nfactor={str(k):int(v) for k,v in factorint(N).items()},point_orders=orders)
 raise RuntimeError(f'could not uniquely determine order: {options}')

if __name__=='__main__':
 from conic_field_order_scan import lam_m
 PRIMES=[(2**31-1,'Mersenne31'),(2**61-1,'Mersenne61'),(2**64-2**32+1,'Goldilocks64')]
 for p,name in PRIMES:
  for eps in (1,-1):
   t=42
   while not lam_m(p,t,eps):t+=1
   lam,m=lam_m(p,t,eps)
   print('START',name,p,eps,t,lam,m,flush=True)
   print(count_certified(p,eps,m),flush=True)
