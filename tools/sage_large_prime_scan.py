#!/usr/bin/env sage -python
"""SageMath search for suitable large-prime conic billiard groups.

Run using: sage -python sage_large_prime_scan.py --trials 30
Uses E_{eps,m}: y^2=x(x-eps)(x-m), 2-isogenous to the lifted
conic-billiard Jacobi intersection, so both have the same group size.
This is a RESEARCH screening tool, not a secure key-generation program.
"""
from __future__ import annotations
import argparse
from sage.all import ZZ, GF, EllipticCurve, set_random_seed

PRIMES={
    'mersenne127': 2**127-1,
    'curve25519': 2**255-19,
    'secp256k1': 2**256-2**32-977,
    'nistp256': 2**256-2**224+2**192+2**96-1,
    'curve448': 2**448-2**224-1,
}

def run(p, eps, a, b, trials):
    p=ZZ(p)
    assert p.is_prime() and p>5
    F=GF(p)
    aa=F(a);bb=F(b);ee=F(eps)
    if aa==0 or bb==0 or ee*bb**2==aa**2:
        raise ValueError('singular or invalid conic parameters')
    print(f'p={p} bits={p.nbits()} eps={eps} trials={trials}',flush=True)
    hits=0
    for i in range(trials):
        t=F.random_element()
        if t==0 or 1+ee*t*t==0:continue
        den=bb**2+aa**2*t**2
        if den==0:continue
        lam=ee*aa**2*bb**2*t**2/den
        if bb**2-ee*lam==0:continue
        m=(ee*bb**2-aa**2)/(bb**2-ee*lam)
        if m in (0,ee):continue
        E=EllipticCurve(F,[F(0),-(ee+m),F(0),ee*m,F(0)])
        N=ZZ(E.cardinality())  # Sage uses optimized point-counting algorithms.
        assert abs(N-(p+1)) <= 2*p.isqrt()+2
        assert N%4==0
        r=N//4
        if r.is_prime():
            # Select a point in the r-primary part. This demonstrates that
            # the group has an r-subgroup. It is not the billiard step itself.
            Q=E.random_point()
            H=4*Q
            if H==E(0) or r*H!=E(0):continue
            hits+=1
            print('FOUND N=4r,', 'i=',i+1,'lambda=',ZZ(lam),'m=',ZZ(m),
                  'N=',N,'r=',r,'B=',(ZZ(aa*(1-ee*t*t)/(1+ee*t*t)),
                                  ZZ(bb*2*t/(1+ee*t*t))),flush=True)
            break
        elif (i+1)%5==0:
            print('checked',i+1,'curve candidates; no N=4r yet',flush=True)
    print('found=',hits,'(lack of a hit does not prove none exist)',flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--prime',choices=list(PRIMES),default='curve25519')
    ap.add_argument('--trials',type=int,default=20)
    ap.add_argument('--a',type=int,default=5)
    ap.add_argument('--b',type=int,default=3)
    ap.add_argument('--kind',choices=['ellipse','hyperbola','both'],default='both')
    args=ap.parse_args()
    set_random_seed(20261118)
    for eps in ([1,-1] if args.kind=='both' else [1 if args.kind=='ellipse' else -1]):
        run(PRIMES[args.prime],eps,args.a,args.b,args.trials)
