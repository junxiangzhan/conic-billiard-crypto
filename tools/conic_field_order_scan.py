#!/usr/bin/env python3
"""Finite-field conic mirror billiards: associated curve order scan.

Research/educational only; not a secure DH implementation.

Conic: (x/a)^2 + epsilon*(y/b)^2 = 1, epsilon=+1 or -1.
Choose A=(a,0), B by rational slope t, and recover the caustic parameter.
The lifted orbit is a twisted Jacobi intersection:
    c^2 + epsilon*s^2 = 1, d^2 + m*s^2 = 1
with m=(epsilon*b*b-a*a)/(b*b-epsilon*lambda).
Its Jacobian is isogenous to y^2=x(x-epsilon)(x-m), so has the same
finite-field group cardinality (for m not 0 or epsilon).
"""
from __future__ import annotations
from collections import Counter
from math import isqrt
from random import Random
from sympy import factorint, isprime
import csv, json, argparse


def chitable(p):
    qr=bytearray(p)
    for y in range(1,(p+1)//2): qr[(y*y)%p]=1
    return qr


def card_legendre(p, eps, m, qr):
    s=0
    for x in range(p):
        f=(x*(x-eps)*(x-m))%p
        if f: s += 1 if qr[f] else -1
    return p+1+s


def base_point(p,t,eps,a=5,b=3):
    den=(1+eps*t*t)%p
    if den==0: return None
    i=pow(den,-1,p)
    return (a*(1-eps*t*t)*i%p,2*b*t*i%p)


def lam_m(p,t,eps,a=5,b=3):
    b2=b*b; a2=a*a
    den=(b2+a2*t*t)%p
    if den==0: return None
    lam=eps*a2*b2*t*t*pow(den,-1,p)%p
    den2=(b2-eps*lam)%p
    if den2==0: return None
    m=(eps*b2-a2)*pow(den2,-1,p)%p
    if m in (0,eps%p):return None
    return lam,m


def bounce(p,P,B,eps=1,a=5,b=3):
    q1=pow(a*a,-1,p); q2=eps*pow(b*b,-1,p)%p
    vx=(B[0]-P[0])%p;vy=(B[1]-P[1])%p
    nx=q1*B[0]%p; ny=q2*B[1]%p
    norm=(nx*nx+ny*ny)%p
    if norm==0:return None
    h=2*(vx*nx+vy*ny)*pow(norm,-1,p)%p
    wx=(vx-h*nx)%p;wy=(vy-h*ny)%p
    denom=(q1*wx*wx+q2*wy*wy)%p
    if denom==0:return None
    step=-2*(q1*B[0]*wx+q2*B[1]*wy)*pow(denom,-1,p)%p
    if step==0:return None
    return ((B[0]+step*wx)%p,(B[1]+step*wy)%p)


def jadd(p,eps,m,U,V):
    s,c,d=U; t,e,f=V
    D=(1-eps*m*s*s*t*t)%p
    if D==0:raise ZeroDivisionError
    inv=pow(D,-1,p)
    return ((s*e*f+t*c*d)*inv%p,
            (c*e-eps*s*t*d*f)*inv%p,
            (d*f-m*s*t*c*e)*inv%p)


def jmul(p,eps,m,n,U):
    R=(0,1,1)
    while n:
        if n&1:R=jadd(p,eps,m,R,U)
        n>>=1
        if n:U=jadd(p,eps,m,U,U)
    return R


def initial_step(p,t,eps,m,a=5,b=3):
    A=(a%p,0); B=base_point(p,t,eps,a,b)
    if B is None: return None
    C=bounce(p,A,B,eps,a,b)
    if C is None:return None
    ib=pow(b,-1,p);ia=pow(a,-1,p)
    s=B[1]*ib%p;c=B[0]*ia%p;s2=C[1]*ib%p
    denominator=(2*s*c)%p
    if denominator==0:return None
    d=(s2*(1-eps*m*pow(s,4,p))*pow(denominator,-1,p))%p
    if (c*c+eps*s*s-1)%p or (d*d+m*s*s-1)%p:return None
    G=(s,c,d)
    try:
        U2=jadd(p,eps,m,G,G)
    except ZeroDivisionError:return None
    if (a*U2[1]%p,b*U2[0]%p)!=C:return None
    return G


def gen_order(p,eps,m,G,N):
    assert jmul(p,eps,m,N,G)==(0,1,1)
    r=N
    for prime,exp in factorint(N).items():
        for _ in range(exp):
            try:
                if jmul(p,eps,m,r//int(prime),G)==(0,1,1):r//=int(prime)
                else:break
            except ZeroDivisionError:return None
    return r


def scan(p, eps, attempts=200, exhaustive=False, seed=20261118):
    assert isprime(p) and p>5 and p not in (2,3,5)
    qr=chitable(p)
    if exhaustive or p<1500:
        ts=range(1,p)
    else:
        r=Random(seed+17*p+eps)
        chosen=set()
        while len(chosen)<min(attempts,p-1):chosen.add(r.randrange(1,p))
        ts=sorted(chosen)
    card_cache={};records=[]; invalid=Counter()
    for t in ts:
        B=base_point(p,t,eps)
        if B is None:invalid['param_pole']+=1;continue
        res=lam_m(p,t,eps)
        if res is None:invalid['singular_caustic']+=1;continue
        lam,m=res
        if m not in card_cache:card_cache[m]=card_legendre(p,eps,m,qr)
        N=card_cache[m]
        facts={int(k):int(v) for k,v in factorint(N).items()}
        big=max(facts)
        cof=N//big
        rec=dict(p=p,kind='ellipse' if eps==1 else 'hyperbola',t=int(t),
                 Bx=B[0],By=B[1],lam=int(lam),m=int(m),
                 group_order=int(N),largest_prime=int(big),cofactor=int(cof),
                 factorization=' * '.join(f'{k}^{v}' if v>1 else str(k) for k,v in sorted(facts.items())))
        records.append(rec)
    records.sort(key=lambda x:(-x['largest_prime'],-x['group_order'],x['t']))
    for rec in records[:min(len(records),12)]:
        G=initial_step(p,rec['t'],eps,rec['m'])
        if G is None:rec['step_order']='not verified (affine exception)';continue
        try:
            ordG=gen_order(p,eps,rec['m'],G,rec['group_order'])
            rec['step_order']=int(ordG) if ordG is not None else 'exception'
            if ordG is not None:
                rec['step_prime_part']=int(rec['largest_prime'] if ordG%rec['largest_prime']==0 else 1)
        except (ZeroDivisionError, AssertionError):rec['step_order']='exception'
    result={'p':int(p),'kind':'ellipse' if eps==1 else 'hyperbola',
            'tested_initial_directions':len(records),'unique_curves':len(card_cache),
            'invalid':dict(invalid),'full_group_order_divisible_by_4':all(rec['group_order']%4==0 for rec in records),
            'best':records[:10]}
    return result,records


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--primes',nargs='+',type=int,default=[103,127,257,1009,1019,65521,65537,131071])
    parser.add_argument('--attempts',type=int,default=100)
    parser.add_argument('--json',default='')
    parser.add_argument('--csv',default='')
    args=parser.parse_args()
    summaries=[];allrows=[]
    for p in args.primes:
        for eps in (1,-1):
            summary,records=scan(p,eps,attempts=args.attempts)
            summaries.append(summary);allrows.extend(records)
            if not records:continue
            b=records[0]
            print(f"p={p}, {summary['kind']:9}, sampled={len(records):5}, curves={summary['unique_curves']:5}, "+
                  f"best N={b['group_order']}={b['factorization']}, r={b['largest_prime']}, cofactor={b['cofactor']}, "+
                  f"t={b['t']}, lambda={b['lam']}, m={b['m']}, step-order={b.get('step_order')}",flush=True)
    if args.json:
        with open(args.json,'w') as f:json.dump({'settings':{'a':5,'b':3},'summaries':summaries},f,ensure_ascii=False,indent=2)
    if args.csv:
        with open(args.csv,'w',newline='') as f:
            dw=csv.DictWriter(f,fieldnames=['p','kind','t','Bx','By','lam','m','group_order','largest_prime','cofactor','factorization','step_order','step_prime_part'],extrasaction='ignore')
            dw.writeheader();dw.writerows(allrows)

if __name__=='__main__':main()
