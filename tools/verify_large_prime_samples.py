"""Verify four large-prime conic-billiard sample parameters.

Uses the original affine mirror reflection and the generalized Jacobi lift.
Point-count certificates:
  * 127-bit CM D=-19; trace sign checked against group-law points.
  * Curve25519: 2-isogeny of standard Montgomery curve (RFC 7748).
  * secp256k1-field CM D=-7; trace sign checked with a point.
  * NIST P-256-field CM D=-67; trace sign checked with a point.

These are *not* secure implementations. No claim of collision resistance,
constant-time behavior, or safe production key agreement is made.
"""
from __future__ import annotations
from math import isqrt, gcd
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from conic_billiard.core import Conic, inv
from conic_billiard.jacobi import ellipse_lift, project_ellipse

SAMPLES = [
    dict(label="Mersenne-127", prime=2**127 - 1, cm_D=19, r=8507059173023461587867128875462120321, cofactor=20,
         N=170141183460469231757342577509242406420,
         t=55579003044651936960779927009516534349,
         m=151326783926863011065916301464956629435),
    dict(label="Curve25519-field", prime=2**255 - 19, cm_D=None,
         r=2**252 + 27742317777372353535851937790883648493, cofactor=8,
         N=57896044618658097711785492504343953926856930875039260848015607506283634007912,
         t=4604094445982880512620096533988940857801848594535175639138443919056596161201,
         m=20800338683988658368647408995589388737092878452977063003340006470870624536394),
    dict(label="secp256k1-field", prime=2**256 - 2**32 - 977, cm_D=7, r=None, cofactor=None,
         N=115792089237316195423570985008687907852589613375810908478214554728305722413224,
         t=45713552117797890483722919253250213011584153888536609223224039558556197264329,
         m=91621925755784148632510482811437537317404925924607100258843383050386844822452),
    dict(label="NIST-P256-field", prime=2**256 - 2**224 + 2**192 + 2**96 - 1, cm_D=67, r=None, cofactor=None,
         N=115792089210356248762697446949407573530765394402118123589611579317542265087204,
         t=4928864831036071905185367580936228903682628659365973224999965418864697784127,
         m=95246762070185159368959048539282756182645115378201331967935434245525146677337),
]

def add_weierstrass(P,Q,p,a2,a4):
    """Affine addition on y^2=x^3+a2*x^2+a4*x (for verification only)."""
    if P is None:return Q
    if Q is None:return P
    x,y=P;u,v=Q
    if x==u and (y+v)%p==0:return None
    if x==u:
        slope=(3*x*x+2*a2*x+a4)*inv(2*y,p)%p
    else:
        slope=(v-y)*inv((u-x)%p,p)%p
    xx=(slope*slope-a2-x-u)%p
    yy=(slope*(x-xx)-y)%p
    return xx,yy

def mul_weierstrass(n,P,p,a2,a4):
    R=None
    while n:
        if n&1:R=add_weierstrass(R,P,p,a2,a4)
        P=add_weierstrass(P,P,p,a2,a4)
        n>>=1
    return R

def cm_sign_witness(p,m,N,D):
    # Any E in this CM isogeny class has trace ±u, where u²+Dv²=4p.
    u=abs(p+1-N)
    vv=(4*p-u*u)
    assert vv>0 and vv%D==0
    v=isqrt(vv//D)
    assert v*v*D+u*u==4*p
    alternate=2*(p+1)-N
    a2=-(1+m)%p;a4=m%p
    for x in range(2,300):
        yy=x*(x*x+a2*x+a4)%p
        if yy==0:continue
        if pow(yy,(p-1)//2,p)!=1:continue
        assert p%4==3
        y=pow(yy,(p+1)//4,p)
        P=(x,y)
        if mul_weierstrass(N,P,p,a2,a4) is None and mul_weierstrass(alternate,P,p,a2,a4) is not None:
            return P, u, v
    raise AssertionError('CM trace sign not distinguished by tested points')

def is_prime_64(n):
    if n<2:return False
    for p in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n%p==0:return n==p
    s=0;d=n-1
    while d%2==0:d//=2;s+=1
    for a in (2,325,9375,28178,450775,9780504,1795265022):
        if a%n==0:continue
        x=pow(a,d,n)
        if x in (1,n-1):continue
        for _ in range(s-1):
            x=x*x%n
            if x==n-1:break
        else:return False
    return True

def certify_127_subgroup_prime(r):
    factors={2:7,3:4,5:1,107:1,179:1,6703:1,20789:1,10463:1,5876498209681471:1}
    assert __import__('math').prod(q**e for q,e in factors.items())==r-1
    for q in factors:
        assert is_prime_64(q)
        for a in range(2,1000):
            if pow(a,r-1,r)==1 and gcd(pow(a,(r-1)//q,r)-1,r)==1:
                break
        else:raise AssertionError('no Lucas primality certificate witness')
    return True

def main():
    certify_127_subgroup_prime(SAMPLES[0]['r'])
    a,b=5,3
    for s in SAMPLES:
        p,N,t,m=s['prime'],s['N'],s['t'],s['m']
        C=Conic.ellipse(p,a,b)
        B=(a*(1-t*t)*inv(1+t*t,p)%p, 2*b*t*inv(1+t*t,p)%p)
        assert C.contains(B)
        lam,J,G=ellipse_lift(C,a,b,(a,0),B)
        assert J.m==m
        assert J.mul(N,G)==J.identity()
        X,Y=(a,0),B
        for n in range(2,65):
            Z=C.next_point(X,Y)
            assert Z==project_ellipse(J.mul(n,G),a,b,p)
            X,Y=Y,Z
        if s['cm_D'] is not None:
            witness,u,v=cm_sign_witness(p,m,N,s['cm_D'])
            assert u*u+s['cm_D']*v*v==4*p
        else:
            # Two-isogenous to Curve25519; verify the Legendre parameter.
            assert m==(486660*inv(486664,p))%p
            assert pow(486664,(p-1)//2,p)==1
        if s['r'] is not None:
            H=J.mul(s['cofactor'],G)
            assert H!=J.identity()
            assert J.mul(s['r'],H)==J.identity()
            status='verified prime-order translation subgroup'
        else:
            status='exact CM group order; subgroup prime factor unverified'
        print(f"{s['label']}: field_bits={p.bit_length()} N={N} status={status}")
        print(f'  B={B} lambda={lam} m={m} first 64 impacts OK')
    print('ALL FOUR SAMPLES VERIFIED')

if __name__=='__main__':main()
