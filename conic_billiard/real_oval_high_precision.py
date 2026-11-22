"""Optional mpmath verification of real-oval reversibility at 100 digits.

Install separately: pip install mpmath
Run: python -m conic_billiard.real_oval_high_precision
"""
from __future__ import annotations


def run(digits: int=100, bounces: int=80):
    import mpmath as mp
    mp.mp.dps=digits

    def oval(t):
        q=mp.sin(t);r=mp.cos(t)
        return -q*q,q*r*mp.sqrt(1+q*q)

    def hit(A,B):
        bx,by=B
        vx,vy=bx-A[0],by-A[1]
        nx,ny=1-3*bx*bx,2*by
        scale=2*(nx*vx+ny*vy)/(nx*nx+ny*ny)
        wx,wy=vx-scale*nx,vy-scale*ny
        a=-wx**3
        b=wy**2-3*bx*wx**2
        c=2*by*wy-(3*bx**2-1)*wx
        eps=mp.power(10,-digits+20)
        if abs(a)<eps*max(1,abs(b),abs(c)):
            roots=[-c/b]
        else:
            disc=b*b-4*a*c
            if disc<0:
                raise ValueError('Negative discriminant')
            sd=mp.sqrt(disc)
            roots=[(-b+sd)/(2*a),(-b-sd)/(2*a)]
        t=min(t for t in roots if t>eps)
        P=bx+t*wx,by+t*wy
        assert P[0]>-1-eps and P[0]<eps
        return P

    A=oval(mp.mpf('0.8651862621171743'))
    B=oval(mp.mpf('2.3356020830597064'))
    seq=[A,B]
    for _ in range(bounces):
        seq.append(hit(seq[-2],seq[-1]))
    rev=[seq[-1],seq[-2]]
    for _ in range(bounces):
        rev.append(hit(rev[-2],rev[-1]))
    max_err=max(max(abs(P[0]-Q[0]),abs(P[1]-Q[1]))
                for P,Q in zip(rev,seq[::-1]))
    return {'precision_digits':digits,'bounces_each_direction':bounces,
            'max_reverse_error':mp.nstr(max_err,18),
            'final_point':[mp.nstr(z,25) for z in seq[-1]]}

if __name__=='__main__':
    import json
    print(json.dumps(run(),indent=2))
