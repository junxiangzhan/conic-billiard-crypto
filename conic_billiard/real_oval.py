"""Real specular billiards inside the bounded oval of y^2 = x^3 - x.

This is a Euclidean, not a finite-field, billiard. The oval is the
connected component with -1 <= x <= 0; the other component has x >= 1.

For incoming A -> B, the reflected ray B+t*w (t>0) first meets the oval
at the unique positive-parameter oval intersection. The other intersection
of the reflected LINE with the cubic belongs to the unbounded component
(or lies at the projective point at infinity for a vertical line).

Calculations use double precision; near tangencies require higher precision.
"""
from __future__ import annotations
import math
from dataclasses import dataclass
from typing import Optional

Point = tuple[float, float]


@dataclass(frozen=True)
class Hit:
    point: Point
    t: float
    other: Optional[Point]
    other_t: Optional[float]
    direction: Point
    roots: tuple[float, ...]


def oval(theta: float) -> Point:
    """One-to-one parametrization for theta in [0,pi), except endpoints."""
    q = math.sin(theta)
    r = math.cos(theta)
    return (-q*q, q*r*math.sqrt(1+q*q))


def cubic_residual(point: Point) -> float:
    x,y=point
    return y*y-x*x*x+x


def roots_for_ray(A: Point, B: Point, tol: float=1e-12):
    bx,by=B
    vx,vy=(bx-A[0],by-A[1])
    # Gradient of y^2 - x^3 + x
    nx,ny=(1-3*bx*bx,2*by)
    norm=nx*nx+ny*ny
    proj=2*(nx*vx+ny*vy)/norm
    wx,wy=(vx-proj*nx,vy-proj*ny)
    aa=-wx**3
    bb=wy*wy-3*bx*wx*wx
    cc=2*by*wy-(3*bx*bx-1)*wx
    # Coefficients can be tiny for nearly vertical rays. This relative
    # cutoff is fine for ordinary samples, but not a symbolic treatment.
    scale=max(abs(bb),abs(cc),1e-30)
    if abs(aa) < tol*scale:
        if abs(bb) < tol*max(abs(cc),1e-30):
            raise ArithmeticError('degenerate line / tangency')
        roots=[-cc/bb]
    else:
        discriminant=bb*bb-4*aa*cc
        if discriminant < -tol*max(bb*bb,abs(4*aa*cc),1e-30):
            raise ArithmeticError('no real remaining intersections')
        discriminant=max(discriminant,0.0)
        root_disc=math.sqrt(discriminant)
        q=-0.5*(bb+math.copysign(root_disc,bb))
        if q==0.0:
            roots=[-bb/(2*aa)]
        else:
            roots=[q/aa, cc/q]
    return (wx,wy), sorted(roots)


def next_hit(A: Point, B: Point) -> Hit:
    w,roots=roots_for_ray(A,B)
    all_hits=[(t,(B[0]+t*w[0],B[1]+t*w[1])) for t in roots]
    # Strict convexity and the inward reflection guarantee one oval hit.
    # Use a membership test (a second, independent selection method) to
    # verify that it is the least positive intersection parameter.
    positive=[(t,C) for t,C in all_hits if t>1e-10]
    if not positive:
        raise ArithmeticError('no forward hit')
    t,C=min(positive, key=lambda z:z[0])
    if not (-1-1e-7 <= C[0] <= 1e-7):
        raise ArithmeticError('first forward hit is not on the oval')
    if abs(cubic_residual(C))>1e-6:
        raise ArithmeticError('computed hit misses cubic')
    others=[(s,D) for s,D in all_hits if abs(s-t)>1e-8*max(1,abs(t))]
    other_t,other=others[0] if others else (None,None)
    if other is not None and other[0] < 1-1e-7:
        raise ArithmeticError('extra intersection not on unbounded component')
    return Hit(C,t,other,other_t,w,tuple(roots))


def forward_chain(A: Point, B: Point, steps: int) -> list[Point]:
    chain=[A,B]
    for _ in range(steps):
        chain.append(next_hit(chain[-2],chain[-1]).point)
    return chain


def distance(A: Point,B: Point) -> float:
    return math.hypot(A[0]-B[0],A[1]-B[1])
