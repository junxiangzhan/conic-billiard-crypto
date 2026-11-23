"""A narrow elliptical-billiard/Jacobi fast-jump *reference implementation*.

This is a separate ECC-type reference baseline. It does NOT derive a fast
formula directly from the state-dependent 4x4 reflection matrix, and it does
NOT claim to support a hyperbolic boundary. Affine formulas have exceptions.
"""
from __future__ import annotations
from dataclasses import dataclass
from .core import Conic, Point, inv

JacobiPoint = tuple[int, int, int]  # (s,c,d)


@dataclass(frozen=True)
class JacobiGroup:
    p: int
    m: int

    def identity(self) -> JacobiPoint:
        return 0, 1, 1

    def on_curve(self, U: JacobiPoint) -> bool:
        s, c, d = U
        return (s*s+c*c-1) % self.p == 0 and (d*d+self.m*s*s-1) % self.p == 0

    def add(self, U: JacobiPoint, V: JacobiPoint) -> JacobiPoint:
        """sn/cn/dn algebraic addition in the displayed affine chart."""
        if not self.on_curve(U) or not self.on_curve(V):
            raise ValueError('invalid Jacobi point')
        p, m = self.p, self.m
        s, c, d = U
        t, e, f = V
        denominator = (1-m*s*s*t*t) % p
        if denominator == 0:
            raise ZeroDivisionError('exceptional affine Jacobi denominator')
        z = inv(denominator, p)
        return ((s*e*f+t*c*d)*z % p,
                (c*e-s*t*d*f)*z % p,
                (d*f-m*s*t*c*e)*z % p)

    def mul(self, scalar: int, U: JacobiPoint) -> JacobiPoint:
        """Binary double-and-add; O(log scalar) field group operations."""
        if scalar < 0:
            s, c, d = U
            return self.mul(-scalar, (-s % self.p, c, d))
        R = self.identity()
        while scalar:
            if scalar & 1:
                R = self.add(R, U)
            U = self.add(U, U)
            scalar >>= 1
        return R


def ellipse_lift(curve: Conic, a: int, b: int,
                 A: Point, B: Point) -> tuple[int, JacobiGroup, JacobiPoint]:
    """Choose a Jacobi lift for A=(a,0), from B and its successor.

    Valid only on this explicit nondegenerate chart. Raises on exceptional
    values rather than guessing a lift or extending the rational map.
    """
    p = curve.p
    if curve.label != 'ellipse' or curve.qx != inv(a*a, p) or curve.qy != inv(b*b, p):
        raise ValueError('expected matching axis-aligned ellipse')
    if A != (a % p, 0) or not curve.contains(B) or A == B:
        raise ValueError('expected A=(a,0), B on the same ellipse and B != A')
    u, v = (B[1]-A[1]) % p, (A[0]-B[0]) % p
    w = (A[0]*B[1]-A[1]*B[0]) % p
    lam = (a*a*u*u+b*b*v*v-w*w) * inv(u*u+v*v, p) % p
    m = (b*b-a*a)*inv(b*b-lam, p) % p
    if m in (0, 1):
        raise ValueError('singular Jacobi model')
    s, c = B[1] * inv(b,p) % p, B[0] * inv(a,p) % p
    C = curve.next_point(A,B)
    s2 = C[1] * inv(b,p) % p
    d = s2*(1-m*pow(s,4,p))*inv(2*s*c,p) % p
    group = JacobiGroup(p,m)
    G = s,c,d
    if not group.on_curve(G):
        raise AssertionError('bad Jacobi lift')
    if project_ellipse(group.mul(2,G),a,b,p) != C:
        raise AssertionError('doubling failed to match reflected impact')
    return lam, group, G


def project_ellipse(U: JacobiPoint, a: int, b: int, p: int) -> Point:
    s, c, _ = U
    return a*c % p, b*s % p


def bsgs(group: JacobiGroup, base: JacobiPoint,
         target: JacobiPoint, order: int) -> int | None:
    """Simple baby-step/giant-step attack, small groups only."""
    from math import isqrt
    m = isqrt(order)+1
    baby = {}
    R = group.identity()
    for i in range(m):
        baby.setdefault(R, i)
        R = group.add(R,base)
    mG = group.mul(m,base)
    neg_mG = (-mG[0]%group.p,mG[1],mG[2])
    Q = target
    for j in range(m+1):
        if Q in baby:
            n = (j*m + baby[Q]) % order
            if group.mul(n,base) == target:
                return n
        Q = group.add(Q, neg_mG)
    return None
