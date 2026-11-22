"""Algebraic reflection on affine conics over an odd prime field.

NOT real-angle geometry: this is the rational algebraic extension of the
usual Euclidean reflection formula with dot product x1*y1 + x2*y2.
All arithmetic is in F_p. Undefined transitions are NOT silently extended.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import TypeAlias

Point: TypeAlias = tuple[int, int]
State: TypeAlias = tuple[Point, Point]
Matrix4: TypeAlias = tuple[tuple[int, int, int, int], ...]


class UndefinedReflection(ValueError):
    """The rational affine reflection map is undefined at this state."""

    def __init__(self, reason: str):
        super().__init__(reason)
        self.reason = reason


@dataclass(frozen=True)
class Conic:
    """F(x,y) = qx*x² + qy*y² + lx*x + ly*y + c = 0 over F_p.

    Supports diagonal Q only, which is sufficient for axis-aligned ellipses,
    hyperbolas, and the parabola y=x². The normal is Q*B + ell/2.
    """

    p: int
    qx: int
    qy: int
    lx: int = 0
    ly: int = 0
    c: int = -1
    label: str = 'conic'

    def __post_init__(self):
        if self.p < 3 or self.p % 2 == 0 or not is_probable_prime(self.p):
            raise ValueError('p must be an odd (probable) prime')
        for attr in ('qx', 'qy', 'lx', 'ly', 'c'):
            object.__setattr__(self, attr, getattr(self, attr) % self.p)

    @classmethod
    def ellipse(cls, p: int, a: int = 5, b: int = 3) -> 'Conic':
        if a % p == 0 or b % p == 0:
            raise ValueError('nonzero semiaxes required')
        return cls(p, inv(a * a, p), inv(b * b, p), label='ellipse')

    @classmethod
    def hyperbola(cls, p: int, a: int = 5, b: int = 3) -> 'Conic':
        """(x/a)^2 - (y/b)^2 = 1 over F_p, with p odd."""
        if a % p == 0 or b % p == 0:
            raise ValueError('nonzero semiaxes required')
        return cls(p, inv(a * a, p), -inv(b * b, p), label='hyperbola')

    @classmethod
    def circle(cls, p: int, radius: int = 1) -> 'Conic':
        if radius % p == 0:
            raise ValueError('nonzero radius required')
        return cls(p, inv(radius*radius,p), inv(radius*radius,p), label='circle')

    @classmethod
    def parabola(cls, p: int) -> 'Conic':
        """y = x², or x² - y = 0."""
        return cls(p, 1, 0, 0, -1, 0, 'parabola')

    def value(self, P: Point) -> int:
        x, y = P
        return (self.qx*x*x + self.qy*y*y + self.lx*x + self.ly*y + self.c) % self.p

    def contains(self, P: Point) -> bool:
        return self.value(P) == 0

    def boundary_points(self) -> list[Point]:
        """Exhaustive affine points; O(p^2), toy primes only."""
        p = self.p
        return [(x, y) for x in range(p) for y in range(p) if self.contains((x, y))]

    def normal(self, B: Point) -> Point:
        p = self.p
        inv2 = inv(2, p)
        return ((self.qx * B[0] + self.lx * inv2) % p,
                (self.qy * B[1] + self.ly * inv2) % p)

    def _data(self, A: Point, B: Point):
        """Compute H, reflected direction w, and line parameter tau."""
        p = self.p
        if not self.contains(A) or not self.contains(B):
            raise ValueError('A and B must both lie on the conic')
        if A == B:
            raise UndefinedReflection('coincident_input_points')
        v = ((B[0] - A[0]) % p, (B[1] - A[1]) % p)
        nx, ny = self.normal(B)
        norm = (nx * nx + ny * ny) % p
        if norm == 0:
            raise UndefinedReflection('isotropic_normal')
        scale = 2 * inv(norm, p) % p
        h00 = (1 - scale * nx * nx) % p
        h01 = (-scale * nx * ny) % p
        h11 = (1 - scale * ny * ny) % p
        H = ((h00, h01), (h01, h11))
        w0 = (h00 * v[0] + h01 * v[1]) % p
        w1 = (h01 * v[0] + h11 * v[1]) % p
        quad = (self.qx * w0 * w0 + self.qy * w1 * w1) % p
        if quad == 0:
            raise UndefinedReflection('next_intersection_at_infinity')
        derivative = (2 * (self.qx * B[0] * w0 + self.qy * B[1] * w1)
                      + self.lx * w0 + self.ly * w1) % p
        tau = -derivative * inv(quad, p) % p
        if tau == 0:
            raise UndefinedReflection('repeated_tangent_intersection')
        return H, (w0, w1), tau

    def next_point(self, A: Point, B: Point) -> Point:
        """Return C for ordered pair (A,B) -> (B,C)."""
        _, (w0, w1), tau = self._data(A, B)
        p = self.p
        C = ((B[0] + tau * w0) % p, (B[1] + tau * w1) % p)
        if C == B or not self.contains(C):
            raise AssertionError('internal reflection/second-intersection error')
        return C

    def transfer_matrix(self, A: Point, B: Point) -> Matrix4:
        """State-dependent 4x4 M taking [Bx,By,Ax,Ay] -> [Cx,Cy,Bx,By].

        M = [[I+tau H, -tau H], [I, 0]], where tau and H depend on A,B.
        This is NOT a fixed matrix for repeated squaring.
        """
        H, _, tau = self._data(A, B)
        p = self.p
        h00, h01 = H[0]
        h11 = H[1][1]
        return (
            ((1 + tau*h00) % p, tau*h01 % p, -tau*h00 % p, -tau*h01 % p),
            (tau*h01 % p, (1 + tau*h11) % p, -tau*h01 % p, -tau*h11 % p),
            (1, 0, 0, 0),
            (0, 1, 0, 0),
        )

    def step(self, state: State) -> State:
        A, B = state
        return B, self.next_point(A, B)

    def trace(self, A: Point, B: Point, steps: int) -> State:
        """Return state after 'steps' transitions from (A,B). O(steps)."""
        if steps < 0:
            raise ValueError('steps must be nonnegative')
        state: State = (A, B)
        for _ in range(steps):
            state = self.step(state)
        return state


def inv(x: int, p: int) -> int:
    if x % p == 0:
        raise ZeroDivisionError('zero has no inverse in F_p')
    return pow(x, -1, p)


def is_probable_prime(n: int) -> bool:
    """Deterministic for n < 2^64; probable-prime test above that.

    NOT a certified prime generator or cryptographic parameter validation.
    Large fields support individual reflections, never exhaustive scans.
    """
    if n < 2:
        return False
    small = (2,3,5,7,11,13,17,19,23,29,31,37)
    if n in small:
        return True
    if any(n % q == 0 for q in small):
        return False
    s = 0
    d = n - 1
    while d % 2 == 0:
        d >>= 1
        s += 1
    if n < 2**64:
        bases = (2, 325, 9375, 28178, 450775, 9780504, 1795265022)
    else:
        bases = (2,3,5,7,11,13,17,19,23,29,31,37,41,43,47)
    for a in bases:
        a %= n
        if a == 0:
            continue
        x = pow(a, d, n)
        if x in (1, n-1):
            continue
        for _ in range(s-1):
            x = x*x % n
            if x == n-1:
                break
        else:
            return False
    return True


def matvec4(M: Matrix4, state: State, p: int) -> State:
    A, B = state
    # state is (previous, current) but column vector stacks (current, previous)
    vec = (B[0], B[1], A[0], A[1])
    out = [sum(M[i][j] * vec[j] for j in range(4)) % p for i in range(4)]
    return (out[2], out[3]), (out[0], out[1])


def determinant4(M: Matrix4, p: int) -> int:
    """Determinant via modular Gaussian elimination."""
    R = [list(row) for row in M]
    value = 1
    for j in range(4):
        k = next((i for i in range(j, 4) if R[i][j] % p), None)
        if k is None:
            return 0
        if k != j:
            R[j], R[k] = R[k], R[j]
            value = -value
        pivot = R[j][j] % p
        value = (value * pivot) % p
        pi = inv(pivot, p)
        for i in range(j+1, 4):
            scale = R[i][j] * pi % p
            for t in range(j, 4):
                R[i][t] = (R[i][t] - scale * R[j][t]) % p
    return value % p
