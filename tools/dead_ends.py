"""Where the affine reflection formula stops, the lift keeps going.

For every second impact B of the conic x^2/25 + eps*y^2/9 = 1 (start (5, 0)),
this walks the translation orbit k*G on the Weierstrass model
y^2 = x(x-eps)(x-m) -- where only O is at infinity, so the walk never fails --
and maps each state back to the state curve J and to an impact point.  Then:

  * wherever core.Conic.next_point is defined, it must agree with the lift;
  * wherever it raises UndefinedReflection, we record what the lift does there
    (the next state is a point at infinity of the projective conic, or an
    ordinary point where the affine formula merely divides by zero).

The orbit always closes ([N]G = O), so every dead end is a removable gap of
the affine formula, not a real end of the dynamics.

It also checks the two turning-point rules that the Poncelet step (the other
tangent from the impact to the caustic, docs/model.md) predicts at a dead end
of the state (A, B):

  * impact at infinity: B, inf, -B, -A, ...  The two tangents to the caustic
    in the asymptotic direction are mirror images through the centre, so the
    ball re-enters at -B and the orbit is point-symmetric about that step.
  * isotropic normal: B, B, A, ...  The tangent of E at B is isotropic, hence
    tangent to every confocal conic, so the ball stays one step and then
    retraces its path.

Small primes only.

    python tools/dead_ends.py
    python tools/dead_ends.py --primes 10009      # hyperbola with isotropic normals
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from conic_billiard.compression import sqrt_mod  # noqa: E402
from conic_billiard.core import Conic, UndefinedReflection, inv  # noqa: E402

A_AXIS, B_AXIS = 5, 3


def lift(p: int, eps: int, B):
    """(m, G) for the orbit (a,0) -> B, as in tools/poncelet_pencil.py (D(G) = 0)."""
    a, b = A_AXIS, B_AXIS
    A_, B_ = a * a % p, eps * b * b % p
    u, v, w = B[1] % p, (a - B[0]) % p, a * B[1] % p
    lam = (A_ * u * u + B_ * v * v - w * w) * inv((u * u + v * v) % p, p) % p
    m = eps * (B_ - A_) * inv((B_ - lam) % p, p) % p
    K = (a * a + b * b * (m - eps)) % p
    s, c = B[1] * inv(b, p) % p, B[0] * inv(a, p) % p
    return m, (s, c, (-(a * a * m * (c - 1)) * inv(eps * K % p, p) - 1) % p)


class Weierstrass:
    """y^2 = x(x-eps)(x-m) with the chord-tangent law; None is O."""

    def __init__(self, p: int, eps: int, m: int):
        self.p, self.e, self.m = p, eps % p, m % p

    def from_J(self, U):
        p, (s, c, d) = self.p, U
        return ((1 + c) * (1 + d) * inv(s * s % p, p) % p,
                (1 + c) * (1 + d) * (c + d) * inv(pow(s, 3, p), p) % p)

    def to_J(self, P):
        """Back to (s, c, d); ('inf', alpha) when the J-point is at infinity (impact at infinity)."""
        p, e, m = self.p, self.e, self.m
        if P is None:
            return (0, 1, 1)
        x, y = P
        if y == 0:
            return {m: (0, 1, p - 1), e: (0, p - 1, 1), 0: (0, p - 1, p - 1)}[x]
        al, be = y * inv((x - m) % p, p) % p, y * inv((x - e) % p, p) % p     # (1+c)/s, (1+d)/s
        den = (al * al + e) % p
        if den == 0:
            return ('inf', al)
        s = 2 * al * inv(den, p) % p
        return (s, (al * s - 1) % p, (be * s - 1) % p)

    def add(self, P, Q):
        p = self.p
        a2, a4 = -(self.e + self.m) % p, self.e * self.m % p
        if P is None or Q is None:
            return Q if P is None else P
        (x1, y1), (x2, y2) = P, Q
        if x1 == x2 and (y1 + y2) % p == 0:
            return None
        lam = ((3 * x1 * x1 + 2 * a2 * x1 + a4) * inv(2 * y1 % p, p) if P == Q
               else (y2 - y1) * inv((x2 - x1) % p, p)) % p
        x3 = (lam * lam - a2 - x1 - x2) % p
        return x3, (lam * (x1 - x3) - y1) % p


def turning_point_ok(pts: list, k: int, reason: str, p: int) -> bool:
    """Dead end of the state (pts[k-1], pts[k]); None marks an impact at infinity.

    Checks the Poncelet turning-point rule on the whole period:
      next_intersection_at_infinity: pts[k+1] is at infinity and P_{k+1+j} = -P_{k+1-j};
      isotropic_normal:              pts[k+1] = pts[k]       and P_{k+1+j} = P_{k-j}.
    """
    N = len(pts)
    neg = lambda Q: None if Q is None else ((-Q[0]) % p, (-Q[1]) % p)
    if reason == 'next_intersection_at_infinity':
        return pts[(k + 1) % N] is None and all(
            pts[(k + 1 + j) % N] == neg(pts[(k + 1 - j) % N]) for j in range(N))
    if reason == 'isotropic_normal':
        return all(pts[(k + 1 + j) % N] == pts[(k - j) % N] for j in range(N))
    return True                     # coincident_input_points follows an isotropic step


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--primes', type=int, nargs='+', default=[103, 1009, 1019, 2003])
    args = ap.parse_args()
    wrong = 0
    for kind, eps in (('hyperbola', -1), ('ellipse', 1)):
        for p in args.primes:
            E = Conic.hyperbola(p, A_AXIS, B_AXIS) if eps == -1 else Conic.ellipse(p, A_AXIS, B_AXIS)
            st = Counter()
            for x in range(p):
                y = sqrt_mod((1 - E.qx * x * x) * inv(E.qy, p) % p, p)
                if not y or (x, y) == (A_AXIS, 0):
                    continue
                try:
                    m, G = lift(p, eps, (x, y))
                except ZeroDivisionError:
                    st['skipped: lift formula undefined for this start'] += 1
                    continue
                if m in (0, eps % p) or G[0] == 0:
                    st['skipped: singular state curve or 2-torsion step'] += 1
                    continue
                W = Weierstrass(p, eps, m)
                Gw, P, states = W.from_J(G), None, [(0, 1, 1)]
                while True:
                    P = W.add(P, Gw)
                    states.append(W.to_J(P))
                    if P is None:
                        break
                impact = lambda U: (A_AXIS * U[1] % p, B_AXIS * U[0] % p)
                finite = lambda U: U[0] != 'inf'
                pts = [impact(U) if finite(U) else None for U in states[:-1]]
                st['orbits'] += 1
                hit = False
                for k in range(1, len(states) - 1):
                    U0, U1, U2 = states[k - 1], states[k], states[k + 1]
                    if not (finite(U0) and finite(U1)):
                        st['steps not checked (an input impact is at infinity)'] += 1
                        continue
                    try:
                        nxt = E.next_point(impact(U0), impact(U1))
                        if finite(U2):
                            st['reflection agrees with lift'] += nxt == impact(U2)
                            wrong += nxt != impact(U2)
                        else:
                            st['WRONG: reflection defined but lift goes to infinity'] += 1
                            wrong += 1
                    except UndefinedReflection as ex:
                        hit = True
                        if not finite(U2):
                            what = 'through an impact at infinity'
                        elif impact(U2) == impact(U1):
                            what = 'with the same impact again (the ball stays one step)'
                        else:
                            what = 'to an ordinary impact'
                        st[f'dead end "{ex.reason}": lift continues {what}'] += 1
                        if turning_point_ok(pts, k, ex.reason, p):
                            if ex.reason != 'coincident_input_points':
                                st[f'  ... and follows the turning-point rule for "{ex.reason}"'] += 1
                        else:
                            st[f'WRONG: turning-point rule fails for "{ex.reason}"'] += 1
                            wrong += 1
                st['orbits with a dead end (all of them close)'] += hit
                if hit:
                    N = len(pts)
                    n_inf = pts.count(None)
                    n_iso = sum(1 for k in range(N) if pts[k] is not None and pts[k] == pts[(k + 1) % N])
                    st[f'  turning points per such orbit: {n_inf} at infinity + {n_iso} isotropic, '
                       f'period {"even" if N % 2 == 0 else "odd"}'] += 1
            print(f'{kind} p={p}:')
            for k, v in st.items():
                print(f'    {k}: {v}')
    print('RESULT:', 'reflection and lift agree wherever reflection is defined' if not wrong else f'{wrong} DISAGREE')
    return 1 if wrong else 0


if __name__ == '__main__':
    raise SystemExit(main())
