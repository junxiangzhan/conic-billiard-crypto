"""Numerical checks for docs/security_dictionary.md (billiard <-> elliptic-curve conditions).

    python tools/dictionary_check.py

Every section compares a proved statement with the independent step-by-step
reflection of conic_billiard/core.py (or with an exact point count):

  [1] circle: C = 2(A.B)/(B.B) B - A, and every orbit period divides p - (-1/p)
  [2] orbits through a focus pair: chords alternate between the two foci and the
      period divides 2(p-1) or 2(p+1), as the fixed points on the focal axis are
      rational or not
  [3] parabola y = x^2: x_{k+1} = h x_k - x_{k-1} with a constant h, and the period
      divides p - (k(1+k)/p)
  [4] rectangular hyperbola a = b: m is not 0 or eps, and the lift equals reflection
  [5] point count: #J = 2*outside + on and #twist = 2*inside + on, over all p+1
      projective impacts (impacts at infinity included)
  [6] axis and centre mirrors act on an orbit as translation by the 2-torsion points

Small primes only.  Exits non-zero on any mismatch.
"""
from __future__ import annotations

import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from conic_billiard.compression import sqrt_mod  # noqa: E402
from conic_billiard.core import Conic, UndefinedReflection, inv  # noqa: E402
from poncelet_pencil import Setup, legendre  # noqa: E402

BAD: list[str] = []
rng = random.Random(1)


def report(name: str, ok: int, total: int, extra: str = '') -> None:
    flag = 'ok' if ok == total else 'FAIL'
    print(f'  [{flag}] {name}: {ok}/{total} {extra}')
    if ok != total:
        BAD.append(name)


def period(E: Conic, A, B, limit: int):
    """Period of the state (A, B) under reflection, or None if a dead end is hit."""
    state, n = (A, B), 0
    try:
        while True:
            state, n = E.step(state), n + 1
            if state == (A, B):
                return n
            if n > limit:
                return None
    except UndefinedReflection:
        return None


def second_hit(E: Conic, X, F):
    """Second intersection of the line X F with E (X on E)."""
    p = E.p
    dx, dy = (F[0] - X[0]) % p, (F[1] - X[1]) % p
    den = (E.qx * dx * dx + E.qy * dy * dy) % p
    if den == 0:
        return None
    t = -2 * (E.qx * X[0] * dx + E.qy * X[1] * dy) * inv(den, p) % p
    return None if t == 0 else ((X[0] + t * dx) % p, (X[1] + t * dy) % p)


def collinear(P, Q, R, p) -> bool:
    return ((Q[0] - P[0]) * (R[1] - P[1]) - (Q[1] - P[1]) * (R[0] - P[0])) % p == 0


def weierstrass_order(p: int, eps: int, m: int) -> int:
    return 1 + sum(1 + legendre(x * (x - eps) * (x - m), p) for x in range(p))


def circle() -> None:
    print('[1] circle x^2 + y^2 = 25')
    for p in (103, 1009, 1019, 2017):
        E = Conic.circle(p, 5)
        pts = E.boundary_points()
        order = p - legendre(-1, p)
        ok_step = ok_per = total = 0
        for _ in range(200):
            A, B = rng.sample(pts, 2)
            dot, bb = (A[0] * B[0] + A[1] * B[1]) % p, (B[0] * B[0] + B[1] * B[1]) % p
            C = tuple((2 * dot * inv(bb, p) * B[i] - A[i]) % p for i in range(2))
            total += 1
            ok_step += E.next_point(A, B) == C
            n = period(E, A, B, 2 * p + 4)
            ok_per += n is not None and order % n == 0
        report(f'p={p}: reflection = mirror in line OB', ok_step, total)
        report(f'p={p}: period divides p - (-1/p) = {order}', ok_per, total)


def focus_orbits() -> None:
    print('[2] orbits through a focus pair (a=5, b=3)')
    primes = (103, 1009, 1013, 1019, 1021, 2003, 2017, 2027, 2029)
    for kind, eps in (('ellipse', 1), ('hyperbola', -1)):
        for p in primes:
            E = Conic.ellipse(p, 5, 3) if eps == 1 else Conic.hyperbola(p, 5, 3)
            A_, B_ = 25 % p, eps * 9 % p
            pairs = []
            f = sqrt_mod((A_ - B_) % p, p)                      # lambda = B: foci on the x-axis
            if f:
                pairs.append(('x-axis foci (lambda=B)', ((f, 0), (p - f, 0)), 1))   # fixed points (+-a, 0)
            f = sqrt_mod((B_ - A_) % p, p)                      # lambda = A: foci on the y-axis
            if f:
                pairs.append(('y-axis foci (lambda=A)', ((0, f), (0, p - f)), legendre(B_, p)))
            pts = [X for X in E.boundary_points() if X[1] != 0 and X[0] != 0]
            for name, (F1, F2), split in pairs:
                bound = 2 * (p - split)
                ok_alt = ok_per = total = closed = 0
                for X in rng.sample(pts, min(60, len(pts))):
                    Y = second_hit(E, X, F1)
                    if Y is None or Y == X:
                        continue
                    seq, state, n = [X, Y], (X, Y), None
                    try:                                         # walk one period, or up to a dead end
                        while len(seq) < 4 * p + 8:
                            state = E.step(state)
                            if state == (X, Y):
                                n = len(seq) - 1
                                break
                            seq.append(state[1])
                    except UndefinedReflection:
                        pass
                    total += 1
                    ok_alt += all(collinear(seq[i], seq[i + 1], (F1, F2)[i % 2], p) for i in range(len(seq) - 1))
                    if n is not None:
                        closed += 1
                        ok_per += bound % n == 0
                if total:
                    report(f'{kind} p={p} {name}: every chord passes through the other focus', ok_alt, total)
                    report(f'{kind} p={p} {name}: period | {bound}', ok_per, closed,
                           f'({total - closed} orbits reach an affine dead end, not counted)')


def parabola() -> None:
    print('[3] parabola y = x^2')
    for p in (103, 1009, 1019, 2017):
        E = Conic.parabola(p)
        ok_rec = ok_per = total = closed = 0
        for _ in range(200):
            u, v = rng.sample(range(p), 2)
            if (1 + 4 * u * v) % p == 0:
                continue
            k = (u - v) ** 2 * inv(1 + 4 * u * v, p) % p
            if k * (1 + k) % p == 0:
                continue
            h = (2 + 4 * k) % p
            A, B = (u, u * u % p), (v, v * v % p)
            try:
                seq, state = [u, v], (A, B)
                for _ in range(6):
                    state = E.step(state)
                    seq.append(state[1][0])
            except UndefinedReflection:
                continue
            total += 1
            ok_rec += all(seq[i + 1] == (h * seq[i] - seq[i - 1]) % p for i in range(1, len(seq) - 1))
            n = period(E, A, B, 2 * p + 4)
            if n is not None:
                closed += 1
                ok_per += (p - legendre(k * (1 + k), p)) % n == 0
        report(f'p={p}: fixed linear recurrence', ok_rec, total)
        report(f'p={p}: period divides p - (k(1+k)/p)', ok_per, closed,
               f'({total - closed} orbits reach an affine dead end, not counted)')


def setups(kind: str, p: int, a: int, b: int, limit: int):
    E = Conic.ellipse(p, a, b) if kind == 'ellipse' else Conic.hyperbola(p, a, b)
    out = []
    for x in range(p):
        y = sqrt_mod((1 - E.qx * x * x) * inv(E.qy, p) % p, p)
        if not y or (x, y) == (a % p, 0):
            continue
        try:
            out.append(Setup(kind, p, a, b, (x, y)))
        except (ValueError, UndefinedReflection, ZeroDivisionError, RuntimeError):
            continue
        if len(out) >= limit:
            break
    return out


def rectangular() -> None:
    print('[4] rectangular hyperbola x^2/25 - y^2/25 = 1')
    for p in (1009, 2003, 10007):
        S_all = setups('hyperbola', p, 5, 5, 3)
        ok = sum(S.m not in (0, p - 1) and all(S.P[k] == S.pi(S.J.mul(k, S.G)) for k in range(S.N))
                 for S in S_all)
        report(f'p={p}: non-singular and lift = reflection', ok, len(S_all),
               f'(periods {[S.N for S in S_all]})')


def counting() -> None:
    print('[5] #J = 2*outside + on, #twist = 2*inside + on (all p+1 impacts)')
    for kind, eps in (('ellipse', 1), ('hyperbola', -1)):
        for p in (1009, 1019, 2017, 2027):
            E = Conic.ellipse(p, 5, 3) if eps == 1 else Conic.hyperbola(p, 5, 3)
            A_, B_ = 25 % p, eps * 9 % p
            affine = E.boundary_points()
            r = sqrt_mod((-B_ * inv(A_, p)) % p, p)             # impacts at infinity (1 : +-r : 0)
            infinite = [] if not r else ([(1, r)] if r == 0 else [(1, r), (1, p - r)])
            assert len(affine) + len(infinite) == p + 1
            ok = total = 0
            for P1 in rng.sample([X for X in affine if X != (5, 0)], 6):
                u, v, w = P1[1], (5 - P1[0]) % p, 5 * P1[1] % p
                if (u * u + v * v) % p == 0:
                    continue
                lam = (A_ * u * u + B_ * v * v - w * w) * inv((u * u + v * v) % p, p) % p
                if lam in (0, A_, B_):
                    continue
                m = eps * (B_ - A_) * inv((B_ - lam) % p, p) % p
                q0, q1 = inv((A_ - lam) % p, p), inv((B_ - lam) % p, p)
                cnt = {1: 0, 0: 0, -1: 0}
                for X in affine:                                 # tangent directions from X
                    Cx = (q0 * X[0] ** 2 + q1 * X[1] ** 2 - 1) % p
                    form = lambda d0, d1: ((q0 * X[0] * d0 + q1 * X[1] * d1) ** 2 - Cx * (q0 * d0 ** 2 + q1 * d1 ** 2)) % p
                    Af, Cf = form(1, 0), form(0, 1)
                    Bf = (form(1, 1) - Af - Cf) * inv(2, p) % p
                    cnt[legendre(Bf * Bf - Af * Cf, p)] += 1
                for (x, y) in infinite:                          # parallel lines n.X = k tangent to C
                    cnt[legendre((A_ - lam) * y * y + (B_ - lam) * x * x, p)] += 1
                N = weierstrass_order(p, eps, m)
                total += 1
                ok += 2 * cnt[1] + cnt[0] == N and 2 * cnt[-1] + cnt[0] == 2 * p + 2 - N
            report(f'{kind} p={p} ({len(infinite)} impacts at infinity)', ok, total)


def mirrors() -> None:
    print('[6] mirrors = translation by 2-torsion')
    for kind, p in (('ellipse', 1019), ('hyperbola', 2003), ('ellipse', 10007), ('hyperbola', 10007)):
        for S in setups(kind, p, 5, 3, 2):
            J = S.J
            table = {'x-axis': ((0, 1, p - 1), lambda P: (P[0], -P[1] % p)),
                     'centre': ((0, p - 1, 1), lambda P: (-P[0] % p, -P[1] % p)),
                     'y-axis': ((0, p - 1, p - 1), lambda P: (-P[0] % p, P[1]))}
            ok = sum(all(f(S.P[k]) == S.pi(J.add(J.mul(k, S.G), T)) for k in range(S.N))
                     for T, f in table.values())
            report(f'{kind} p={p} N={S.N}', ok, 3)


if __name__ == '__main__':
    circle()
    focus_orbits()
    parabola()
    rectangular()
    counting()
    mirrors()
    print('RESULT:', 'all checks pass' if not BAD else f'{len(BAD)} FAILED: {BAD}')
    sys.exit(1 if BAD else 0)
