"""Poncelet-pencil experiment: "add a state" without the Jacobi addition law.

Research script, small primes only.  The algebraic
identities behind it are proved in tools/poncelet_proof.py; this script checks
them numerically on actual orbits and exercises the edge cases.

Setting (eps = +1 ellipse, eps = -1 hyperbola), over F_p:
  boundary  E : x^2/A + y^2/B = 1  with A = a^2, B = eps*b^2,  start P_0 = (a, 0)
  caustic   C : x^2/(A-lam) + y^2/(B-lam) = 1, fixed by the chord P_0 P_1
  pencil    C_mu : E + mu*C   (mu = 'inf' is C itself)
  lift      J : c^2 + eps*s^2 = 1, d^2 + m*s^2 = 1,  m = eps*(B-A)/(B-lam),
            pi(s,c,d) = (a*c, b*s), and P_k = pi([k]G) where G is the point over
            P_1 with mu(G) = inf.

Checks:
  [0] the lift: P_k = pi([k]G) on the whole orbit (for eps=+1 also: G equals
      conic_billiard.jacobi.ellipse_lift)
  [1] mu = inf (the caustic) reproduces reflection
  [2] chain P_0, P_b, P_2b, ... with the correct pencil member equals the orbit
  [3] the two pencil members tangent to chord P_0 P_b correspond to (s,c,d) and (s,c,-d)
  [4] rational (square-root-free) ladder vs sequential reflection; exceptions
      must lie inside the predicted set of degenerate cases
  [5] conic type: A*B*q0*q1 of C_mu(W) is x(W) - eps up to squares, and
      W -> x(W) - eps is a homomorphism J(F_p) -> F_p^*/F_p^*2

    python tools/poncelet_pencil.py
    python tools/poncelet_pencil.py --p 10007 --B 9 3304
    python tools/poncelet_pencil.py --kind hyperbola --p 1019
"""
from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from conic_billiard.compression import sqrt_mod  # noqa: E402
from conic_billiard.core import Conic, UndefinedReflection, inv  # noqa: E402

INF = 'inf'
Point = tuple[int, int]
JPoint = tuple[int, int, int]


def legendre(x: int, p: int) -> int:
    x %= p
    return 0 if x == 0 else (1 if pow(x, (p - 1) // 2, p) == 1 else -1)


class JacobiEps:
    """c^2 + eps*s^2 = 1, d^2 + m*s^2 = 1 with the Jacobi addition law (points (s,c,d))."""

    def __init__(self, p: int, eps: int, m: int):
        self.p, self.eps, self.m = p, eps, m % p
        self.O: JPoint = (0, 1, 1)
        self.T: JPoint = (0, 1, p - 1)

    def on(self, U: JPoint) -> bool:
        s, c, d = U
        p = self.p
        return (c * c + self.eps * s * s - 1) % p == 0 and (d * d + self.m * s * s - 1) % p == 0

    def add(self, U: JPoint, V: JPoint) -> JPoint:
        p, e, m = self.p, self.eps, self.m
        (s1, c1, d1), (s2, c2, d2) = U, V
        den = (1 - e * m * s1 * s1 * s2 * s2) % p
        if den == 0:
            raise ZeroDivisionError('exceptional affine Jacobi denominator')
        z = inv(den, p)
        return ((s1 * c2 * d2 + s2 * c1 * d1) * z % p,
                (c1 * c2 - e * s1 * s2 * d1 * d2) * z % p,
                (d1 * d2 - m * s1 * s2 * c1 * c2) * z % p)

    def neg(self, U: JPoint) -> JPoint:
        return (-U[0] % self.p, U[1], U[2])

    def mul(self, k: int, U: JPoint) -> JPoint:
        if k < 0:
            return self.mul(-k, self.neg(U))
        R = self.O
        while k:
            if k & 1:
                R = self.add(R, U)
            U = self.add(U, U)
            k >>= 1
        return R


class Pencil:
    """Diagonal conics q0 x^2 + q1 y^2 + q2 = 0 spanned by E and C (mu = INF is C)."""

    def __init__(self, p: int, A: int, B: int, lam: int):
        self.p = p
        self.e = (inv(A, p), inv(B, p), p - 1)
        self.c = (inv((A - lam) % p, p), inv((B - lam) % p, p), p - 1)

    def member(self, mu) -> tuple[int, int, int]:
        if mu == INF:
            return self.c
        return tuple((self.e[i] + mu * self.c[i]) % self.p for i in range(3))

    def tangency(self, X: Point, Y: Point) -> tuple[int, int, int]:
        """(A2, B1, C0): line XY is tangent to member mu iff A2 mu^2 + B1 mu + C0 = 0."""
        p, e, c = self.p, self.e, self.c
        l = ((X[1] - Y[1]) % p, (Y[0] - X[0]) % p, (X[0] * Y[1] - Y[0] * X[1]) % p)
        A2 = B1 = C0 = 0
        for i in range(3):
            j, k = [t for t in range(3) if t != i]
            C0 += l[i] ** 2 * e[j] * e[k]
            B1 += l[i] ** 2 * (e[j] * c[k] + c[j] * e[k])
            A2 += l[i] ** 2 * c[j] * c[k]
        return A2 % p, B1 % p, C0 % p

    def tangent_members(self, X: Point, Y: Point) -> list:
        """Both members tangent to chord XY (uses a square root)."""
        p = self.p
        A2, B1, C0 = self.tangency(X, Y)
        if A2 == 0:
            return [INF] + ([(-C0) * inv(B1, p) % p] if B1 else [])
        r = sqrt_mod(B1 * B1 - 4 * A2 * C0, p)
        if r is None:
            return []
        z = inv(2 * A2, p)
        return sorted({(-B1 + r) * z % p, (-B1 - r) * z % p}, key=str)

    def is_tangent(self, mu, X: Point, Y: Point) -> bool:
        A2, B1, C0 = self.tangency(X, Y)
        return A2 == 0 if mu == INF else (A2 * mu * mu + B1 * mu + C0) % self.p == 0

    def common_member(self, X1: Point, Y1: Point, X2: Point, Y2: Point):
        """The member tangent to BOTH chords, by elimination -- no square root.
        Raises if the two tangency quadratics are proportional (both roots shared)."""
        p = self.p
        a1, b1, c1 = self.tangency(X1, Y1)
        a2, b2, c2 = self.tangency(X2, Y2)
        num, den = (a1 * c2 - a2 * c1) % p, (a2 * b1 - a1 * b2) % p       # mu = num/den
        if den:
            return num * inv(den, p) % p
        num2, den2 = (c1 * b2 - c2 * b1) % p, (c2 * a1 - c1 * a2) % p     # mu = num2/den2
        if num2 and den2 == 0:
            return INF
        if den2:
            return num2 * inv(den2, p) % p
        raise UndefinedReflection('pencil_roots_coincide')

    def step(self, mu, prev: Point, X: Point) -> Point:
        """From X: the tangent to member mu that is NOT the line X-prev, second hit with E."""
        p, q, e = self.p, self.member(mu), self.e
        if q[0] * q[1] * q[2] % p == 0:
            raise UndefinedReflection('degenerate_pencil_member')

        def form(dv):
            xqd = (q[0] * X[0] * dv[0] + q[1] * X[1] * dv[1]) % p
            xqx = (q[0] * X[0] ** 2 + q[1] * X[1] ** 2 + q[2]) % p
            dqd = (q[0] * dv[0] ** 2 + q[1] * dv[1] ** 2) % p
            return (xqd * xqd - xqx * dqd) % p

        Af, Cf = form((1, 0)), form((0, 1))
        Bf = (form((1, 1)) - Af - Cf) * inv(2, p) % p
        x1, y1 = (prev[0] - X[0]) % p, (prev[1] - X[1]) % p
        if x1 == 0 and y1 == 0:
            raise UndefinedReflection('coincident_input_points')
        # Vieta on the binary quadratic Af u^2 + 2Bf uv + Cf v^2: other root of (x1:y1).
        dv = ((-2 * Bf) % p, Af) if x1 == 0 else ((Cf, (-2 * Bf) % p) if y1 == 0 else (Cf * y1 % p, Af * x1 % p))
        den = (e[0] * dv[0] ** 2 + e[1] * dv[1] ** 2) % p
        if den == 0:
            raise UndefinedReflection('next_intersection_at_infinity')
        tt = -2 * (e[0] * X[0] * dv[0] + e[1] * X[1] * dv[1]) * inv(den, p) % p
        if tt == 0:
            raise UndefinedReflection('repeated_tangent_intersection')
        return (X[0] + tt * dv[0]) % p, (X[1] + tt * dv[1]) % p


class Setup:
    """Orbit, caustic, pencil and lift for one start (a,0) -> B."""

    def __init__(self, kind: str, p: int, a: int, b: int, B: Point):
        self.kind, self.p, self.a, self.b = kind, p, a, b
        self.eps = 1 if kind == 'ellipse' else -1
        eps = self.eps
        self.E = Conic.ellipse(p, a, b) if eps == 1 else Conic.hyperbola(p, a, b)
        self.A_, self.B_ = a * a % p, eps * b * b % p
        P0 = (a % p, 0)
        if not self.E.contains(B) or B == P0:
            raise ValueError('B must be another point of the boundary')
        u, v, w = B[1] % p, (P0[0] - B[0]) % p, P0[0] * B[1] % p
        self.lam = (self.A_ * u * u + self.B_ * v * v - w * w) * inv(u * u + v * v, p) % p
        self.m = eps * (self.B_ - self.A_) * inv((self.B_ - self.lam) % p, p) % p
        if self.m in (0, eps % p):
            raise ValueError('singular state curve (m = 0 or m = eps)')
        self.J = JacobiEps(p, eps, self.m)
        self.pen = Pencil(p, self.A_, self.B_, self.lam)
        self.K = (a * a + b * b * (self.m - eps)) % p
        if self.K == 0:
            raise ValueError('K = a^2 + b^2 (m - eps) = 0: every pencil member is degenerate')
        s, c = B[1] * inv(b, p) % p, B[0] * inv(a, p) % p
        # d_G solves D(G) = a^2 m (c-1) + eps K (1+d) = 0  (mu(G) = infinity)
        self.G = (s, c, (-(a * a * self.m * (c - 1)) * inv(eps * self.K % p, p) - 1) % p)
        P = [P0, B]
        while not (len(P) > 2 and P[-2] == P0 and P[-1] == B):
            P.append(self.E.next_point(P[-2], P[-1]))      # may raise UndefinedReflection
            if len(P) > 4 * p + 10:
                raise RuntimeError('orbit did not close')
        self.P = P[:-2]
        self.N = len(self.P)
        self.P_back = self.E.next_point(B, P0)

    def pi(self, U: JPoint) -> Point:
        return self.a * U[1] % self.p, self.b * U[0] % self.p

    def mu_of(self, U: JPoint):
        """Closed form mu(W) = -(A-B)(m-eps)(c-1) / (a^2 m (c-1) + eps K (1+d))."""
        p, eps, a = self.p, self.eps, self.a
        s, c, d = U
        num = -(self.A_ - self.B_) * (self.m - eps) * (c - 1) % p
        den = (a * a * self.m * (c - 1) + eps * self.K * (1 + d)) % p
        if den == 0:
            return INF if num else None
        return num * inv(den, p) % p

    def ladder(self, n: int) -> Point:
        """P_n from (P_1, P_2) with anchors P_0, P_1, P_{-1}; square-root free."""
        pen, P, Pb = self.pen, self.P, self.P_back
        R0, R1 = P[1], P[2 % self.N]
        for bit in bin(n)[3:]:
            mu_m = pen.common_member(P[0], R0, P[1], R1)        # +m
            if bit == '0':
                R0, R1 = pen.step(mu_m, P[0], R0), pen.step(mu_m, P[1], R1)
            else:
                mu_m1 = pen.common_member(P[0], R1, Pb, R0)     # +(m+1)
                R0, R1 = pen.step(mu_m, P[1], R1), pen.step(mu_m1, P[0], R1)
        return R0

    def ladder_hypotheses_fail(self, n: int) -> bool:
        """True if some hypothesis of Theorems A/C fails for a step of ladder(n)
        (docs/geometric_group_law.md Theorem D; conditions in docs/elliptic_correspondence.md section 2).  Sufficient for correctness when
        False; when True the ladder may raise or may still return the right point."""
        J, G = self.J, self.G
        two_torsion = {J.O, J.T, (0, self.p - 1, 1), (0, self.p - 1, self.p - 1)}

        def bad_member(W: JPoint) -> bool:          # mu(W) undefined or C_mu(W) degenerate
            mu = self.mu_of(W)
            return W in two_torsion or mu is None or (mu != INF and (mu + 1) % self.p == 0)

        same = lambda U, V: U == V or U == J.add(J.T, J.neg(V))       # pi(U) == pi(V)

        def bad_common(W: JPoint, ch1, ch2) -> bool:
            """Chords pi(U1)pi(V1), pi(U2)pi(V2), both of difference W: a chord collapses,
            or their other roots mu(T-(U+V)) coincide (fibre lemma: equal up to sign)."""
            (U1, V1), (U2, V2) = ch1, ch2
            if same(U1, V1) or same(U2, V2) or bad_member(W):
                return True
            o1, o2 = J.add(J.T, J.neg(J.add(U1, V1))), J.add(J.T, J.neg(J.add(U2, V2)))
            return o1 == o2 or o1 == J.neg(o2)

        def bad_step(W: JPoint, X: JPoint) -> bool:  # step X-W -> X -> X+W
            Y, Z = J.add(X, J.neg(W)), J.add(X, W)
            q, (x, y) = self.pen.member(self.mu_of(W)), self.pi(X)
            on_member = (q[0] * x * x + q[1] * y * y + q[2]) % self.p == 0
            return same(Y, Z) or same(X, Y) or same(X, Z) or on_member

        try:
            m = 1
            for bit in bin(n)[3:]:
                Wm, Wm1 = J.mul(m, G), J.mul(m + 1, G)
                if bad_common(Wm, (Wm, J.O), (Wm1, G)):                     # chords P0Pm, P1Pm+1
                    return True
                if bit == '0':
                    if bad_step(Wm, Wm) or bad_step(Wm, Wm1):
                        return True
                    m = 2 * m
                else:
                    if bad_common(Wm1, (Wm1, J.O), (Wm, J.neg(G))) \
                            or bad_step(Wm, Wm1) or bad_step(Wm1, Wm1):    # chords P0Pm+1, P-1Pm
                        return True
                    m = 2 * m + 1
        except ZeroDivisionError:
            return True
        return False


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--kind', choices=['ellipse', 'hyperbola'], default='ellipse')
    ap.add_argument('--p', type=int, default=1009)
    ap.add_argument('--a', type=int, default=5)
    ap.add_argument('--b', type=int, default=3)
    ap.add_argument('--B', type=int, nargs=2, metavar=('X', 'Y'),
                    help='second impact; default: first point with a closed, defined orbit')
    ap.add_argument('--seed', type=int, default=1)
    args = ap.parse_args()
    p = args.p
    rng = random.Random(args.seed)

    def try_setup(B):
        try:
            return Setup(args.kind, p, args.a, args.b, B)
        except (UndefinedReflection, ValueError, ZeroDivisionError, RuntimeError):
            return None

    if args.B:
        S = try_setup((args.B[0] % p, args.B[1] % p))
        if S is None:
            print('this start gives no closed, fully defined orbit with a non-singular lift; pick another --B or --p')
            return 2
    else:
        S = None
        conic = Conic.ellipse(p, args.a, args.b) if args.kind == 'ellipse' else Conic.hyperbola(p, args.a, args.b)
        for x in range(2, p):
            yy = (1 - conic.qx * x * x) * inv(conic.qy, p) % p
            y = sqrt_mod(yy, p)
            if y:
                S = try_setup((x, y))
                if S and S.N > 20:
                    break
        if S is None:
            print('no usable start found')
            return 2
    J, P, N, pen, eps = S.J, S.P, S.N, S.pen, S.eps
    print(f'{args.kind} p={p} a={args.a} b={args.b}  start={P[0]},{P[1]}  lambda={S.lam}  m={S.m}  orbit period={N}')
    wrong = 0

    # [0] the lift
    W, ok = J.O, True
    for k in range(N):
        ok &= J.on(W) and S.pi(W) == P[k]
        W = J.add(W, S.G)
    ok &= W == J.O
    if eps == 1:
        from conic_billiard.jacobi import ellipse_lift
        _, _, G2 = ellipse_lift(S.E, args.a, args.b, P[0], P[1])
        ok &= G2 == S.G
    print(f'[0] lift: P_k = pi([k]G) for all k, [N]G = O' + (', G = ellipse_lift(...)' if eps == 1 else '') + f': {ok}')
    wrong += not ok

    # [1] caustic = reflection
    ok = all(pen.step(INF, P[k - 1], P[k]) == P[(k + 1) % N] for k in range(1, N))
    print(f'[1] mu=inf step reproduces reflection on the whole orbit: {ok}')
    wrong += not ok

    # [2] chains with the closed-form member
    good = bad = exc = 0
    for bb in range(2, min(N - 1, 40)):
        mu = S.mu_of(J.mul(bb, S.G))
        ch = [P[0], P[bb]]
        try:
            for _ in range(10):
                ch.append(pen.step(mu, ch[-2], ch[-1]))
        except (UndefinedReflection, ZeroDivisionError):
            exc += 1
            continue
        if all(pen.is_tangent(mu, P[k], P[(k + bb) % N]) for k in range(N)) and \
                all(ch[k] == P[k * bb % N] for k in range(len(ch))):
            good += 1
        else:
            bad += 1
    print(f'[2] closed-form mu([b]G): all chords P_k P_k+b tangent, chain = orbit: {good} ok, {bad} wrong, {exc} exceptional')
    wrong += bad

    # [3] the two roots <-> (s,c,d), (s,c,-d)
    good = bad = exc = 0
    for bb in range(2, min(N - 1, 40)):
        Q = J.mul(bb, S.G)
        Qf = (Q[0], Q[1], (-Q[2]) % p)
        roots = pen.tangent_members(P[0], P[bb])
        if set(map(str, roots)) != {str(S.mu_of(Q)), str(S.mu_of(Qf))}:      # a double root when Q = Qf
            bad += 1
            continue
        for aa in range(2, 16):
            try:
                outs = set()
                for mu in roots:
                    ch = [P[0], P[bb]]
                    for _ in range(aa - 1):
                        ch.append(pen.step(mu, ch[-2], ch[-1]))
                    outs.add(ch[aa])
            except (UndefinedReflection, ZeroDivisionError):
                exc += 1
                continue
            good += outs == {S.pi(J.mul(aa, Q)), S.pi(J.mul(aa, Qf))}
            bad += outs != {S.pi(J.mul(aa, Q)), S.pi(J.mul(aa, Qf))}
    print(f'[3] roots = mu(s,c,d), mu(s,c,-d); chains = pi([a]Q), pi([a]Q\'): {good} ok, {bad} wrong, {exc} exceptional')
    wrong += bad

    # [4] rational ladder; exceptions must be predicted
    ns = list(range(2, 4 * N)) + [rng.randrange(2, 2 ** 64) for _ in range(300)]
    good = bad = 0
    exc, unpredicted = [], []
    for n in ns:
        try:
            r = S.ladder(n)
        except (UndefinedReflection, ZeroDivisionError):
            exc.append(n)
            if not S.ladder_hypotheses_fail(n):
                unpredicted.append(n)
            continue
        if r == P[n % N]:
            good += 1
        else:
            bad += 1
            if not S.ladder_hypotheses_fail(n):
                unpredicted.append(n)
    predicted = sum(S.ladder_hypotheses_fail(n) for n in ns)
    small = [n for n in exc if n < 4 * N]
    degenerate_below_N = [n for n in range(2, N) if S.ladder_hypotheses_fail(n)]
    big_exc = len([n for n in exc if n >= 4 * N])
    print(f'[4] rational ladder on n=2..{4 * N - 1} and 300 random 64-bit n: {good} ok, {bad} wrong, '
          f'{len(exc)} exceptional ({len(small)} below 4N, {big_exc} of the random ones), '
          f'predicted degenerate: {predicted}, unpredicted failures: {len(unpredicted)}')
    print(f'    exceptional n below 4N: {small[:24]}{" ..." if len(small) > 24 else ""}')
    print(f'    predicted degenerate n in [2, N): {degenerate_below_N}')
    wrong += bad + len(unpredicted)

    # [5] conic type and the 2-descent homomorphism delta(W) = x(W) - eps
    def delta(U: JPoint) -> int:
        s, c, d = U
        if s == 0:      # 2-torsion: standard convention via the other two roots
            e1, e2, e3 = 0, eps % p, S.m
            x0 = {J.O: None, J.T: e3, (0, p - 1, 1): e2, (0, p - 1, p - 1): e1}[U]
            if x0 is None:
                return 1
            if x0 == eps % p:
                return legendre((e2 - e1) * (e2 - e3), p)
            return legendre(x0 - eps, p)
        return legendre((1 + c) * (1 + d) * inv(s * s, p) - eps, p)

    def q0q1_class(mu) -> int:
        q = pen.member(mu)
        return legendre(S.A_ * S.B_ * q[0] * q[1], p)

    agree = total = 0
    for k in range(1, N):
        W = J.mul(k, S.G)
        Wf = (W[0], W[1], (-W[2]) % p)
        mu, muf = S.mu_of(W), S.mu_of(Wf)
        if mu in (None, INF) or muf in (None, INF) or W[0] == 0:
            continue
        total += 1
        agree += q0q1_class(mu) == delta(W) and q0q1_class(muf) == delta(W) * legendre(S.m - eps, p)
    pts = []
    while len(pts) < 60:
        s = rng.randrange(1, p)
        c = sqrt_mod((1 - eps * s * s) % p, p)
        d = sqrt_mod((1 - S.m * s * s) % p, p)
        if c is not None and d is not None:
            pts.append((s, c if rng.random() < .5 else -c % p, d if rng.random() < .5 else -d % p))
    assoc = assoc_n = 0
    for U, V, W3 in zip(pts, pts[1:], pts[2:]):
        try:
            assoc += J.add(J.add(U, V), W3) == J.add(U, J.add(V, W3))
            assoc_n += 1
        except ZeroDivisionError:
            pass
    print(f'[5] addition law is associative on random points of J(F_p): {assoc}/{assoc_n}')
    wrong += assoc_n - assoc
    hom = hom_n = 0
    for U in pts:
        for V in pts[:15]:
            try:
                Wuv = J.add(U, V)
            except ZeroDivisionError:
                continue
            hom_n += 1
            hom += delta(Wuv) == delta(U) * delta(V)
    print(f'[5] type: A*B*q0*q1 ~ delta(W) and delta(W)*(m-eps) for the two roots: {agree}/{total}; '
          f'delta is a homomorphism on random points: {hom}/{hom_n}')
    print(f'    delta(G) = {delta(S.G):+d}, (m-eps|p) = {legendre(S.m - eps, p):+d}, (-1|p) = {legendre(-1, p):+d}')
    wrong += (total - agree) + (hom_n - hom)

    print('RESULT:', 'no wrong values' if wrong == 0 else f'{wrong} WRONG')
    return 1 if wrong else 0


if __name__ == '__main__':
    raise SystemExit(main())
