"""Timing: square-root-free pencil ladder vs Jacobi double-and-add.

Both compute the impact point after n bounces of the ellipse billiard
x^2/25 + y^2/9 = 1 started at (5, 0).  Every timed result is also checked for
equality, so on the 127- and 255-bit primes this doubles as a correctness test
far beyond the reach of orbit enumeration.  Pure Python, variable time.

    python tools/poncelet_bench.py
"""
from __future__ import annotations

import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from conic_billiard.compression import sqrt_mod  # noqa: E402
from conic_billiard.core import Conic, UndefinedReflection, inv  # noqa: E402
from conic_billiard.jacobi import ellipse_lift, project_ellipse  # noqa: E402
from poncelet_pencil import Pencil  # noqa: E402

PRIMES = {'p = 10007': 10007, 'p = 2^127 - 1': 2 ** 127 - 1, 'p = 2^255 - 19': 2 ** 255 - 19}


def setup(p: int, rng: random.Random, a: int = 5, b: int = 3):
    E = Conic.ellipse(p, a, b)
    A = (a, 0)
    while True:
        x = rng.randrange(2, p)
        y = sqrt_mod((1 - E.qx * x * x) * inv(E.qy, p) % p, p)
        if not y:
            continue
        try:
            B = (x, y)
            lam, group, G = ellipse_lift(E, a, b, A, B)
            P2 = E.next_point(A, B)
            Pb = E.next_point(B, A)
            pen = Pencil(p, a * a, b * b, lam)
            return E, group, G, pen, (A, B, P2, Pb)
        except (UndefinedReflection, ValueError, ZeroDivisionError, AssertionError):
            continue


def ladder(pen: Pencil, anchors, n: int):
    P0, P1, P2, Pb = anchors
    R0, R1 = P1, P2
    for bit in bin(n)[3:]:
        mu_m = pen.common_member(P0, R0, P1, R1)
        if bit == '0':
            R0, R1 = pen.step(mu_m, P0, R0), pen.step(mu_m, P1, R1)
        else:
            mu_m1 = pen.common_member(P0, R1, Pb, R0)
            R0, R1 = pen.step(mu_m, P1, R1), pen.step(mu_m1, P0, R1)
    return R0


def main() -> int:
    rng = random.Random(7)
    wrong = 0
    print(f'{"field":<16}{"bits of n":>10}{"Jacobi ms":>12}{"pencil ms":>12}{"ratio":>8}  agree')
    for label, p in PRIMES.items():
        E, group, G, pen, anchors = setup(p, rng)
        ns = [rng.randrange(2 ** (p.bit_length() - 2), 2 ** p.bit_length()) for _ in range(20)]
        tj = tp = 0.0
        agree = 0
        for n in ns:
            t0 = time.perf_counter()
            want = project_ellipse(group.mul(n, G), 5, 3, p)
            t1 = time.perf_counter()
            got = ladder(pen, anchors, n)
            t2 = time.perf_counter()
            tj += t1 - t0
            tp += t2 - t1
            agree += got == want
        wrong += len(ns) - agree
        print(f'{label:<16}{p.bit_length():>10}{1000 * tj / len(ns):>12.3f}{1000 * tp / len(ns):>12.3f}'
              f'{tp / tj:>8.1f}  {agree}/{len(ns)}')
    print('RESULT:', 'all equal' if not wrong else f'{wrong} DIFFER')
    return 1 if wrong else 0


if __name__ == '__main__':
    raise SystemExit(main())
