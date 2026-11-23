"""From an elliptic curve back to a billiard (the converse of the lift).

Given p, eps and m, and a point W = (s, c, d) on the state curve
    J : c^2 + eps*s^2 = 1,  d^2 + m*s^2 = 1     (isomorphic to y^2 = x(x-eps)(x-m)),
take a = 1 and b^2 = -(m(c-1) + eps(1+d)) / (eps (m-eps)(1+d)).  When b^2 is a
non-zero square, the billiard in x^2/a^2 + eps*y^2/b^2 = 1 started at (a, 0)
with second impact (a c, b s) has state curve parameter m and step G = W
(proved in tools/poncelet_proof.py, lemma L11).  Then "count the bounces" on
that billiard is exactly the discrete logarithm to base W on J.

This script tries random (m, W) and checks the claim with the independent
numerical lift of tools/poncelet_pencil.py.  Small primes only.

    python tools/realize_curve.py
"""
from __future__ import annotations

import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from conic_billiard.compression import sqrt_mod  # noqa: E402
from conic_billiard.core import UndefinedReflection, inv  # noqa: E402
from poncelet_pencil import Setup  # noqa: E402


def main() -> int:
    rng = random.Random(3)
    wrong = 0
    print(f'{"p":>6} {"p%4":>4} {"kind":<10}{"points":>8}{"b^2 square":>12}{"m, G match":>12}'
          f'{"orbit has a dead end":>22}')
    for p in (1009, 1019, 2003, 2011):
        for eps in (1, -1):
            kind = 'ellipse' if eps == 1 else 'hyperbola'
            tried = square = match = dead = 0
            while tried < 200:
                m = rng.randrange(2, p)
                if m == eps % p:
                    continue
                s = rng.randrange(1, p)
                c = sqrt_mod((1 - eps * s * s) % p, p)
                d = sqrt_mod((1 - m * s * s) % p, p)
                if c is None or d is None:
                    continue
                c = c if rng.random() < .5 else -c % p
                d = d if rng.random() < .5 else -d % p
                if (1 + d) % p == 0:
                    continue
                tried += 1
                b = sqrt_mod((-(m * (c - 1) + eps * (1 + d))) * inv(eps * (m - eps) * (1 + d) % p, p) % p, p)
                if not b:
                    continue
                square += 1
                try:
                    S = Setup(kind, p, 1, b, (c, b * s % p))
                except (UndefinedReflection, ZeroDivisionError, RuntimeError):
                    dead += 1           # the orbit itself hits an affine dead end (see tools/dead_ends.py)
                    continue
                except ValueError:
                    dead += 1
                    continue
                ok = S.m == m and S.G == (s, c, d)
                match += ok
                wrong += not ok
            print(f'{p:>6} {p % 4:>4} {kind:<10}{tried:>8}{square:>12}{match:>12}{dead:>22}')
    print('RESULT:', 'every constructed billiard has the requested m and G' if not wrong else f'{wrong} MISMATCHES')
    return 1 if wrong else 0


if __name__ == '__main__':
    raise SystemExit(main())
