# Conic Billiard Crypto

Conic billiards over a finite field: a ball bouncing inside an ellipse or a hyperbola modulo $p$. The project studies two questions:

- how much group arithmetic can be done with billiard geometry alone;
- whether "the number of bounces" can serve as a discrete-logarithm-style secret.

> **Research code, not a cryptographic library.** No parameter set in this repository has been validated for deployment, and the small-field DH demo recovers its own private key with BSGS. See [SECURITY.md](SECURITY.md) for the threat model and known limitations.

## Quick start

Requires Python 3.10+. The core package and the tests have no third-party dependencies; the optional tools are managed with [uv](https://docs.astral.sh/uv/).

```bash
python -m unittest discover -s tests -v

python -m conic_billiard scan --p 103 --kind all       # enumerate the state graph, periods and prime factors
python -m conic_billiard trace --kind hyperbola --p 103 --A 5 0 --B 98 0 --steps 8
python -m conic_billiard dh                            # DH demos: F_1009 (broken by BSGS) and a 255-bit group with 32-byte public keys
python -m conic_billiard bench                         # step-by-step reflection vs 4x4 matrices vs fast jumps

python tools/poncelet_pencil.py                        # geometric group law: numerical checks of the step, directions and ladder
python tools/dead_ends.py                              # how exceptional states are crossed (turning rules)
python tools/realize_curve.py                          # build a billiard from an elliptic curve
python tools/dictionary_check.py                       # billiard translation of security conditions (numerical checks)
python tools/poncelet_bench.py                         # geometric ladder vs Jacobi doubling
uv run --group tools python tools/poncelet_proof.py    # computer-algebra proofs (requires sympy)
```

## Core idea: one geometric step

The construction takes three steps ([model.md](docs/model.md)):

1. Start from the **closed real ellipse**: equal-angle reflection orbits are tangent to a confocal caustic.
2. Rewrite reflection as a **Poncelet step**: from the impact point, draw the other tangent to the caustic. This is a projective formula that uses only field operations.
3. **Take it as the definition** and extend it to finite fields, hyperbolas and impacts at infinity.

The boundary is $E: x^2/a^2+\epsilon\,y^2/b^2=1$ ($\epsilon=+1$ for an ellipse, $-1$ for a hyperbola), and the ball starts at $(a,0)$.

1. **Reflection.** The algebraic equal-angle reflection uses only addition, subtraction, multiplication and division, so it works over any odd prime field ([model.md](docs/model.md)).
2. **Caustic and pencil.** Every orbit is tangent to the same confocal caustic $C$. Linear combinations of $E$ and $C$ form a pencil of conics $C_\mu=E+\mu C$.
3. **The step formula.** From an impact $X$, draw a tangent to $C_\mu$ and intersect it with $E$ again. When one tangent is known, Vieta's formulas give the other directly. **Reflection is the step with $\mu=\infty$.**
4. **Addition.** For each step length $b$, one member $C_b$ of the pencil is tangent to every chord $P_kP_{k+b}$: a "$b$-step ruler". The same step along $C_b$ adds $b$ to any state, and the right member can be read off a chord without knowing $b$.
5. **Direction.** Each impact point carries two states, one for each direction of travel. They are the two roots of a quadratic equation.
6. **Ladder.** Doubling plus addition reaches the $n$-th impact in $O(\log n)$ steps, with no square roots.

Proofs and numerical checks: [geometric_group_law.md](docs/geometric_group_law.md) and [proofs.md](docs/proofs.md).

## Two-way correspondence with elliptic curves

- **Billiard → curve:** states live on the genus-1 curve $J: c^2+\epsilon s^2=1,\ d^2+ms^2=1$, reflection is translation by a fixed point $G$, and $J$ is isomorphic to $y^2=x(x-\epsilon)(x-m)$. This holds for both ellipses and hyperbolas. The exceptional states of the affine formula are turning points of the orbit: a ball that flies off to infinity comes back from $-B$ (on the other branch, for a hyperbola), and a ball that meets an isotropic normal pauses for one step and retraces its path. Both rules follow from the Poncelet step and agree numerically with the lift.
- **Curve → billiard:** on a curve with full rational 2-torsion, a point $W$ is the step of some billiard whenever it satisfies a square condition ($\delta(W)\equiv-(m-\epsilon)$). In an odd-order subgroup the condition depends only on the curve.

Counting bounces is therefore the ECDLP on this class of curves, no harder and no easier. In the other direction, elliptic-curve security conditions translate into billiard language: the group order is the number of states sharing a caustic, $4\mid N$ comes from the symmetries of the boundary, and circles, chords through a focus and parabolas degenerate ([security_dictionary.md](docs/security_dictionary.md)). The Jacobi addition formulas are algebraic coordinates for the same group. In Python they are about 3× faster than the geometric ladder, so the DH demo, the compressed encoding and the group-order tools use them. Details: [elliptic_correspondence.md](docs/elliptic_correspondence.md).

## Small-field regression results ($a=5, b=3$, re-checked by the tests)

| $p$ | Conic | Affine points | Longest period | Starting states with an undefined step |
|---:|---|---:|---:|---:|
| 103 | ellipse | 104 | 60 | 0 |
| 103 | hyperbola | 102 | 60 | 196 |
| 103 | parabola | 103 | 104 | 102 |
| 127 | ellipse | 128 | 126 | 0 |
| 127 | hyperbola | 126 | 74 | 244 |
| 127 | parabola | 127 | 128 | 126 |

"Undefined step" means the affine reflection formula is undefined at some step, not that the orbit stops; the Poncelet step gives a turning point there ([model.md](docs/model.md)). The longest period need not carry the largest prime factor: for $p=103$ the longest period, 60, has largest prime factor 5, while an orbit of period 58 carries the factor 29.

## Documentation

| Document | Contents |
|---|---|
| [docs/model.md](docs/model.md) | Reflection formula, ellipses and hyperbolas, exceptional states, point counts |
| [docs/geometric_group_law.md](docs/geometric_group_law.md) | **Main thread**: step formula, directions, addition, ladder, conic types, numerical checks |
| [docs/elliptic_correspondence.md](docs/elliptic_correspondence.md) | Lift, Weierstrass isomorphism, billiards from curves, Jacobi formulas and speed |
| [docs/proofs.md](docs/proofs.md) | Proof method, hypotheses and lemmas L0–L16 |
| [docs/security_dictionary.md](docs/security_dictionary.md) | Security conditions in billiard language: degeneracy, group order = state count, symmetry and $4\mid N$, attack models |
| [docs/group_order_search.md](docs/group_order_search.md) | Group-order screening and 127/255/256-bit samples |
| [docs/compressed_points.md](docs/compressed_points.md) | State encoding as one field element plus a direction bit |
| [docs/open_questions.md](docs/open_questions.md) | Open questions |
| [docs/references.md](docs/references.md) | References |

## Project layout

```
conic_billiard/
  core.py         finite-field reflection, 4x4 state matrices, classification of exceptions
  search.py       state graphs, cycles and prime-factor statistics
  jacobi.py       Jacobi addition on the state curve (algebraic coordinates, faster)
  compression.py  state encoding as one field element + one direction bit
  demos.py        two DH demos and a speed comparison
  __main__.py     CLI: scan / trace / dh / bench
tests/            unit tests and known-value regressions
tools/
  poncelet_pencil.py   numerical checks of the geometric group law (step, directions, rational ladder)
  poncelet_proof.py    computer-algebra proofs (sympy)
  dead_ends.py         exceptional states and the lift
  dictionary_check.py  billiard translation of security conditions (numerical checks)
  realize_curve.py     build a billiard from an elliptic curve
  poncelet_bench.py    geometric ladder vs Jacobi doubling
  conic_field_order_scan.py / hasse_bsgs_cardinality.py / sage_large_prime_scan.py /
  verify_large_prime_samples.py      group-order screening and large-prime samples
results/          group-order candidates and large-prime sample data
```

## Group-order tools (optional)

```bash
uv sync --group tools          # installs sympy (see dependency-groups in pyproject.toml)
uv run --group tools python tools/conic_field_order_scan.py --primes 103 127 257 1009 1019 65537
python tools/verify_large_prime_samples.py
sage -python tools/sage_large_prime_scan.py --prime curve25519 --trials 30
```

## Scope

1. **Hardness equals the ECDLP.** The billiard problem is equivalent to the ECDLP on the associated curve, so it adds no new hardness assumption. Group orders are always multiples of 4.
2. **The geometric proofs cover axis-aligned ellipses and hyperbolas with orbits through $(a,0)$.** Other cases are listed in [open_questions.md](docs/open_questions.md).
3. **The DH demo uses ellipses only.** The hyperbola lift is proved and the group-order screening covers hyperbolas ([group_order_search.md](docs/group_order_search.md)), but `jacobi.py` and the DH demo are not yet wired for $\epsilon=-1$.
4. **The geometric ladder is slower than the Jacobi formulas.** Its role is to show that the group structure lives in the billiard itself.

Earlier experiments with cubic-curve mirrors and real billiards were removed; the full version is kept at the git tag `archive/full-v0.2.0`.

License: [Apache-2.0](LICENSE). Copyright 2026 Jun-Xiang Zhan.
