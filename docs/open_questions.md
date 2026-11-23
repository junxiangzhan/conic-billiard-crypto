# Open questions

Solved questions have moved into the main documents:

- native geometric doubling: [geometric_group_law.md](geometric_group_law.md);
- whether a fixed matrix weakens the DLP: [SECURITY.md](../SECURITY.md);
- the genus-1 lift for ellipses and hyperbolas: [elliptic_correspondence.md](elliptic_correspondence.md).

The question "can the more random appearance be quantified" was dropped: the system is a deterministic group translation, and its statistical appearance has no bearing on security.

## Q1 — Orbits not starting on the major axis, and general conics

The current proofs cover only axis-aligned ellipses and hyperbolas, and orbits must pass through $(a,0)$, because the lift uses that point as the identity.

- Orbits not through $(a,0)$: can another point of the orbit serve as the base, with the same lift and pencil formulas? A tangent from $(a,0)$ to the caustic need not be rational; how should that case be handled?
- Rotated conics, or conics with linear terms: what is the state curve? Circles and parabolas are done; both degenerate to multiplicative groups (Theorem 2 of [security_dictionary.md](security_dictionary.md)).
- A unified degeneracy criterion: is "the boundary and the caustic do not meet in four distinct points" exactly the set of degenerate cases? Circles, chords through a focus and parabolas all fit, but there is no proof.

**Acceptance criteria:** give the explicit state curve and projection, prove that reflection is a translation, and use step-by-step reflection in `core.py` as an independent check.

## Q2 — Complete handling of exceptional states

For orbits through $(a,0)$, the numerical evidence ($a=5,b=3$, four primes) is that every exceptional state of the affine formula is just a gap: the lift crosses it and every orbit closes (section 7 of [geometric_group_law.md](geometric_group_law.md)). The exceptional states are intersections at infinity, isotropic normals, and the coincident input points that immediately follow.

The Poncelet step already gives a projective version that does not go through the lift: the two turning rules of [model.md](model.md), with geometric arguments, agreeing numerically with the lift. Still to do:

- prove that the turning rules agree with the lift, and give the version for general $a,b$;
- prove that every orbit meeting an exception has exactly two turning points of the same kind: rule out odd periods, and both kinds of turn on the same orbit;
- prove the same for **all** starting states (including the orbits of Q1);
- implement the turning rules in `core.py`, so that `trace` and `scan` can cross exceptional states;
- explain where the "starting states with an undefined step" counted by `scan` (for example the 196 for the hyperbola at $p=103$) land in the cycles after continuation.

**Acceptance criteria:** the continuation agrees with the original map on their common domain, is uniquely defined, and introduces no extra branches; where no continuation exists, describe the failure set explicitly.

## Q3 — Screening secure parameters

[group_order_search.md](group_order_search.md) can already screen for $N=4r$ samples, for both ellipses and hyperbolas (including a 61-bit hyperbola sample). Still missing:

- **A hyperbola DH demo:** the group law and the lift are proved, and `tools/` has an $\epsilon=-1$ implementation, but `conic_billiard/jacobi.py` and the DH demo use ellipses only.
- **Choice of normalization:** in a large prime subgroup, whether the curve → billiard construction (section 4 of [elliptic_correspondence.md](elliptic_correspondence.md)) works depends only on whether $-(m-\epsilon)$ is a square. For a curve with full rational 2-torsion, is there always a Weierstrass normalization $(\epsilon,m)$ for which it holds?
- **Dropping the direction bit:** in an odd-order subgroup with $m-\epsilon$ a non-square, the direction can be recovered from the conic type. Can this remove that bit from the public key without introducing a new ambiguity?
- **Standard ECC checks:** small subgroups, anomalous curves, embedding degree, twists, and input validation.

`largest_prime_factor_any_period` is a screening signal for small fields, not a security proof; a long cycle can also consist of small prime factors only.

## Basic reproducible checks

```bash
python -m unittest discover -s tests -v
python -m conic_billiard scan --p 103 --kind all
python tools/poncelet_pencil.py
python tools/dead_ends.py
uv run --group tools python tools/poncelet_proof.py
```

For every new continuation or jump formula, first check it against step-by-step reflection in `core.py` as an **independent reference**, before discussing performance or cryptographic use.
