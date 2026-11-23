# Unreleased — geometric group law

- Documentation is now in English (README, SECURITY.md, docs/), with neutral wording in descriptions and docstrings. Stale notes fixed: circles and parabolas are treated in `docs/security_dictionary.md`, not as open questions.
- License changed from MIT to Apache-2.0 (`LICENSE`, `NOTICE`, `pyproject.toml`). `toy_dh()` renamed to `small_field_dh()`; the `dh` CLI key is now `F1009_breakable`.
- Docs follow one route: closed real ellipse -> Poncelet step (other tangent to the caustic, a projective formula) -> that step taken as the definition, extended to finite fields, hyperbolas and impacts at infinity (`docs/model.md`, README).
- New result: the affine dead ends are turning points of the Poncelet step. Impact at infinity from state (A,B): B, inf, -B, -A, ... (point-symmetric orbit; for a hyperbola -B is on the other branch). Isotropic normal: B, B, A, ... (the orbit retraces itself). Geometric arguments in `docs/model.md`; `tools/dead_ends.py` now checks both rules against the lift (all pass) and counts turning points per orbit (always two of the same kind, even period; not proved in general). `tools/dead_ends.py` takes `--primes`.
- New `docs/security_dictionary.md`: elliptic-curve security conditions in billiard language, with proofs. Singular state curve (m in {0, eps, inf}) iff the boundary is a circle or the chord passes through one of the four foci; circle, focus orbits and parabola reduce to a multiplicative group (rotation, composition of two focal involutions, fixed linear recurrence); #J = 2*outside + on and #twist = 2*inside + on, counting impacts by whether their tangents to the caustic are rational (so Hasse, supersingular = "as many impacts inside as outside", anomalous impossible); the three 2-torsion translations are the two axis mirrors and the central symmetry; small-subgroup, invalid-curve and twist attacks as short orbits, a wrong caustic and impacts inside the caustic. The rectangular hyperbola is not degenerate. Unproved or purely arithmetic items (embedding degree, CM, the unified degeneracy criterion, line pairs, isotropic start chords) are listed as a comparison table.
- New lemmas L13–L16 in `tools/poncelet_proof.py` (tangent discriminant, caustic at an impact, the other two 2-torsion translations, circle and parabola reflection); new `tools/dictionary_check.py` (numerical checks of every theorem in the dictionary).
- Corrections: the group-order scan does cover hyperbolas (README and open questions said otherwise; only the DH demo is ellipse-only); `demos.py` no longer calls the hyperbola lift an open task; `conic_field_order_scan.py` says isomorphic, not isogenous.

- Docs restructured around the geometric method: one formula (tangent to a pencil member, intersect the boundary) gives reflection (mu = infinity), addition and the O(log n) ladder. New `docs/model.md`, `docs/geometric_group_law.md`, `docs/elliptic_correspondence.md`, `docs/proofs.md`; `docs/poncelet_pencil.md` and `docs/hyperbola.md` were merged into them. Jacobi formulas are now presented as the algebraic coordinates (proofs, group orders, ~3x faster implementation).
- Math in Markdown now uses `$...$` / `$$...$$`.
- Open questions: removed the solved geometric-doubling question, the fixed-matrix question (answered in SECURITY.md by the embedding-degree argument) and the statistics question; kept other starts/general conics, complete handling of exceptional states, and parameter screening.
- New proofs: the pencil step at mu = infinity is the reflection (L10); every point W of the state curve with a square b^2 is the step G of a billiard with the same m (L11), and b^2 ~ -(m - eps)*delta(W) mod squares (L12), so in an odd-order subgroup realizability depends only on the curve.
- SECURITY.md: exact embedding degree of the Curve25519 subgroup ((l-1)/6, from the full factorisation of l-1).
- New tools: `tools/dead_ends.py` (wherever the affine reflection is undefined, the lift continues and every orbit closes; reflection and lift agree everywhere else) and `tools/realize_curve.py` (elliptic curve -> billiard construction).

- Added `tools/poncelet_proof.py`: computer-algebra proofs (sympy, exact Groebner reduction) of the identities behind the pencil construction, for ellipses and hyperbolas: reflection = second tangent to the caustic; closed form of the pencil member for a difference W; chord pi(U)pi(V) tangent to the member for U-V; the step G; T - W = (s,c,-d); Weierstrass model y^2 = x(x-eps)(x-m) and the conic-type formula; the other tangent member of a chord is mu(T-(U+V)); mu(W) = -1 only at 2-torsion (with the conic-type formula: for W not O or T, the member is degenerate only at 2-torsion).
- Added `tools/poncelet_pencil.py`: numerical checks on 3 ellipse and 4 hyperbola orbits, and a square-root-free Montgomery-style O(log n) ladder. No wrong values; every ladder exception lies inside the set where the proved sufficient conditions fail (equal to it on 6 of the 7 orbits).
- Added `tools/poncelet_bench.py`: the pencil ladder agrees with Jacobi double-and-add on random n up to 255 bits and is about 3x slower in pure Python.
- New result: hyperbola orbits starting at `(a,0)` lift to `c^2 - s^2 = 1, d^2 + m s^2 = 1`, with reflection = translation by G.
- Corrections: ladder exceptions are not limited to n near the identity unless n is first reduced mod ord(G); whether the conic type replaces the direction bit depends on whether `m - eps` is a square, not on `p mod 4`.
- Packages are managed with uv: sympy is in the `tools` dependency group (`uv sync --group tools`); `uv.lock` added.

# v0.3.0 — simplification

- Removed the cubic-curve mirror, real-oval first-hit and fixed-branch experiments (kept in git tag `archive/full-v0.2.0`).
- Merged `toy_dh.py`, `compressed_dh_demo.py` and `benchmark.py` into `demos.py`; run them with `python -m conic_billiard dh` and `bench`.
- Merged Conic tests into `tests/test_conic.py` and DH tests into `tests/test_dh.py`.
- Merged `docs/security.md` into `SECURITY.md` and `docs/large_prime_samples.md` into `docs/group_order_search.md`; shortened README; removed regenerable files from `results/` and `GITHUB_SETUP.md`.

# v0.2.0 — compressed lifted-state encoding

- Added one-field-element + branch-bit canonical encoding (`conic_billiard/compression.py`) for affine generalized Jacobi states (`epsilon=±1`).
- Added special-chart codes for two points with `c=-1`, nonresidue/pole/padding rejection, and variable-time modular square roots.
- Added a 32-byte encoded-public-key demonstration DH exchange for the Curve25519-related prime subgroup, with a demonstration subgroup-membership check.
- Added exhaustive small-field encoder round-trip tests and large-prime sample regression tests.
- Retained ellipse/hyperbola mirror reflection, 4×4 state matrices, algebraic cubic mirror experiments, trajectory/group-order searches, CI, and large-prime verification scripts.

**Unresolved:** complete projective exceptional-point encoding; production-safe constant-time group arithmetic; proven robust hyperbola lift; real-world security certification.
