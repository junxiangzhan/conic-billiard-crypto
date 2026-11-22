# v0.2.0 — compressed lifted-state encoding

- Added one-field-element + branch-bit canonical encoding (`conic_billiard/compression.py`) for affine generalized Jacobi states (`epsilon=±1`).
- Added special-chart codes for two points with `c=-1`, nonresidue/pole/padding rejection, and variable-time modular square roots.
- Added a 32-byte encoded-public-key toy DH exchange for the Curve25519-related prime subgroup, with a demonstration subgroup-membership check.
- Added exhaustive small-field encoder round-trip tests and large-prime sample regression tests.
- Retained ellipse/hyperbola mirror reflection, 4×4 state matrices, algebraic cubic mirror experiments, trajectory/group-order searches, CI, and large-prime verification scripts.

**Unresolved:** complete projective exceptional-point encoding; production-safe constant-time group arithmetic; proven robust hyperbola lift; real-world security certification. Do **not** deploy for cryptography.
