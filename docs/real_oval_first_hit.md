# Real bounded-oval billiard vs finite-field cubic mirror correspondence

## Motivation

The finite-field *algebraic* mirror relation on a smooth cubic is generically two-valued. It would be incorrect to treat both extra intersections as physical successive collisions. A genuine Euclidean billiard follows a **directed ray** and collides with its boundary at the **first positive-distance intersection**. The order condition exists over the real numbers but not over finite fields.

## Real geometry

Use the smooth real cubic

\[
E: y^2=x^3-x.
\]

Its real locus has two connected components:
- a bounded **strictly convex** oval with \(-1\le x\le 0\);
- an unbounded component with \(x\ge 1\).

A convenient parametrization of **the real oval only** is

\[
x=-\sin^2\theta,\qquad y=\sin\theta\cos\theta\sqrt{1+\sin^2\theta},\quad 0\le\theta<\pi.
\]

This parametrization is used *only to generate test points*, not to define the reflection algorithm. Strict convexity follows, for example, from the concavity of the upper graph \(y=\sqrt{x^3-x}\) for \(-1<x<0\):

\[
y''=\frac{3x}{y}-\frac{(3x^2-1)^2}{4y^3}<0.
\]

For incoming boundary points \(A,B\), with interior chord \(v=B-A\), put

\[
N_B=(1-3x_B^2,2y_B),\qquad
w=v-2\frac{v\cdot N_B}{N_B\cdot N_B}N_B.
\]

The outgoing ray is \(R(t)=B+tw\), **\(t>0\)**. Substituting in the cubic and factoring out \(t=0\) yields

\[
-w_x^3t^2+(w_y^2-3x_Bw_x^2)t+
\bigl(2y_Bw_y-(3x_B^2-1)w_x\bigr)=0.
\]

The **physical** successor is the smallest positive real root, equivalently the unique next intersection with the bounded oval (for non-tangent, interior-directed rays). The remaining finite line intersection belongs to the unbounded real component, or the third projective intersection may be the point at infinity. In exact Euclidean billiard dynamics, reversing velocity reverses the trajectory.

## Numerical experiment

Command:

```bash
python -m conic_billiard.real_oval_experiment
```

Using deterministic seed `20261118`, 10,000 random, non-near-tangent ordered pairs \((A,B)\) on the real oval produced:

| Property | Result |
|---|---:|
| First forward hit on oval found | 10,000 / 10,000 |
| Third intersection behind B | 5,227 |
| Third intersection beyond first hit, in front of B | 4,772 |
| Third finite intersection numerically ill-conditioned / treated as linear | 1 |
| Maximum one-bounce reverse error (double precision) | 1.58e-13 |
| Maximum equation residual | 4.45e-13 |

A further 200 initial pairs each underwent 80 successive bounces. Reverse-roundoff statistics using ordinary double-precision floating point:

| Forward and reverse depth | Median maximum position error over 200 initial pairs | Fraction exceeding 1e-4 |
|---|---:|---:|
| 10 | 4.93e-15 | 0% |
| 20 | 3.22e-14 | 0% |
| 40 | 1.15e-12 | 0% |
| 80 | 1.29e-7 | 30.5% |

Numerical sensitivity does **not** imply mathematical non-invertibility. An independent, optional 100-digit `mpmath` check on a selected unstable orbit achieves a maximum reverse error about `7.16e-83` after 80 bounces in each direction.

```bash
pip install mpmath           # optional, only for the high-precision check
python -m conic_billiard.real_oval_high_precision
```

## Why this does not automatically define a finite-field map

Over \(\mathbb F_p\), there is no order compatible with field operations, so predicates such as \(t>0\), \(t_1<t_2\), or \(-1<x<0\) do not have field-invariant counterparts. The topological real-oval component is likewise not available as a canonical finite-field subset.

A concrete exact witness:

- \(A=(-1,0), B=(0,0)\), reflected direction \(w=(-1,0)\);
- reflected line parameter \(t\) has other roots \(t=1,-1\);
- over \(\mathbb R\), first-hit is \(t=1\), giving \(C=(-1,0)\); the \(t=-1\) point \((1,0)\) is behind the ray;
- over \(\mathbb F_{101}\), the two roots are \(1\) and \(100\). Both are equally legitimate field elements; choosing based on their integer representatives invents an extrinsic convention.

Exhaustive \(\mathbb F_{101}\) algebraic relation statistics for all 10,506 distinct ordered affine pairs: 5,164 pairs have two extra rational intersections; 5,050 no extra affine rational intersection; 192 have one; 100 contact at the current point. Contact is separately categorized. These counts refer to **algebraic line intersections**, not physical first-hit outcomes.

Some missing finite-field roots lie in \(\mathbb F_{p^2}\), not at the unique projective point at infinity. Projective completion alone does not solve this.

## Remaining research question

Find a **precisely defined**, reversible finite-field lift with explicit branch data and no ad-hoc ordering, then analyze its periodic orbits and whether it has a fast-jump law. It must not be assumed to inherit the real physical first-hit rule or elliptic-curve-group security. A candidate that uses an order on integer residues can be a valid algorithm, but it is not a canonical algebraic continuation of specular billiards.

## Sources

- M. Weinreich, *The Dynamical Degree of Billiards in an Algebraic Curve* (2025), https://arxiv.org/abs/2305.14287: cubic algebraic billiards fit a multi-valued algebraic correspondence rather than the real first-hit map.
- Standard convex billiard geometry: the physical collision is the first forward ray intersection with the bounding convex oval.
