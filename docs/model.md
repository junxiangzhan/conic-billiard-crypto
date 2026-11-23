# Model: mirror reflection over a finite field

This page defines the billiard itself: the boundary, the reflection formula, exceptional states and point counts. Ellipses and hyperbolas over a finite field use the same formulas and differ only in a sign $\epsilon$.

## Route: from the closed ellipse to the hyperbola

1. **The closed ellipse.** The elliptic billiard in the real plane: equal-angle reflection, with every trajectory tangent to one conic confocal with the boundary, the caustic (Chasles–Graves).
2. **Extract the Poncelet formula.** Reflection equals "draw the **other** tangent from the impact point to the caustic, then intersect it with the boundary again" ([proofs.md](proofs.md) L0, L10). This description uses only tangents and intersections, not angles, so it is projective; as a formula it uses only field operations and therefore works over any odd prime field.
3. **Take it as the definition and extend.** With the Poncelet step as the definition of the billiard:
   - for $\epsilon=-1$ the formula is unchanged, giving the hyperbola;
   - the affine equal-angle formula is undefined at exceptional states, but the Poncelet step is still defined there, which gives two turning rules (see "Exceptional states" below). The ball can pass through a point at infinity and hit the other branch of the curve.

In the projective plane, ellipses, parabolas and hyperbolas are the same kind of curve; they differ only in whether the line at infinity meets them in 0, 1 or 2 points. The two branches of a hyperbola together with its two points at infinity form one closed loop. The same holds over finite fields; see "Projective view and point counts" at the end. Changing the pencil member in the same step gives addition and fast jumps; see [geometric_group_law.md](geometric_group_law.md).

## Boundary

Take an odd prime $p$, $a,b\in\mathbb F_p^\times$ and $\epsilon\in\{+1,-1\}$, and set $A=a^2$, $B=\epsilon b^2$:

$$
E:\quad \frac{x^2}{A}+\frac{y^2}{B}=1 .
$$

- $\epsilon=+1$ is the ellipse $x^2/a^2+y^2/b^2=1$.
- $\epsilon=-1$ is the hyperbola $x^2/a^2-y^2/b^2=1$.

The code also supports the circle and the parabola $y=x^2$ as control cases. Parabola impacts linearize to a fixed $2\times2$ matrix, so the discrete logarithm falls to a finite field and is weak. The circle also reduces to a multiplicative group; both cases are treated in Theorem 2 of [security_dictionary.md](security_dictionary.md).

## Reflection formula

Write a general conic as $F(X)=X^TQX+\ell^TX+c=0$. For impacts $A,B$ let $v=B-A$ and the normal $N=QB+\ell/2$. When the denominators are invertible:

$$
H_B=I-2\frac{NN^T}{N^TN},\qquad w=H_Bv,\qquad
\tau=-\frac{2B^TQw+\ell^Tw}{w^TQw},\qquad C=B+\tau w .
$$

This is the algebraic continuation of real equal-angle reflection. It uses the inner product $u^Tv$ rather than any notion of "angle" in the finite field. The whole formula uses only field operations, so it works over any field of odd characteristic.

**Geometric meaning.** In the real plane, the trajectory of equal-angle reflection is always tangent to one conic confocal with $E$, the caustic $C$. This remains true over finite fields, and reflection is exactly "the other tangent from the impact to the caustic"; see L0 in [proofs.md](proofs.md). Reflection is therefore a special case of a more general geometric step; see [geometric_group_law.md](geometric_group_law.md).

## States and 4×4 matrices

A state is a pair of consecutive impacts $X_k=(P_{k-1},P_k)$. One reflection step can be written as

$$
\begin{bmatrix} P_{k+1}\\P_k\end{bmatrix}
=\underbrace{\begin{bmatrix}I+\tau_kH_k & -\tau_kH_k\\I&0\end{bmatrix}}_{M(X_k)}
\begin{bmatrix}P_k\\P_{k-1}\end{bmatrix}.
$$

The matrix $M$ depends on the state, so $n$ steps cannot be computed as $M^n$. Fast jumps come from the pencil construction ([geometric_group_law.md](geometric_group_law.md)); for why no fixed small matrix exists, see [SECURITY.md](../SECURITY.md).

## Exceptional states

The affine formula cannot be evaluated in the following cases. The code raises `UndefinedReflection` instead of filling in a value.

| Code | Condition | Meaning |
|---|---|---|
| `coincident_input_points` | $A=B$ | two coincident points do not define an incoming line |
| `isotropic_normal` | $N^TN=0$ | the normal is isotropic and the mirror matrix has a zero denominator |
| `next_intersection_at_infinity` | $w^TQw=0$ | the other intersection of the reflected line with the boundary is at infinity |
| `repeated_tangent_intersection` | $\tau=0$ | the reflected line is tangent to the boundary at the impact |

**These are gaps in the affine formula, not true endpoints of the dynamics.** The Poncelet step is still defined there, and the two main exceptions are **turning points** of the orbit. Suppose the exception occurs at the state $(A,B)$, that is, the ball travels from $A$ to $B$:

| Exception | Impact sequence given by the Poncelet step | Symmetry of the orbit |
|---|---|---|
| intersection at infinity | $B,\ \infty,\ -B,\ -A,\ \dots$ | point symmetry about the step at infinity: $P_{k+1+j}=-P_{k+1-j}$ |
| isotropic normal | $B,\ B,\ A,\ \dots$ | the orbit retraces itself about the repeated step: $P_{k+1+j}=P_{k-j}$ |

The coincident input points $(B,B)$ that immediately follow an isotropic step are part of the same turn; the incoming line is the tangent to $E$ at $B$.

**Intersection at infinity.** The reflected direction $w$ satisfies $w^TQw=0$, so $w$ is an asymptotic direction of $E$, and the reflected line $\ell$ meets $E$ at $B$ and at the point at infinity $Q_\infty=[w:0]$.

1. The two tangents from $Q_\infty$ to the caustic are the two parallel lines with direction $w$ tangent to the caustic. The caustic is centrally symmetric, so they are $\ell$ and $-\ell$ (the tangency condition on the intercept is $k^2=\text{const}$, whose roots are negatives of each other).
2. $-\ell$ meets $E$ at $Q_\infty$ and $-B$, so the ball comes back from $-B$. For a hyperbola, $-B$ is on the other branch.
3. The state $(Q_\infty,-B)$ is the image of $(Q_\infty,B)$ under the central symmetry, and continuing from $(Q_\infty,B)$ walks the orbit backwards. The central symmetry preserves $E$, the caustic and the reflection formula, so the sequence continues $-A,-P_{k-2},\dots$.

**Isotropic normal.** When $N^TN=0$, the tangent $T_B$ to $E$ at $B$ (perpendicular to $N$) is also an isotropic line through one of the circular points $I$ or $J$. This happens only when $-1$ is a square, that is, $p\equiv1\pmod4$.

1. The tangents from a circular point to $E$ are common tangents of the whole confocal family, so $T_B$ is also tangent to the caustic.
2. The two tangents from $B$ to the caustic are therefore $BA$ and $T_B$. Intersecting $T_B$ with $E$ again gives $B$ (the point of tangency counts twice), so the ball stays for one step.
3. The incoming line of the state $(B,B)$ is $T_B$, and the other tangent is $BA$, so the ball returns to $A$.

**Turning points come in pairs (numerical observation).** In every orbit checked by `tools/dead_ends.py` that meets an exception, there are exactly two turning points, of the same kind, and the period is even. Part of the reason:

- The turning rules give the orbit a time-reversal symmetry, and a cycle of period $N$ has two positions on the symmetry axis under such a reflection.
- When $N$ is even, the second position is of the same kind as the first. An axis impact of the infinity kind satisfies $P=-P$, so it must be the point at infinity. An axis position of the isotropic kind is a repeated impact, and the only common tangents of $E$ and the caustic are the four isotropic ones, so it is also an isotropic step.
- Two symmetries of the same kind would compose to a period shorter than $N$, so there are no more.

The cases of odd $N$, or both kinds of turn on the same orbit, have not been ruled out by proof.

The argument above describes what the Poncelet step gives. Its agreement with the lift to the genus-1 curve (where every step is a group translation, defined everywhere) is currently a numerical check: `tools/dead_ends.py` checked all orbits with $a=5,b=3$ starting at $(5,0)$, and every exception follows the rules in the table above; numbers are in section 7 of [geometric_group_law.md](geometric_group_law.md). A general proof, and an implementation of the continuation in `core.py`, are listed as [open_questions.md](open_questions.md) Q2.

## Projective view and point counts

Let $[X:Y:Z]\in\mathbb P^2(\mathbb F_p)$; the projective closure is $X^2/A+Y^2/B=Z^2$. Every nondegenerate projective conic over $\mathbb F_p$ has $p+1$ rational points. The number of affine points is $p+1$ minus the number of points at infinity.

Points at infinity satisfy $Y^2=-(B/A)X^2$, so whether they are rational depends on whether $-B/A=-\epsilon b^2/a^2$ is a square, that is, whether $-\epsilon$ is a square:

| | $p\equiv3\pmod4$ | $p\equiv1\pmod4$ |
|---|---|---|
| ellipse ($\epsilon=+1$) | 0 points at infinity, $p+1$ affine points | 2, with $p-1$ affine points |
| hyperbola ($\epsilon=-1$) | 2, with $p-1$ affine points | 2, with $p-1$ affine points |

**Over a finite field, "ellipse" and "hyperbola" differ only by this sign; it is not a classification as in real geometry.** Both have genus-1 state curves; see [elliptic_correspondence.md](elliptic_correspondence.md).

## Small-field example

$p=103$, $a=5$, $b=3$. On the hyperbola, $A=(5,0)$ and $B=(98,0)=(-5,0)$ form an orbit of period 2. Of all $102\cdot101=10{,}302$ ordered starting states, 10,106 can compute the next step and 196 stop because the second intersection lies at infinity; the longest affine cycle has length 60.

```bash
python -m conic_billiard scan --p 103 --kind hyperbola
python -m conic_billiard trace --kind hyperbola --p 103 --A 5 0 --B 98 0
```
