# Geometric group law: one formula for reflection, addition and fast jumps

The core of this project is one purely geometric step: **from an impact point, draw a tangent to some conic, then intersect it with the boundary again.** Reflection is a special case. With a different conic, the same formula "adds an arbitrary state", and together with doubling it reaches the $n$-th impact in $O(\log n)$ steps. The whole process uses only field operations: no square roots and no curve coordinates.

The background is Poncelet's classical theorem. This project proves every identity it needs in its own coordinates by computer algebra ([proofs.md](proofs.md)) and checks the boundary cases numerically.

```bash
python tools/poncelet_pencil.py                                   # numerical checks, p=1009 ellipse (default)
python tools/poncelet_pencil.py --kind hyperbola --p 2003 --B 2 1289
python tools/dead_ends.py                                         # how exceptional states are crossed
uv run --group tools python tools/poncelet_proof.py               # computer-algebra proofs (requires sympy)
```

## 1. Notation

The model is in [model.md](model.md): $\epsilon=\pm1$, $A=a^2$, $B=\epsilon b^2$, boundary $E: x^2/A+y^2/B=1$, starting point $P_0=(a,0)$, second impact $P_1$.

| Name | Definition |
|---|---|
| caustic $C$ | $x^2/(A-\lambda)+y^2/(B-\lambda)=1$, with $\lambda$ fixed by "the chord $P_0P_1$ is tangent to it" |
| pencil | $C_\mu:\ E+\mu\,C=0$ (sum of the two quadratic forms, both with constant term $-1$); $\mu=\infty$ stands for $C$ itself |
| member coefficients | $C_\mu$ written as $q_0x^2+q_1y^2+q_2=0$ |

## 2. The step formula

**The step $\mathrm{step}_\mu(Y,X)$.** Given a point $X$ on the boundary and a line from $X$, tangent to $C_\mu$, that leads to $Y$: take the **other** tangent from $X$ to $C_\mu$ and return its second intersection $Z$ with $E$.

Concretely, the line through $X$ with direction $v=(v_1,v_2)$ is tangent to $C_\mu$ if and only if

$$
(X^TQ_\mu v)^2-(X^TQ_\mu X)(v^TQ_\mu v)=A_fv_1^2+2B_fv_1v_2+C_fv_2^2=0 .
$$

One root is known, $(x_1,y_1)=Y-X$; by Vieta's formulas the other is $d=(C_f\,y_1,\ A_f\,x_1)$. Intersect along $d$ with $E$:

$$
t=-2\,\frac{X^TQ_Ed}{d^TQ_Ed},\qquad Z=X+t\,d .
$$

In matrix form, $Z=X+t\,D\,(Y-X)$ with $D=\begin{pmatrix}0&C_f\\A_f&0\end{pmatrix}$. This has the same shape as the reflection matrix $M(X)$ of [model.md](model.md):

$$
\begin{bmatrix}Z\\X\end{bmatrix}=\begin{bmatrix}I-tD&tD\\I&0\end{bmatrix}\begin{bmatrix}X\\Y\end{bmatrix}.
$$

**Reflection is the step with $\mu=\infty$.** For $\mu=\infty$ (the caustic itself), $d$ is parallel to the reflected direction $w=H_Bv$, and both use the same scale-free formula for the second intersection. So $\mathrm{step}_\infty(P_{k-1},P_k)=P_{k+1}$, and ordinary reflection is an instance of this formula ([proofs.md](proofs.md) L0, L10).

## 3. States, directions and "differences"

**Direction.** There are two tangents from an impact to the caustic: one leads to the next impact, the other to the previous one. So an impact point alone does not describe a state; we also need to record which tangent to take. A pair of consecutive impacts $(P_{k-1},P_k)$ carries exactly this information.

**Reading the "difference" off a chord.** For each $b$ there is a pencil member $C_b$ tangent to **all** chords $P_kP_{k+b}$. The two tangents from any impact to $C_b$ lead to the impacts $b$ steps ahead and $b$ steps back, so $C_b$ is a "$b$-step ruler" that does not distinguish sign. From the chord $P_0P_b$ alone, "tangent to $C_\mu$" is a quadratic equation in $\mu$, and solving it does not require knowing $b$. Both roots are rational functions of the state coordinates $(c,d)$ (L1). With only an impact and no direction, solving this quadratic takes one square root, and the sign of that root is the direction bit. The ladder avoids this square root through Theorem C in section 4.

**The two roots are the two directions.** The other root corresponds to the same impact with the opposite direction. Taking $a$ steps along the wrong root gives:

- $a$ odd: the same impact, opposite direction;
- $a$ even: the impact mirrored in the $x$-axis.

(L2, L6; in group language this is $[a]\,(T-Q)$, see [elliptic_correspondence.md](elliptic_correspondence.md).) When the two roots have different conic types is covered in section 6.

## 4. Addition

**Theorem A (one step is a translation).** Let $Y, X, Z$ be impacts on one orbit, $b$ steps apart, pairwise distinct, with $C_b$ nondegenerate and $X$ not on $C_b$. Then the two tangents from $X$ to $C_b$ lead exactly to $Y$ and $Z$, so $\mathrm{step}_{C_b}(Y,X)=Z$. (L3, L8; the precise conditions in group language are in [elliptic_correspondence.md](elliptic_correspondence.md).)

**Theorem C (choosing $C_b$ without square roots).** Two chords with "difference" $b$ share one root, $C_b$; their other roots are generally different (L9 and the fibre lemma). Eliminating $\mu^2$ from the two quadratics $A_i\mu^2+B_i\mu+C_i=0$ gives

$$
\mu=\frac{A_1C_2-A_2C_1}{A_2B_1-A_1B_2},\qquad\text{or, when the denominator is zero,}\qquad \mu=\frac{C_1B_2-C_2B_1}{C_2A_1-C_1A_2}\ (\text{possibly }\infty).
$$

## 5. Doubling and the ladder

The state $(P_m,P_{m+1})$ walks through the binary digits of $n$, using only the fixed $P_0, P_1$ and $P_{-1}$, obtained by one backward reflection:

- $C_m$: the common member of the chords $P_0P_m$ and $P_1P_{m+1}$;
- $C_{m+1}$: the common member of the chords $P_0P_{m+1}$ and $P_{-1}P_m$ (needed only for bit 1);
- bit 0: $P_{2m}=\mathrm{step}_{C_m}(P_0,P_m)$, $P_{2m+1}=\mathrm{step}_{C_m}(P_1,P_{m+1})$;
- bit 1: $P_{2m+1}=\mathrm{step}_{C_m}(P_1,P_{m+1})$, $P_{2m+2}=\mathrm{step}_{C_{m+1}}(P_0,P_{m+1})$.

This has the same structure as the Montgomery ladder on elliptic curves. Every step is the step formula of section 2, with a different conic; reflection itself is the case $C_1=C_\infty$.

**Theorem D (when the ladder is guaranteed correct).** If every $m$ visited satisfies the conditions of Theorems A and C, the ladder outputs $P_n$. This is a sufficient condition: when it fails, the ladder may raise an exception or may still be correct. The main failure cases are:

- the other roots of the two chords coincide;
- the two endpoints of a chord coincide;
- the corresponding member is degenerate.

In the group these are conditions such as $[2m+2]G=O$, $[2m]G=O$ or points of order 2; see [elliptic_correspondence.md](elliptic_correspondence.md).

**Practical rule (argued, not formally proved): reduce $n$ below the orbit period $N$ first.** Then only $n\in\{N-2,N-1\}$ for even $N$, and a very few $n$ related to points of order 2, degenerate; for odd $N$ none do. Without reduction, random large $n$ often degenerate: for $N=134$, 197 of 300 random 64-bit $n$ raise an exception.

## 6. Conic types and the direction bit

A diagonal conic $q_0x^2+q_1y^2+q_2=0$ over $\mathbb F_p$ has one of two types: with $-q_0q_1$ a non-square it has no points at infinity (elliptic type); with $-q_0q_1$ a square it has two (hyperbolic type). By L7:

- The two roots have different types **if and only if $m-\epsilon$ is a non-square** ($m$ is the state-curve parameter, fixed by the caustic). In that case the type alone selects the right root, with no extra direction bit.
- If the peer's point is known to lie in an odd-order subgroup, the right root satisfies $q_0q_1\equiv\epsilon$ (modulo squares). It is therefore of elliptic type if and only if $-\epsilon$ is a non-square.

This does not depend on $p \bmod 4$. The proof uses 2-descent; see [elliptic_correspondence.md](elliptic_correspondence.md).

## 7. Numerical checks

**Ladder and directions** (`tools/poncelet_pencil.py`, all with $a=5,b=3$, starting at $(5,0)$). The ladder was tested on $n=2,\dots,4N-1$ and on 300 random 64-bit $n$:

| Conic | $p$ | $P_1$ | $N$ | Lift | Chain | Two roots ↔ two directions | Ladder: correct / wrong / exception | Exceptions all predicted | Type |
|---|---:|---|---:|---|---|---|---|---|---|
| ellipse | 1009 | (730, 932) | 134 | ✓ | 38/38 | 532/532 | 613 / 0 / 221 | ✓ | 130/130 |
| ellipse | 2003 | (7, 1491) | 64 | ✓ | 37 (1 exception) | 518 (14 exceptions) | 258 / 0 / 296 | ✓ | 60/60 |
| ellipse | 10007 | (9, 3304) | 4992 | ✓ | 38/38 | 532/532 | 20236 / 0 / 30 | ✓ | 4988/4988 |
| hyperbola | 1009 | (2, 966) | 510 | ✓ | 38/38 | 517 (15 exceptions) | 2261 / 0 / 77 | ✓ | 506/506 |
| hyperbola | 1019 | (2, 620) | 104 | ✓ | 33 (5 exceptions) | 451 (81 exceptions) | 369 / 0 / 345 | ✓ | 98/98 |
| hyperbola | 2003 | (2, 1289) | 498 | ✓ | 38/38 | 532/532 | 2194 / 0 / 96 | ✓ | 494/494 |
| hyperbola | 10007 | (2, 8564) | 843 | ✓ | 38/38 | 532/532 | 3637 / 0 / 33 | ✓ | 840/840 |

- **Wrong values:** none.
- **Exceptions vs prediction:** every exception falls where the conditions of Theorem D fail. On six orbits the two sets are identical; on the hyperbola with $p=1019$, 12 additional $n$ are predicted, where the conditions fail but the ladder is still correct.
- **After reducing to $n<N$:**
  - the actual failures for the first three ellipses are $\{132,133\}$, $\{62,63\}$ and $\{4990,4991\}$;
  - for the hyperbola with $p=1019$ they are $\{100,\dots,103\}$;
  - the hyperbola with $p=10007$ ($N$ odd) has none.

**Exceptional states are crossed** (`tools/dead_ends.py`, all orbits starting at $(5,0)$):

| Conic | $p$ | Orbits | Reflection agrees with lift | Orbits meeting an exception | Kinds and counts of exceptional states |
|---|---:|---:|---:|---:|---|
| hyperbola | 103 | 50 | 1,733 | 14 | intersection at infinity 28 |
| hyperbola | 1009 | 502 | 143,520 | 161 | intersection at infinity 322 |
| hyperbola | 1019 | 507 | 179,827 | 253 | intersection at infinity 506 |
| hyperbola | 2003 | 1000 | 672,820 | 342 | intersection at infinity 684 |
| ellipse | 1009 | 501 | 108,203 | 287 | isotropic normal 448, coincident input points 448, intersection at infinity 126 |
| ellipse | 103 / 1019 / 2003 | 51 / 509 / 1001 | all agree | 0 | — |

Wherever reflection is defined, it agrees with the lift, without exception. At exceptional states the next step of the lift is always defined, and every exception follows the turning rules of [model.md](model.md) (consequences of the Poncelet step):

- intersection at infinity: impact sequence $B,\infty,-B,-A,\dots$, with the whole orbit point-symmetric about this step;
- isotropic normal: impact sequence $B,B,A,\dots$, where the ball stays for one step and retraces its path. The coincident input points that follow are part of the same turn.

Every orbit eventually closes. Every orbit that meets an exception has exactly two turning points of the same kind, and an even period. One further run with `--primes 10009`: the hyperbola there also shows isotropic normals, and all follow the rules.

These are numerical checks, restricted to orbits with $a=5, b=3$ starting at $(5,0)$. The tool skips two cases: the very few starting points where the lift formula itself fails (0 to 2 orbits per case), and steps whose input impact is at infinity and so cannot be checked against affine reflection.

## 8. Speed and where the code lives

The purely geometric ladder is correct, but in pure Python it is about 3× slower than doubling in Jacobi coordinates. Both compute the same group; numbers and the comparison are in section 5 of [elliptic_correspondence.md](elliptic_correspondence.md). The DH demo and the group-order tools currently use the faster Jacobi implementation; the geometric method lives in `tools/`.

## 9. Scope

- **Classical background and what is done here.** Poncelet's general theorem and the Chasles–Graves lemma are classical. This project proves the needed identities in its own coordinates, implements the rational ladder, and classifies the exceptions and directions.
- **The starting point is fixed at $(a,0)$.** Orbits that do not pass through $(a,0)$ are not treated.
- **Only axis-aligned ellipses and hyperbolas.** Circles and parabolas appear as degenerate cases in [security_dictionary.md](security_dictionary.md); rotated and general conics are not treated.
- **Security.** This is another algorithm for the same group; see [SECURITY.md](../SECURITY.md).
