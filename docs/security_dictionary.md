# Security conditions in billiard language

"Counting bounces" and the ECDLP on the state curve reduce to each other ([elliptic_correspondence.md](elliptic_correspondence.md)), so every elliptic-curve security condition can be restated in billiard terms. This page:

- states as theorems, with proofs, the translations that can be proved;
- lists in a comparison table at the end those that are not proved or have only arithmetic meaning.

Notation follows [model.md](model.md) and [geometric_group_law.md](geometric_group_law.md): $A=a^2$, $B=\epsilon b^2$, starting point $P_0=(a,0)$, second impact $P_1$, caustic parameter $\lambda$, state curve $J$ with parameter $m=\epsilon(B-A)/(B-\lambda)$. $\chi$ is the quadratic character of $\mathbb F_p$ ($+1$ on squares, $-1$ on non-squares, $0$ at zero).

```bash
python tools/dictionary_check.py                         # numerical checks of every theorem on this page
uv run --group tools python tools/poncelet_proof.py      # algebraic lemmas L13–L16 (and L0–L12)
```

## Summary

| Elliptic-curve condition | Billiard statement | Status |
|---|---|---|
| singular curve ($m\in\{0,\epsilon,\infty\}$) | the boundary is a circle, or the chord passes through one of the four foci | Theorem 1 |
| singular group becomes multiplicative | circle: a fixed rotation; through a focus: composition of two focal involutions; parabola: a fixed linear recurrence | Theorem 2 |
| group order $\#J$ | number of states sharing the caustic $=2\cdot\text{outside}+\text{on}$ | Theorem 3 |
| group order of the quadratic twist | $2\cdot\text{inside}+\text{on}$ | Theorem 3 |
| Hasse bound | $\lvert\text{outside}-\text{inside}\rvert\le2\sqrt p$ | Corollary 3.1 |
| supersingular (MOV, embedding degree $\le2$) | as many impacts inside the caustic as outside | Corollary 3.2 |
| anomalous ($\#J=p$, Smart's attack) | cannot occur | Corollary 3.3 |
| $4\mid\#J$, full rational 2-torsion | the two axis mirrors and the central symmetry of the boundary | Theorem 4 |
| point orders, small subgroups | orbit periods, short-period billiards | Theorem 5 |
| invalid-curve attack | a received state not on the agreed caustic | Theorem 5 |
| twist attack (impact only) | an impact inside the caustic | Theorem 5 |
| embedding degree, CM, generic attacks | arithmetic of the state count only | comparison table |

## 1. Degeneracy

**Theorem 1 (singular cases).** Suppose $P_1\neq P_0$ and the chord $P_0P_1$ is not isotropic ($u^2+v^2\neq0$). Then:

| $m$ | Condition | Billiard statement |
|---|---|---|
| $m=0$ | $A=B$ | the boundary is the circle $x^2+y^2=a^2$ |
| $m=\epsilon$ | $\lambda=A$ | the chord passes through $(0,\pm f')$, $f'^2=B-A$ |
| $m=\infty$ | $\lambda=B$ | the chord passes through $(\pm f,0)$, $f^2=A-B$ |

$(\pm f,0)$ and $(0,\pm f')$ are the four foci of the confocal family. For a real ellipse only the first two are real points; in $\mathbb F_p$, which pair is rational depends on whether $A-B$ and $B-A$ are squares.

*Proof.* All three cases can be read off $m=\epsilon(B-A)/(B-\lambda)$: $m=0$ if and only if $A=B$, $m=\epsilon$ if and only if $\lambda=A$, and the denominator vanishes if and only if $\lambda=B$.

$\lambda$ is defined by tangency: the line $ux+vy=w$ is tangent to $x^2/(A-\lambda)+y^2/(B-\lambda)=1$ if and only if $(A-\lambda)u^2+(B-\lambda)v^2=w^2$. For $u^2+v^2\neq0$ this gives a unique $\lambda$.

- For $\lambda=B$ the condition becomes $w^2=(A-B)u^2$, i.e. $w=\pm fu$, i.e. the line passes through $(\pm f,0)$. When $f\notin\mathbb F_p$, only $u=w=0$ remains, which is the major axis itself.
- The case $\lambda=A$ is the same with the axes swapped. $\square$

**Remarks.**

- **The "hyperbola" with $\epsilon=-1$ and $b^2=-a^2$ is a circle**, since then $B=A$. It exists only when $-1$ is a square.
- **The rectangular hyperbola ($a=b$, $\epsilon=-1$) is not degenerate.** Here $B=-A\neq A$, so $m\neq0$; it degenerates only on chords through a focus, like any other conic. Numerical check in `dictionary_check.py` [4]: for orbits with $p=1009,2003,10007$, lift and reflection agree, and periods such as 61 and 211 do not divide $p\pm1$.

  $x^2-y^2=1$ is the "unit circle" of the Lorentz metric, but this model uses the Euclidean inner product. When $-1$ is a square, the substitution $y\mapsto iy$ turns "circle + Euclidean metric" into "rectangular hyperbola + Lorentz metric", and "rectangular hyperbola + Euclidean metric" is a different combination.
- **The hypothesis $K\neq0$ holds automatically.** It simplifies to $K=a^2+b^2(m-\epsilon)=\lambda(B-A)/(B-\lambda)$. $\lambda=0$ would mean the chord is tangent to $E$ itself, which a chord through two distinct impacts cannot be.

**Theorem 2 (dynamics after degeneration).** In each of the following cases, counting bounces becomes a discrete logarithm in a subgroup of $\mathbb F_p^\times$ or $\mathbb F_{p^2}^\times$. Index calculus solves those in subexponential time, so all of them are weak.

**(a) The circle $x^2+y^2=r^2$.**

1. The next impact is the mirror image of $P_{k-1}$ in the line $OP_k$ (L15).
2. View a point $(x,y)$ as the element $z=x+iy$ of $R=\mathbb F_p[X]/(X^2+1)$. When $-1$ is a square, $R\cong\mathbb F_p\times\mathbb F_p$; otherwise $R=\mathbb F_{p^2}$. The norm is $z\bar z=x^2+y^2=r^2$.
3. The mirror in the line $Oz_k$ is $z\mapsto z_k^2\bar z/r^2$, and $\bar z_{k-1}=r^2/z_{k-1}$, so $z_{k+1}=z_k^2/z_{k-1}$.
4. Hence the ratio $\rho=z_k/z_{k-1}$ is invariant, $z_k=z_0\rho^k$, and $\rho$ lies in the norm-1 group $U$ with $\lvert U\rvert=p-\chi(-1)$.

Moreover, the affine formula is defined everywhere on the circle: the normal is parallel to $B$, and $B\cdot B=r^2\neq0$. If a chord from $B$ were isotropic, it would be perpendicular to $B$, which in two dimensions forces the chord to have length zero.

**(b) Chords through a focus ($\lambda=B$ or $\lambda=A$).**

1. By L0, the reflected line is tangent to the same degenerate member, i.e. it passes through one of the two foci. It is a different line from the incoming one, both pass through the impact, and the foci are not on $E$, so the reflected line passes through the **other** focus.
2. The orbit therefore alternates through $F_1,F_2$, and every two steps apply $\tau=\sigma_{F_2}\circ\sigma_{F_1}$, where $\sigma_F$ is "join to $F$ and intersect $E$ again", a projective involution of $E\cong\mathbb P^1$.
3. $F_1\neq F_2$, so $\tau\neq\mathrm{id}$. Its fixed points are the two intersections of the line through the foci with $E$: $(\pm a,0)$ for $\lambda=B$, $(0,\pm\sqrt B)$ for $\lambda=A$. A non-identity element of $\mathrm{PGL}_2$ has at most two fixed points, so these are all.
4. The two fixed points are distinct, so $\tau$ is diagonalizable:
   - when the fixed points are rational, $\tau$ is conjugate to $t\mapsto\kappa t$ with $\kappa\in\mathbb F_p^\times$;
   - when they are conjugate, $\kappa$ lies in the norm-1 subgroup of $\mathbb F_{p^2}^\times$.

So the state period divides $2(p-1)$ (for $\lambda=B$, since $(\pm a,0)$ are always rational) or $2(p-\chi(B))$ (for $\lambda=A$). Over the reals, $\tau$ is hyperbolic, which is why in the real billiard "orbits through a focus converge to the major axis".

**(c) The parabola $y=x^2$.**

1. The $x$-coordinates of the impacts satisfy $x_{k+1}=hx_k-x_{k-1}$, where $h=2+4k$ and $k=(x_0-x_1)^2/(1+4x_0x_1)$ is invariant along the orbit (L16).
2. The characteristic equation $z^2-hz+1=0$ has discriminant $h^2-4=16k(1+k)$. For $k(1+k)\neq0$ the two roots $z,z^{-1}$ are distinct and $x_k=\alpha z^k+\beta z^{-k}$.
3. So the period divides $p-\chi(k(1+k))$.

Reflection uses only the inner product and is unchanged when the inner product is scaled by a constant. So every parabola whose axis direction is not isotropic can be moved to $y=x^2$ by a reflection-preserving similarity.

**Numerical checks** (`dictionary_check.py` [1]–[3]):
- Circle: $p=103,1009,1019,2017$, 200 sampled states; mirror formula and periods all match.
- Through a focus: ellipses and hyperbolas with $a=5,b=3$, 9 primes; every chord passes through the other focus. Periods divide the predicted values; both the split and non-split cases are covered.
- Parabola: the recurrence holds throughout, and the periods all match.

Orbits that meet affine exceptions (impacts at infinity, isotropic normals) do not yield a full period and are left out of the period checks. The proofs themselves work on the projective line and are not affected by these exceptions.

## 2. Group order is the state count

For an impact $X$ on the boundary (including impacts at infinity), the tangents from $X$ to the caustic fall into three cases:

- **outside**: two rational tangents, i.e. two possible directions of departure;
- **on**: $X$ lies on the caustic, with one tangent;
- **inside**: the two tangents are conjugate over $\mathbb F_{p^2}$, with no rational direction of departure.

**Theorem 3.** Suppose $m\notin\{0,\epsilon,\infty\}$ and let $J^t$ be the quadratic twist of $J$. Then

$$
\#J(\mathbb F_p)=2\cdot\#\text{outside}+\#\text{on},\qquad
\#J^t(\mathbb F_p)=2\cdot\#\text{inside}+\#\text{on},
$$

and the three kinds of impact total $p+1$. In other words, the group order of the elliptic curve is **the total number of billiard states sharing this caustic**.

*Proof.*

1. **Tangent discriminant.** The directions of tangents from $X$ to $q_0x^2+q_1y^2=1$ form a binary quadratic form with discriminant $\Delta(X)=q_0q_1\,C(X)$ (L13, first half). For the caustic, $q_0=1/(A-\lambda)$ and $q_1=1/(B-\lambda)$. For an impact $X=\pi(s,c)=(ac,bs)$ we also have $C(X)=\lambda(1-ms^2)/(A-\lambda)$ (L13, second half). Hence
   $$\Delta(X)\equiv\lambda(B-\lambda)\,(1-ms^2)\pmod{\text{squares}}.$$
2. **The constant factor is a square.** At $P_0$ ($s=0$) the chord $P_0P_1$ is a rational tangent, so the quadratic has a rational root, the other root is rational too, and $\Delta(P_0)$ is a square. Also $\lambda\neq0$ (remark to Theorem 1), so $\lambda(B-\lambda)$ is a nonzero square and $\chi(\Delta(X))=\chi(1-ms^2)$.
3. **Affine part.** When $J$ is projected to $(s,c)$, the values of $d$ above each impact satisfy $d^2=1-ms^2$, which has $1+\chi(1-ms^2)$ rational solutions; $\pi$ is a bijection from $c^2+\epsilon s^2=1$ to $E$.
4. **Impacts at infinity.** These exist only when $\chi(-\epsilon)=1$.
   - On the state curve: the projective closure $c^2+\epsilon s^2=z^2$, $d^2+ms^2=z^2$ at $z=0$ gives $d^2=-ms^2$, with $1+\chi(-m)$ points.
   - On the billiard: the parallel lines with an asymptotic direction $(x,y)$ tangent to the caustic are rational if and only if $(A-\lambda)y^2+(B-\lambda)x^2$ is a square. Substituting $y^2=-Bx^2/A$, this value is $x^2\lambda(B-A)/A$.
   - The two agree: after multiplying by the square $\lambda(B-\lambda)$, both sides become $(B-A)(B-\lambda)$, up to the factor $-\epsilon$, and $\chi(-\epsilon)=1$.
5. **Group order.** For $m\notin\{0,\epsilon\}$ the intersection of these two quadrics is a smooth projective curve isomorphic to the Weierstrass model (L7). Summing gives $\#J=\sum_X(1+\chi(\Delta(X)))=2\cdot\#\text{outside}+\#\text{on}$.
6. **Twist.** Replacing $d^2$ by $g\,d^2$ ($g$ a non-square) gives a curve isomorphic to $J$ over $\mathbb F_{p^2}$, its quadratic twist. Its point count is $\sum_X(1-\chi(\Delta(X)))=2\cdot\#\text{inside}+\#\text{on}$, and the two add up to exactly $2(p+1)$. $\square$

**Corollary 3.1 (Hasse).** outside $-$ inside $=\#J-(p+1)$, so $\lvert\text{outside}-\text{inside}\rvert\le2\sqrt p$.

**Corollary 3.2 (supersingular).** For $p>3$, $J$ is supersingular if and only if $\#J=p+1$, that is, **as many impacts lie inside the caustic as outside**. The embedding degree is then at most 2, and the MOV attack moves the discrete logarithm to $\mathbb F_{p^2}$, so such a caustic is weak.

**Corollary 3.3 (no anomalous case).** $4\mid\#J$ (Theorem 4) and $p$ is odd, so $\#J\neq p$ and Smart's attack never applies.

**Numerical checks** (`dictionary_check.py` [5]): ellipses and hyperbolas, $p=1009,1019,2017,2027$, 6 caustics each. All $p+1$ impacts, including impacts at infinity, are counted one by one, and both formulas hold in every case.

## 3. Why the group order is a multiple of 4

**Theorem 4.** All three points of order 2 of $J$ are rational, and their translations act on impacts as the three symmetries of the boundary:

| Point of order 2 | Translation | Action on impacts |
|---|---|---|
| $T=(0,1,-1)$ | $(s,c,d)\mapsto(-s,c,-d)$ | mirror in the $x$-axis |
| $(0,-1,1)$ | $(s,c,d)\mapsto(-s,-c,d)$ | central symmetry |
| $(0,-1,-1)$ | $(s,c,d)\mapsto(s,-c,-d)$ | mirror in the $y$-axis |

So for any orbit $P_k=\pi([k]G)$, the mirrored orbit $\sigma(P_k)=\pi([k]G+T_\sigma)$ is the orbit starting from the state $T_\sigma$. By Lagrange's theorem $4\mid\#J$. In billiard terms: the Klein four-group $\mathbb Z/2\times\mathbb Z/2$ of boundary symmetries acts freely on the states sharing one caustic.

*Proof.* The translation formulas are L6 and L14. Substituting into $\pi(s,c,d)=(ac,bs)$ gives the right-hand column, and this holds for every point of $J$, so the mirror image of a whole orbit is the translated orbit. $\square$

**Numerical checks** (`dictionary_check.py` [6]): 8 orbits on ellipses and hyperbolas, including orbits of period 4; all three symmetries match.

## 4. Attack models

**Theorem 5.**

1. **Point order is orbit period.** Suppose $\mathrm{ord}(G)>2$.
   - If $[k]G=O$, the state clearly returns to the start.
   - Conversely, if the state returns to the start, then from $P_k=P_0$ and $\pi(U)=\pi(V)\iff V\in\{U,T-U\}$ (L6) we get $[k]G\in\{O,T\}$.
   - If $[k]G=T$, then $P_{k-1}=P_{-1}$ gives $T-G\in\{-G,\ T+G\}$; the first gives $T=O$ and the second $2G=O$, both impossible. So $[k]G=O$.

   Hence the orbit period equals $\mathrm{ord}(G)$. By L11, every point $Q$ satisfying the square condition is the step of some billiard, and $[k]Q$ is that billiard after $k$ bounces. **A small-subgroup attack is sending a short-period billiard.**
2. **Invalid-curve attack.** The caustic determines $m$, and $m$ determines the curve. If a received state is not on the agreed caustic, the computation happens on another curve $J_{m'}$, whose group order is given by Theorem 3 for that caustic and may be smooth. Input validation therefore means: **check that the state's caustic is the agreed one.**
3. **Twist attack.** When only the impact is sent, or the compressed encoding (one field element plus a direction bit, see [compressed_points.md](compressed_points.md)) is used, decoding must compute $d=\pm\sqrt{1-ms^2}$. By step 2 of the proof of Theorem 3, a solution exists if and only if the impact lies outside or on the caustic. **An impact inside the caustic has no rational state**; it belongs to the twist $J^t$. Without this check, the peer gets the computation done on the twist, and what leaks is governed by the twist's group order $2\cdot\#\text{inside}+\#\text{on}$. The decoder in `compression.py` rejects such inputs here ("quartic value is not a quadratic residue"). $\square$

## 5. Comparison table (not proved, or arithmetic only)

| Condition | Billiard statement | Why it is only a comparison |
|---|---|---|
| embedding degree $k$ (MOV / Frey–Rück) | the large prime factor $r$ of the state count divides $p^k-1$ | arithmetic of $N$ and $p$ only. The one case with a geometric statement is the supersingular case of Corollary 3.2 |
| CM discriminant, endomorphism ring | none found | no known billiard meaning |
| generic attacks (Pollard rho, BSGS) | about $\sqrt r$ group operations | generic group algorithms, independent of the geometry |
| unified degeneracy criterion | the boundary and the caustic do not meet in four distinct points (over the algebraic closure) | the circle (meeting only at the circular points), chords through a focus (degenerate caustic) and the parabola (tangent at infinity) all fit, but "if and only if" is not proved |
| pairs of lines (degenerate boundary) | the ball alternates between two lines, and its direction is multiplied by the fixed rotation $H_2H_1$ every two steps | argued only, not checked with a tool; such boundaries are outside the `Conic` model |
| isotropic starting chord ($u^2+v^2=0$) | no confocal caustic | the tangency condition does not involve $\lambda$; such orbits are outside the model and `dead_ends.py` skips them |
| orbits not starting at $(a,0)$ | — | see [open_questions.md](open_questions.md) Q1 |

## Conclusion

The translation works in both directions, and most conditions have proved billiard statements:

- **degeneracy**: circles, chords through a focus, parabolas;
- **group order**: counting impacts inside and outside the caustic;
- **multiple of 4**: the symmetry group of the boundary;
- **input validation**: the right caustic, and impacts outside it.

What does not translate are purely arithmetic conditions such as the embedding degree and CM. They still have to be checked on the group order $N$, which now has a billiard meaning.
