# Two-way correspondence with elliptic curves

Behind the geometric group law ([geometric_group_law.md](geometric_group_law.md)) is an elliptic curve. This page covers both directions:

1. **Billiard → curve:** every billiard orbit lifts to a genus-1 curve, and reflection is adding a fixed point.
2. **Curve → billiard:** conversely, every point satisfying a square condition on an elliptic curve with full rational 2-torsion is the step of some billiard.

So "counting bounces" and the discrete logarithm on this class of curves are the same problem. Proofs are in [proofs.md](proofs.md); numerical checks are given in each section.

## 1. State curve and lift

Notation as in [geometric_group_law.md](geometric_group_law.md) ($\epsilon,A,B,\lambda$). Let

$$
J:\quad c^2+\epsilon s^2=1,\qquad d^2+ms^2=1,\qquad m=\frac{\epsilon(B-A)}{B-\lambda},
$$

| Name | Definition |
|---|---|
| projection | $\pi(s,c,d)=(a\,c,\ b\,s)$ |
| identity | $O=(0,1,1)$, projecting to the starting point $(a,0)$ |
| point of order 2 | $T=(0,1,-1)$ |
| addition | $(s_1,c_1,d_1)+(s_2,c_2,d_2)=\dfrac{\big(s_1c_2d_2+s_2c_1d_1,\ c_1c_2-\epsilon s_1s_2d_1d_2,\ d_1d_2-ms_1s_2c_1c_2\big)}{1-\epsilon ms_1^2s_2^2}$ |
| negation | $-(s,c,d)=(-s,c,d)$, mirroring the impact in the $x$-axis |

For $\epsilon=+1$ the addition is `conic_billiard/jacobi.py`. The $\epsilon=-1$ formula comes from the substitution $s\mapsto is$, $m\mapsto-m$. This is an isomorphism of curves over $\mathbb F_p(i)$ fixing $O$, so it preserves the group law (Silverman III.4.8); the coefficients of the formula lie in $\mathbb F_p$, so $J(\mathbb F_p)$ is a subgroup.

**The step $G$.** For the second impact $P_1$, take $G=(s,c,d)$ with $\pi(G)=P_1$, choosing $d$ so that
$$
a^2m(c-1)+\epsilon K(1+d)=0,\qquad K=a^2+b^2(m-\epsilon).
$$
This $d$ is unique and $G$ lies on $J$ (L4). For $\epsilon=+1$ it is the point returned by `jacobi.ellipse_lift`.

**Lift theorem.** By L0, L3, L4 and L10 of [proofs.md](proofs.md): as long as every step satisfies the nondegeneracy conditions below, the reflection orbit is $P_k=\pi([k]G)$. This holds **for both ellipses and hyperbolas**. At the exceptional states where the affine formula is undefined, the lift is still defined; numerically it agrees with reflection everywhere and crosses every exceptional state (section 7 of [geometric_group_law.md](geometric_group_law.md)).

**Pencil members.** The member tangent to chords with "difference" $W$ is

$$
\mu(W)=\frac{-(A-B)(m-\epsilon)(c-1)}{a^2m(c-1)+\epsilon K(1+d)},
$$

so $\mu(G)=\infty$, the caustic.

| Fact | Lemma |
|---|---|
| the chord $\pi(U)\pi(V)$ is tangent to $C_{\mu(U-V)}$ and $C_{\mu(T-(U+V))}$ | L3, L9 |
| $\mu(W)=\mu(W')$ if and only if $W'=\pm W$ (fibre lemma) | L5 |
| the member is degenerate if and only if $W$ has order 2 | L8 |

## 2. Precise conditions for the geometric theorems

| Theorem | Conditions |
|---|---|
| A (step = translation by $W$) | $W\notin J[2]$; $\pi(X-W),\pi(X),\pi(X+W)$ pairwise distinct; $\pi(X)\notin C_{\mu(W)}$; nonzero denominators in the formula |
| B (the wrong root) | the result is $[a]\,(T-Q)$: $T-[a]Q$ for odd $a$ (same impact, opposite direction), $-[a]Q$ for even $a$ (mirror image) |
| C (common root) | the other roots of the two chords differ, $\mu(T-(U_1+V_1))\neq\mu(T-(U_2+V_2))$, i.e. $T-(U_1+V_1)\neq\pm(T-(U_2+V_2))$ |
| D (ladder) | every prefix $m$ satisfies A and C; the main failures are $[2m+2]G=O$ or $[2m]G=O$, $[m]G$ or $[m+1]G$ of order 2, and coincident chord endpoints; when $T\in\langle G\rangle$ there are a few further cases involving $T$ |

`ladder_hypotheses_fail()` in `tools/poncelet_pencil.py` implements these conditions one by one.

## 3. Weierstrass model

$$
x=\frac{(1+c)(1+d)}{s^2},\qquad y=\frac{(1+c)(1+d)(c+d)}{s^3}
$$

satisfy $y^2=x(x-\epsilon)(x-m)$ (L7). This is an **isomorphism**, not merely an isogeny:

- Near $O$: $s$ is a local parameter, $x\approx4/s^2$, $y\approx8/s^3$.
- At the other three points with $s=0$: $1+c\approx\epsilon s^2/2$ or $1+d\approx ms^2/2$, so $x,y$ are finite.
- At the points at infinity of $J$: $c,d$ have the same order as $s$, so $x,y$ are finite.

So $x$ has its only pole, of order 2, at $O$, and $y$ its only pole, of order 3, at $O$; since $0,\epsilon,m$ are pairwise distinct, Silverman III.3.1 makes this an isomorphism sending $O$ to the point at infinity. [group_order_search.md](group_order_search.md) computes group orders with this model.

This also shows that the group order of $J$ is always divisible by 4: all three points of order 2 are rational.

## 4. From curves back to billiards

**Which curves.** Suppose all 2-torsion of $E/\mathbb F_p$ is rational, i.e. $y^2=(x-e_1)(x-e_2)(x-e_3)$. Translating one root to 0 and scaling $x$ by a square brings it to $y^2=x(x-\epsilon)(x-m)$. The condition is that some difference $e_j-e_i$ lies in $\pm(\mathbb F_p^\times)^2$:

- for $p\equiv3\pmod4$ this always holds;
- for $p\equiv1\pmod4$ one of the three differences must be a square.

Curves without full rational 2-torsion (for example the prime-order secp256k1 and P-256) cannot be realized this way.

**Which point.** For a point $W=(s,c,d)$ on $J$ with $1+d\neq0$, take $a=1$ and

$$
b^2=-\frac{m(c-1)+\epsilon(1+d)}{\epsilon(m-\epsilon)(1+d)}=-\frac{c+d}{(1+c)(m-\epsilon)} .
$$

If this is a nonzero square, the billiard starting at $(a,0)$ with second impact $(a\,c,\ b\,s)$ has state-curve parameter $m$ and step $G=W$ (L11).

By L12 and L7, $b^2\equiv-(m-\epsilon)\,\delta(W)$ (modulo squares), where $\delta$ is the 2-descent map of section 6. So the square condition is

$$
\delta(W)\equiv-(m-\epsilon)\pmod{\text{squares}} .
$$

- $\delta$ is a homomorphism and is invariant under odd multiples. $\#J$ is a multiple of 4, so any $u$ coprime to the group order is odd, and **passing to a multiple $[u]P$ does not change the condition**.
- On an odd-order subgroup $\delta\equiv1$. So either every point of the subgroup is realizable or none is, depending only on whether $-(m-\epsilon)$ is a square.
- If not, the only option is another Weierstrass normalization, i.e. another root as origin or another $\epsilon$, giving another pair $(\epsilon, m)$. Whether some normalization always works is not yet proved ([open_questions.md](open_questions.md) Q3).

**Numerical checks** (`tools/realize_curve.py`). For $p=1009,1019,2003,2011$ and $\epsilon=\pm1$, 200 random $(m,W)$ per case, about half of which satisfy the square condition. This reflects $\delta$ splitting roughly evenly on random points. Every case that satisfies the condition and can be checked by step-by-step reflection gives the same $m$ and $G$, with no mismatches. The remaining orbits meet affine exceptional states (isotropic normal or intersection at infinity) and were not checked by step-by-step reflection; L11 itself is an algebraic identity and does not depend on this check.

**Reduction.** For curves with full rational 2-torsion that satisfy the condition above, the ECDLP with base $W$ and the bounce-counting problem of the corresponding billiard reduce to each other. Three points matter in an actual reduction:

1. **An impact determines the state only up to $T-Q$.** This is unique when $T\notin\langle G\rangle$; otherwise one extra trial is needed.
2. **The inverse map from a Weierstrass point back to $J$**: $\alpha=y/(x-m)=(1+c)/s$, $\beta=y/(x-\epsilon)=(1+d)/s$, $s=2\alpha/(\alpha^2+\epsilon)$, $c=\alpha s-1$, $d=\beta s-1$; $\alpha^2+\epsilon=0$ gives an impact at infinity. Implemented in `tools/dead_ends.py`.
3. **"Counting bounces" is defined through the lift.** Many realized billiards (most hyperbolas, and ellipses with $p\equiv1\pmod4$) meet affine exceptional states. There the bounce count is defined by the continuation of the Poncelet step, i.e. the turning rules of [model.md](model.md), which agree numerically with the continuation of the lift.

Security implications: [SECURITY.md](../SECURITY.md).

## 5. What the algebraic coordinates are for: Jacobi formulas and speed

The geometric group law and the Jacobi addition formulas compute the same group. The algebraic coordinates have three practical uses:

- **Proofs:** every identity in [proofs.md](proofs.md) is a polynomial statement on $J$.
- **Group orders:** existing point-counting algorithms apply through the Weierstrass model; see [group_order_search.md](group_order_search.md).
- **Speed:** about 3× faster than the geometric ladder in pure Python.

`tools/poncelet_bench.py`: 20 random $n$ per prime, both methods computed and compared; all equal.

| Field | Bits of $n$ | Jacobi doubling (ms) | Geometric ladder (ms) | Ratio |
|---|---:|---:|---:|---:|
| $p=10007$ | 14 | 0.031 | 0.177 | 5.6× |
| $p=2^{127}-1$ | 127 | 2.03 | 6.52 | 3.2× |
| $p=2^{255}-19$ | 255 | 8.49 | 25.97 | 3.1× |

Timings are pure-Python, variable-time measurements; only the ratios are meaningful. Over the 127- and 255-bit fields the two methods agree value for value, which is also a correctness check far beyond the range of orbit enumeration.

**Compressed encoding.** A state can be stored as one field element plus a branch bit, and that bit is the direction. See [compressed_points.md](compressed_points.md).

## 6. Conic types and 2-descent

Let $\delta(W)=x(W)-\epsilon$, taken modulo squares, with $\delta(O)=1$ and, at the point of order 2 with $x=\epsilon$, namely $(0,-1,1)$, the value $\epsilon(\epsilon-m)$. This is the 2-descent (Kummer) map of $y^2=x(x-\epsilon)(x-m)$, a group homomorphism $J(\mathbb F_p)\to\mathbb F_p^\times/\mathbb F_p^{\times2}$ (Silverman X.1.4; the homomorphism property holds over any field of characteristic not 2), and $\delta(T)=m-\epsilon$.

By L7, the coefficients of a pencil member satisfy

$$
A\,B\,q_0q_1\big(C_{\mu(W)}\big)\equiv\delta(W)\pmod{\text{squares}},
$$

and the other root gives $\delta(T-W)\equiv(m-\epsilon)\,\delta(W)$. Also $AB\equiv\epsilon$, and $C_\mu$ is of elliptic type if and only if $-q_0q_1$ is a non-square. This gives the two conclusions of section 6 of [geometric_group_law.md](geometric_group_law.md). In addition, $\delta$ is identically $+1$ on elements of odd order.

Check [5] of `tools/poncelet_pencil.py` verifies the types of the two roots on seven orbits (all consistent), and the homomorphism property of $\delta$ and associativity of addition on about 900 random pairs of points.
