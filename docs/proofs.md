# Proofs: method, hypotheses and lemmas

This page lists the lemmas proved by `tools/poncelet_proof.py` and explains why these computer-algebra checks count as proofs.

```bash
uv sync --group tools
uv run --group tools python tools/poncelet_proof.py      # about 3–4 minutes; prints RESULT at the end
```

## Method

Each lemma is written as an identity of rational functions. After clearing denominators, the numerator is reduced modulo the ideal

$$
\big(c_i^2+\epsilon s_i^2-1,\ \ d_i^2+ms_i^2-1\big)
$$

in lex order with $c_i,d_i$ before $s_i$.

- **A zero remainder is a proof.** The leading terms $c_i^2, d_i^2$ of these generators are pairwise coprime, so they form a Gröbner basis (Buchberger's first criterion), and a zero remainder proves that the numerator lies in the ideal.
- **Valid over every $\mathbb F_p$.** The generators are monic in $c_i,d_i$, so the reduction never divides by an integer. An identity found over $\mathbb Q$ therefore holds over $\mathbb Z$, and hence over $\mathbb F_p$ for every odd prime $p$, provided the cleared denominators are nonzero.
- **No reliance on simplification.** The script uses `factor()` only for printing; every decision uses `expand`/`together`/`cancel` plus Gröbner reduction.
- **Negative controls.** An independent review deliberately broke the formulas, for example $D\to D+1$ or $\mu(U-V)\to\mu(U+V)$, and the checks failed.

## Hypotheses (nondegeneracy conditions)

$p$ is an odd prime; $a,b\neq0$; $A\neq B$; $\lambda\notin\{A,B\}$; $m\notin\{0,\epsilon\}$; $K=a^2+b^2(m-\epsilon)\neq0$; and the denominators cleared by each lemma are nonzero:

- $u^2+v^2$ for a chord;
- $n\cdot n$ for a normal (isotropic normal);
- $1+\epsilon t^2$ for the parametrization;
- $s$, $1+c$, $1+d$, $D$;
- $1-\epsilon ms_1^2s_2^2$ for addition.

When one of these denominators is zero, the numerical tools raise `UndefinedReflection` or `ZeroDivisionError` instead of silently patching the value.

## Lemmas (valid for both $\epsilon=+1$ and $\epsilon=-1$)

| Lemma | Statement | Used for |
|---|---|---|
| L0 | Take two boundary points by a rational parameter (covering every point except $(-a,0)$). If the caustic is tangent to the chord $X_0X_1$, the algebraic reflection of the ray at $X_1$ is also tangent to it; when the reflected line differs from the incoming line, it is "the other tangent". | step formula, lift |
| L10 | The direction $d=(C_fy_1, A_fx_1)$ of the pencil step at $\mu=\infty$ is parallel to the reflected direction; both use the same scale-free formula for the second intersection, so **this step is reflection**. (Proved for the generic branch $x_1,y_1\neq0$; for $x_1=0$ or $y_1=0$ the code uses the other form of the Vieta formula, checked only numerically.) | step formula |
| L1 | For the tangency quadratic of the chord $\pi(O)\pi(W)$, the discriminant on $J$ equals $(R\,d)^2$, so both roots are rational functions. | reading the difference off a chord |
| L2 | $\mu(W)=N/D$ is one root; the other is $d\to-d$; $\mu$ does not involve $s$, so $\mu(-W)=\mu(W)$. | directions |
| L3 | For any $U,V\in J$, the chord $\pi(U)\pi(V)$ is tangent to $C_{\mu(U-V)}$. | Theorem A |
| L9 | The other tangent member of that chord is $C_{\mu(T-(U+V))}$. | Theorem C |
| L4 | For any boundary point $P_1$, the $G$ with $\pi(G)=P_1$ and $D(G)=0$ is unique and lies on $J$ ($\mu(G)=\infty$). | lift |
| L11 | Conversely, for $W$ on $J$, taking $b^2=-a^2\,\frac{m(c-1)+\epsilon(1+d)}{\epsilon(m-\epsilon)(1+d)}$ gives a billiard whose parameter is still $m$. This $b^2$ is exactly the solution of $D(W)=0$, so its $G$ is $W$. | curve → billiard |
| L12 | The $b^2$ above equals $-a^2(c+d)/((1+c)(m-\epsilon))$; with L7, $b^2\equiv-(m-\epsilon)\,\delta(W)$ (modulo squares). | realizability condition |
| L5 | $N, D$ are linear in $(c,d)$, and both vanish at $(1,-1)$ (the image of $T$). | fibre lemma |
| L6 | $T+(s,c,d)=(-s,c,-d)$, so $T-W=(s,c,-d)$; $2T=O$. | directions, Theorem B |
| L7 | $x=(1+c)(1+d)/s^2$, $y=(1+c)(1+d)(c+d)/s^3$ satisfy $y^2=x(x-\epsilon)(x-m)$; and $AB\,q_0q_1(C_{\mu(W)})=(x-\epsilon)\cdot\big(Ks(1+d)/((1+c)D)\big)^2$. | Weierstrass model, conic types |
| L8 | $N+D=\epsilon K(c+d)$ and $c^2-d^2=(m-\epsilon)s^2$, so $\mu(W)=-1$ if and only if $s=0$. With L7: for $W\notin\{O,T\}$, the member is degenerate if and only if $W$ has order 2. | Theorems A, D |
| L13 (second half) | On $c^2+\epsilon s^2=1$, the value of the caustic at an impact is $C(ac,bs)=\lambda(1-ms^2)/(A-\lambda)$. | group order = state count |
| L14 | $(0,-1,1)+(s,c,d)=(-s,-c,d)$ and $(0,-1,-1)+(s,c,d)=(s,-c,-d)$, both points of order 2; with L6, the translations by the three points of order 2 are the mirror in the $x$-axis, the central symmetry and the mirror in the $y$-axis. | symmetry and $4\mid\#J$ |

## Lemmas (independent of $\epsilon$)

| Lemma | Statement | Used for |
|---|---|---|
| L13 (first half) | The quadratic form in the direction of tangents from $X$ to $q_0x^2+q_1y^2=1$ has discriminant $q_0q_1\,C(X)$. | group order = state count |
| L15 | Circle $x^2+y^2=r^2$: the next impact is $2(X_0\cdot X_1)/(X_1\cdot X_1)\,X_1-X_0$, the mirror image of $X_0$ in the line $OX_1$. | circle degeneracy |
| L16 | Parabola $y=x^2$: from $(u,u^2),(v,v^2)$ the next impact has $x$-coordinate $(2+4k)v-u$ with $k=(u-v)^2/(1+4uv)$, and $k$ is unchanged for the next chord. | parabola degeneracy |

L13–L16 are used in [security_dictionary.md](security_dictionary.md).

**Fibre lemma** (from L5): $\mu=\kappa$ is a line through $(1,-1)$ in the $(c,d)$ plane. $(c,d)$ lies on the nondegenerate conic $\epsilon(d^2-1)=m(c^2-1)$, which passes through $(1,-1)$, so the line meets it in at most one further point; and $(c,d)$ determines $W$ up to sign. Hence for $W,W'\neq T$, $\mu(W)=\mu(W')$ if and only if $W'=\pm W$.

## Parts supported by argument or numerical checks only

- **Degenerate set after reduction:** "after reducing $n$ below the period only a few degenerate $n$ remain" is argued, not formally proved.
- **Ladder exception set:** the comparison of actual exceptions with the predicted set is numerical only (seven orbits).
- **Translation of security conditions:** Theorems 1–5 of [security_dictionary.md](security_dictionary.md) follow from L0–L16 plus standard facts (fixed points of $\mathrm{PGL}_2$ elements, the quadratic-character point-count formula, the Hasse bound); the items in its comparison table are not proved.
- **Exceptional states:** the turning rules (intersection at infinity: $B,\infty,-B$; isotropic normal: $B,B,A$) have geometric arguments ([model.md](model.md)). Their agreement with the lift, and exactly two turning points per orbit, are checked numerically only, on orbits starting at $(5,0)$ (`tools/dead_ends.py`).
- **Associativity of the group law:** for $\epsilon=-1$ it follows from the substitution argument, and it is also checked numerically on random points.

## Classical sources

The Chasles–Graves lemma (L0), Poncelet's general theorem (the geometric content of L3), Silverman III.3.1 and III.4.8 (isomorphism and group law), X.1.4 (2-descent). See [references.md](references.md).
