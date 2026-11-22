# Cubic mirror billiards (三次橢圓曲線的真正鏡面反射)

**Status: exploratory algebraic dynamics, not a DH construction.**

We use the *smooth cubic boundary* rather than elliptic-curve chord-and-tangent addition:

\[
E: F(x,y)=y^2-x^3-ax-b=0\quad\text{over }\mathbb F_p,\quad p>3,
\]

with \(4a^3+27b^2\ne0\). For the experiments, \(a=-1,b=0\), so \(E:y^2=x^3-x\). The projective completion adds the unique point \(\mathcal O=[0:1:0]\). It **does not** magically create affine rational intersection points where the remaining roots lie in a field extension.

## Mirror reflection and the 2-valued successor relation

For two **distinct affine rational** points \(A,B\), let \(v=B-A\) and

\[
N_B=(-3x_B^2-a,2y_B),\quad H_B=I-2\frac{N_BN_B^T}{N_B^TN_B},\quad w=H_Bv.
\]

If \(N_B^TN_B=0\), the formula is undefined: the Euclidean-like bilinear form has **isotropic** directions in many finite fields, even when the curve is nonsingular. We do not fill in these states arbitrarily.

A possible next intersection has the form \(C=B+tw\) for \(t\ne0\), with

\[
\boxed{-w_x^3t^2+(w_y^2-3x_Bw_x^2)t+
 (2y_Bw_y-(3x_B^2+a)w_x)=0.}
\]

There are 0, 1 or 2 **distinct affine rational** next points. The \(t=0\) contact is excluded from the list of *new* hits. When \(w_x=0\) there is a projective intersection at \(\mathcal O\) and generally one other finite root. Reflection itself is reversible as a **relation**: if \((A,B)\to(B,C)\), then \((C,B)\to(B,A)\) is among the reversed candidates, when the formulas are defined. This is not a one-to-one dynamical map without further choices.

Finite fields do not order ray parameters, so a rule such as "nearest positive intersection" is not naturally available. Lexicographically selecting the smaller/larger candidate is a valid software experiment but is neither geometrically canonical nor generally invertible.

## Verified exhaustive experiments

| Quantity | \(p=101\) | \(p=1009\) |
|---|---:|---:|
| Affine points | 103 | 1,039 |
| Ordered distinct \((A,B)\) | 10,506 | 1,078,482 |
| Two next points | 5,164 | 550,626 |
| No affine rational next point | 5,050 | 520,550 |
| One next point | 192 | 2,122 |
| Contact at B (reported separately) | 100 | 1,032 |
| Isotropic normal | 0 | 4,152 |
| Reverse-validated directed transitions | 10,620 | 1,104,406 |

Status categories are mutually exclusive, but contact-at-B cases may also have a different valid next point.

For \(p=101\), \(A=(0,0),B=(2,39)\) yields reflected \(w=(98,1)\) and \(27t^2+48t+10=0\), with nonsquare discriminant 12: no \(\mathbb F_{101}\) next hit. With \(A=(1,0),B=(2,39)\), the two next points are \((50,35)\) and \((96,48)\).

For \(p=1009\), four points of the cubic have an isotropic normal; each makes the reflection formula fail for all 1,038 distinct starting points A.

## Some invariant tests and the failure of naive fast jumping

Every valid transition preserves the Legendre-symbol square class of chord length:

\[
\chi\big((B-A)\cdot(B-A)\big)=\chi\big((C-B)\cdot(C-B)\big),
\]

because \(C-B=tH_B(B-A)\) and \(H_B^TH_B=I\). We checked all directed transitions for both primes. This is a **one-trit** invariant (0/±1), not an elliptic caustic or a fast-jump group law.

We also solved \(I(A,B)=I(B,C)\) for all polynomials \(I\) in four affine state coordinates of total degree at most 4, using 220 sampled valid transition equations over each field. The sampled linear-system nullities for degrees \(1,2,3,4\) were \(1,1,3,11\) for **both** primes, precisely the dimensions spanned by constants and the obvious equations \(F(A)=F(B)=0\) and their allowed degree-4 multiples. Since these are always global invariants and the sample puts an upper bound on the global invariant space, this proves **no nontrivial polynomial first integral of total degree ≤4** in this particular algebraic state representation and for these two test fields. This does not rule out rational, high-degree, or algebraic invariants.

If we forcibly define a deterministic map by selecting the lexicographically **minimum** candidate after each hit, the longest cycle in \(\mathbb F_{101}\) is 4; choosing the **maximum** yields a longest cycle of 14. These are *arbitrary selection policies*; neither is a canonical billiard group or DH candidate.

The reflection relation is generically 2-valued, so \(T^n(X)\) is not even a single point before a branch policy is defined. No direct \(O(\log n)\) fast jump is established here. A reversible, globally defined branch-lift should be designed **before** claiming a useful cyclic group.

## How to run

```bash
python -m conic_billiard.cubic_experiment --out results/cubic_mirror.json
python -m unittest discover -s tests -v
```

## References

- Max Weinreich, *The Dynamical Degree of Billiards in an Algebraic Curve*, Journal of Geometric Analysis 35 (2025). https://arxiv.org/abs/2305.14287 — general algebraic billiards are \((d-1)\)-to-\((d-1)\) correspondences in degree \(d\).
- Max Weinreich, *Algebraic Billiards in the Fermat Hyperbola*, https://arxiv.org/abs/2406.16172 — algebraic dynamics and multivalued correspondences on non-quadratic curves.
