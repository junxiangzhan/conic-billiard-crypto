# Finite-Field Conic Billiards: Group-Order Screening

Reproducible screening of group orders for the billiard state curves. Screening is the first stage of parameter selection; the remaining checks are listed under "Larger cryptographic primes".

## Frozen conic mirror model

Let `p` be an odd prime, `epsilon = +1` (ellipse) or `-1` (hyperbola), and `a=5`, `b=3` in the experimental scans. The boundary is

$$
C_\epsilon:\quad (x/a)^2 + \epsilon(y/b)^2=1 \pmod p.
$$

Choose boundary start `A=(a,0)` and initial second impact (rational parametrization)

$$
B(t)=\left(a\frac{1-\epsilon t^2}{1+\epsilon t^2},\;b\frac{2t}{1+\epsilon t^2}\right)
$$

when denominators are nonzero. The orbit's caustic parameter and the generalized Jacobi parameter are

$$
\lambda=\frac{\epsilon a^2b^2t^2}{b^2+a^2t^2},\qquad
m=\frac{\epsilon b^2-a^2}{b^2-\epsilon\lambda}.
$$

The **lifted** fixed-caustic state space is a generalized Jacobi intersection

$$
c^2+\epsilon s^2=1,\qquad d^2+m s^2=1.
$$

Its Weierstrass model (with the **same finite-field group cardinality**; in fact isomorphic via `x=(1+c)(1+d)/s^2`, `y=(1+c)(1+d)(c+d)/s^3`, see [elliptic_correspondence.md](elliptic_correspondence.md), section 3) is

$$
E_{\epsilon,m}:\quad y^2=x(x-\epsilon)(x-m).
$$

Both `m=0` and `m=epsilon` are singular cases and excluded. The model is **not** the original boundary conic, which has genus zero. The map from a reflection orbit into the lifted group requires tracking the extra branch `d`. For the original billiard, singular affine denominators must be handled via appropriate alternative charts; they are not resolved by blindly dividing.

The curve `E_{epsilon,m}` has full rational 2-torsion, so `4` divides `N=#E(F_p)`. We prioritize `N=4r` for a **prime** `r`, then verify that the actual reflection translation `G` has a nontrivial `r`-component (`[4]G` has order `r`). Counting `N` is *not* the same as proving the specific step has order `N`.

## Exact experimental orders

All values in the CSV/JSON are exact for the specified `p`, `epsilon` and `t`, not estimated from Hasse bounds. For `p <= 524287`, orders were computed with the quadratic-character point-count sum and factors checked. For 31- and 61-bit Mersenne primes and one Goldilocks-64 direction, orders were determined using elliptic-curve BSGS inside the Hasse interval; the possible order was made unique by combining certified orders of rational points with the known four rational 2-torsion points.

| p | ellipse `N` | ellipse largest prime factor | hyperbola `N` | hyperbola largest prime factor |
| ---: | ---: | ---: | ---: | ---: |
| 103 | 116 | 29 | 116 | 29 |
| 127 | 148 | 37 | 148 | 37 |
| 257 | 284 | 71 | 284 | 71 |
| 1009 | 1052 | 263 | 1052 | 263 |
| 1019 | 1076 | 269 | 1076 | 269 |
| 65521 | 65788 | 16447 | 65068 | 16267 |
| 65537 | 66032 | 4127 | 65644 | 16411 |
| 131071 | 131548 | 32887 | 131324 | 32831 |
| 524287 | 525244 | 131311 | 525284 | 131321 |
| 2147483647 | 2147539868 | 536884967 | 2147496172 | 536874043 |
| 2305843009213693951 | 2305843011732043748 | 576460752933010937 | 2305843009511141252 | 576460752377785313 |

The displayed candidates are selected *per conic* from an exhaustive sweep of the rational boundary slope `t` for primes through 1019 and deterministic finite samples for larger primes. They are the best found by this sweep, not a search over every starting state, every conic parameterization, or every prime-field curve.

At the Goldilocks prime `p=18446744069414584321`, **only a fixed `t=42` instance** is reported: ellipse `N=18446744067019653600`, max prime factor `754400599`; hyperbola `N=18446744066777287312`, max prime factor `10774967328725051`. Neither is a `4r` instance; other values of `t` at this prime have not been searched.

## Reproducible exact candidate at 61 bits

Using `p=2305843009213693951` and `a=5,b=3`:

- **Ellipse**: `t=1982817399921699654`, `lambda=189542885390157605`, `m=681520321761746604`; `N=2305843011732043748=4*576460752933010937`, actual reflected translation has order `1152921505866021874=2r`.
- **Hyperbola**: `t=1964008736653934424`, `lambda=1394708209538118479`, `m=1338682032702777952`; `N=2305843009511141252=4*576460752377785313`, actual reflected translation has order `576460752377785313=r`.

The first 64 reflected impacts were compared against the corresponding generalized Jacobi group multiples for each of the two 61-bit instances, with exact agreement.

## Larger cryptographic primes

Standard field primes, including `2^127-1`, `2^255-19` (Curve25519 field), `2^256-2^32-977` (secp256k1 field), and NIST P-256's prime, have **not** been assigned an exact order for this conic family by these computations. The **original named ECC curve's** subgroup order must never be copied to a new `E_{epsilon,m}` merely because the modulus is the same.

Use SageMath's SEA/other cardinality algorithms to compute each candidate's `N`, check `N=4r` with prime `r`, validate the actual billiard step, subgroup membership, exceptional points, pairing/embedding-degree risks, twists, and implementation side-channel properties. The included `tools/sage_large_prime_scan.py` performs the first point-counting and `4r` screening stages only. At a 255-bit field prime, `r` for a `4r` candidate can have at most about 253 bits, so its generic rho work factor is about `2^126.5` group operations, subject to constants.

## Large-prime samples (2026-11-18)

**Scope.** Original equal-angle mirror reflection on the finite-field axis-aligned ellipse

$$x^2/25+y^2/9=1\pmod p.$$

The group order is that of a **lifted fixed-caustic genus-one state curve**, not the genus-zero boundary conic. Four specific sample values of the initial rational parameter `t` were tested, one per field. All four samples passed: 64 direct mirror impacts against fast Jacobi-coordinate group multiplication, the documented curve-model parameter checks, and `[N]G=O`.

| Field | Exact group order N | Established prime subgroup | Status |
|---|---:|---:|---|
| Mersenne-127 | `170141183460469231757342577509242406420` | `8507059173023461587867128875462120321` (cofactor 20) | verified |
| Curve25519-field | `57896044618658097711785492504343953926856930875039260848015607506283634007912` | `7237005577332262213973186563042994240857116359379907606001950938285454250989` (cofactor 8) | verified |
| secp256k1-field | `115792089237316195423570985008687907852589613375810908478214554728305722413224` | Not yet established | N known; N/4 composite |
| NIST-P256-field | `115792089210356248762697446949407573530765394402118123589611579317542265087204` | Not yet established | N known; N/4 composite |

### Methods and certificates

- **127-bit Mersenne prime:** Selected a CM point-count example with discriminant −19. Its trace satisfies `4p = u² + 19v²`. The sign was uniquely fixed by an actual elliptic-curve point (one possible group order annihilates the point and the other does not). The candidate `r=N/20` is prime by a self-contained Lucas certificate generated from the complete factorization of `r−1`. The actual mirror translation `G` satisfies `[20]G != O` and `[r]([20]G)=O`.
- **Curve25519 prime field:** Started with the standard Montgomery curve `y²=x³+486662x²+x`, which has order `8r` in RFC 7748. Its 2-isogenous curve has equation `y²=x(x−486660)(x−486664)`. Because `486664` is a square in this field, it is isomorphic to `y²=x(x−1)(x−m)` with `m=486660/486664 (mod p)`. Therefore the lifted billiard curve has the same order. The actual billiard translation satisfies `[8]G != O` and `[r]([8]G)=O`, with `r` the standard Curve25519 252-bit prime.
- **secp256k1 and NIST P-256 *fields*:** These are *not* the named secp256k1 or NIST P-256 elliptic curves. The chosen billiard curves were selected with CM discriminants −7 and −67 respectively; exact orders were obtained from `4p = u²+Dv²`, then the Frobenius trace sign was confirmed against group-law point witnesses. In both sampled cases `N/4` is composite, so neither sample is established to have a 254-bit prime-order subgroup. Their largest prime factors have **not** been computed, and no conclusion about other parameter choices follows.

### Notes

1. Finite-field `j` values with special CM allow exact point counts without general SEA. The samples are therefore chosen by CM, not by an exhaustive search of caustics, curves, or prime-order subgroups.
2. Large prime subgroup security depends on implementation, twist/subgroup validation, embedding degree, endomorphisms, and cryptographic protocol design. The 127-bit example has a 123-bit prime subgroup, about `2^61` generic work, below current recommendations.
3. In the 255-bit example, security is inherited from a standard 2-isogenous elliptic-curve group, not from a new dynamical hardness assumption.
4. All tests use affine formula charts. The code does not implement constant-time or complete exceptional-case handling.
5. The **actual mirror translation** was checked (not merely the existence of a group of order N). This is critical when selecting a parameter for billiard fast-jump operations.

### Reproduction

```bash
python tools/verify_large_prime_samples.py
python -m unittest discover -s tests -v
```

For full unrestricted candidate searches at 255/256 bits, use `tools/sage_large_prime_scan.py` under SageMath. The samples above use special CM and published isogeny orders rather than generic SEA point counts.

References: RFC 7748 (Curve25519 parameters); the CM trace norm equation (e.g. MIT 18.783 Lecture 21); standard 2-isogeny formulas (Vélu).
