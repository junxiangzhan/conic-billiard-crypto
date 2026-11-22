# Finite-Field Conic Billiards: Group-Order Screening

This is a reproducible research screening report, **not cryptographically secure curve parameter generation**.

## Frozen conic mirror model

Let `p` be an odd prime, `epsilon = +1` (ellipse) or `-1` (hyperbola), and `a=5`, `b=3` in the experimental scans. The boundary is

\[
C_\epsilon:\quad (x/a)^2 + \epsilon(y/b)^2=1 \pmod p.
\]

Choose boundary start `A=(a,0)` and initial second impact (rational parametrization)

\[
B(t)=\left(a\frac{1-\epsilon t^2}{1+\epsilon t^2},\;b\frac{2t}{1+\epsilon t^2}\right)
\]

when denominators are nonzero. The orbit's caustic parameter and the generalized Jacobi parameter are

\[
\lambda=\frac{\epsilon a^2b^2t^2}{b^2+a^2t^2},\qquad
m=\frac{\epsilon b^2-a^2}{b^2-\epsilon\lambda}.
\]

The **lifted** fixed-caustic state space is a generalized Jacobi intersection

\[
c^2+\epsilon s^2=1,\qquad d^2+m s^2=1.
\]

Its isogenous Weierstrass model (with the **same finite-field group cardinality**) is

\[
E_{\epsilon,m}:\quad y^2=x(x-\epsilon)(x-m).
\]

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

The displayed candidates are selected *per conic* from an exhaustive sweep of the rational boundary slope `t` for primes through 1019 and deterministic finite samples for larger primes. They do **not** claim to be the optimal global choice or exhaustive over every starting state, every conic parameterization, or every prime-field curve.

At the Goldilocks prime `p=18446744069414584321`, **only a fixed `t=42` instance** is reported: ellipse `N=18446744067019653600`, max prime factor `754400599`; hyperbola `N=18446744066777287312`, max prime factor `10774967328725051`. Neither is a `4r` instance. This does **not** imply the prime cannot support a good instance; further search is needed.

## Reproducible exact candidate at 61 bits

Using `p=2305843009213693951` and `a=5,b=3`:

- **Ellipse**: `t=1982817399921699654`, `lambda=189542885390157605`, `m=681520321761746604`; `N=2305843011732043748=4*576460752933010937`, actual reflected translation has order `1152921505866021874=2r`.
- **Hyperbola**: `t=1964008736653934424`, `lambda=1394708209538118479`, `m=1338682032702777952`; `N=2305843009511141252=4*576460752377785313`, actual reflected translation has order `576460752377785313=r`.

The first 64 reflected impacts were compared against the corresponding generalized Jacobi group multiples for each of the two 61-bit instances, with exact agreement.

## Larger cryptographic primes: work still required

Standard field primes, including `2^127-1`, `2^255-19` (Curve25519 field), `2^256-2^32-977` (secp256k1 field), and NIST P-256's prime, have **not** been assigned an exact order for this conic family by these computations. The **original named ECC curve's** subgroup order must never be copied to a new `E_{epsilon,m}` merely because the modulus is the same.

Use SageMath's SEA/other cardinality algorithms to compute each candidate's `N`, check `N=4r` with prime `r`, validate the actual billiard step, subgroup membership, exceptional points, pairing/embedding-degree risks, twists, and implementation side-channel properties. The included `tools/sage_large_prime_scan.py` performs the first point-counting and `4r` screening stages only. At a 255-bit field prime, `r` for a `4r` candidate can have at most about 253 bits, so its generic rho work factor is about `2^126.5` group operations, subject to constants.

**This project is a mathematical experiment. Do not deploy it for real cryptographic key exchange.**
