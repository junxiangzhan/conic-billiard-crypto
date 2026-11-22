# Large-prime conic mirror-screening results (2026-11-18)

**Scope.** Original equal-angle mirror reflection on the finite-field axis-aligned ellipse

\[x^2/25+y^2/9=1\pmod p.\]

The group order is that of a **lifted fixed-caustic genus-one state curve**, not the genus-zero boundary conic. Four specific sample values of the initial rational parameter `t` were tested, one per field. All four samples passed: 64 direct mirror impacts against fast Jacobi-coordinate group multiplication, the documented curve-model parameter checks, and `[N]G=O`.

| Field | Exact group order N | Established prime subgroup | Status |
|---|---:|---:|---|
| Mersenne-127 | `170141183460469231757342577509242406420` | `8507059173023461587867128875462120321` (cofactor 20) | verified |
| Curve25519-field | `57896044618658097711785492504343953926856930875039260848015607506283634007912` | `7237005577332262213973186563042994240857116359379907606001950938285454250989` (cofactor 8) | verified |
| secp256k1-field | `115792089237316195423570985008687907852589613375810908478214554728305722413224` | Not yet established | N known; N/4 composite |
| NIST-P256-field | `115792089210356248762697446949407573530765394402118123589611579317542265087204` | Not yet established | N known; N/4 composite |

## Methods and certificates

- **127-bit Mersenne prime:** Selected a CM point-count example with discriminant −19. Its trace satisfies `4p = u² + 19v²`. The sign was uniquely fixed by an actual elliptic-curve point (one possible group order annihilates the point and the other does not). The candidate `r=N/20` is prime by a self-contained Lucas certificate generated from the complete factorization of `r−1`. The actual mirror translation `G` satisfies `[20]G != O` and `[r]([20]G)=O`.
- **Curve25519 prime field:** Started with the standard Montgomery curve `y²=x³+486662x²+x`, which has order `8r` in RFC 7748. Its 2-isogenous curve has equation `y²=x(x−486660)(x−486664)`. Because `486664` is a square in this field, it is isomorphic to `y²=x(x−1)(x−m)` with `m=486660/486664 (mod p)`. Therefore the lifted billiard curve has the same order. The actual billiard translation satisfies `[8]G != O` and `[r]([8]G)=O`, with `r` the standard Curve25519 252-bit prime.
- **secp256k1 and NIST P-256 *fields*:** These are *not* the named secp256k1 or NIST P-256 elliptic curves. The chosen billiard curves were selected with CM discriminants −7 and −67 respectively; exact orders were obtained from `4p = u²+Dv²`, then the Frobenius trace sign was confirmed against group-law point witnesses. In both sampled cases `N/4` is composite, so neither sample is established to have a 254-bit prime-order subgroup. Their largest prime factors have **not** been computed, and no conclusion about other parameter choices follows.

## Mathematical caveats

1. Finite-field `j` values with special CM allow exact point counts without general SEA. This is **not** an exhaustive search of caustics, curves, or prime-order subgroups.
2. Large prime subgroup security depends on implementation, twist/subgroup validation, embedding degree, endomorphisms, and cryptographic protocol design. These are research demos only; the 127-bit example is far below modern security strength.
3. In the 255-bit example, security is inherited from a standard 2-isogenous elliptic-curve group, not from a new dynamical hardness assumption.
4. All tests use affine formula charts. The code does not implement constant-time or complete exceptional-case handling. Do not use these parameters in deployed key exchange.
5. The **actual mirror translation** was checked (not merely the existence of a group of order N). This is critical when selecting a parameter for billiard fast-jump operations.

## Reproduction

```bash
python tools/verify_large_prime_samples.py
python -m unittest discover -s tests -v
```

For full unrestricted candidate searches at 255/256 bits, use `tools/sage_large_prime_scan.py` under SageMath. Generic SEA point counts were **not** performed in the present environment; special CM and published isogeny orders were used instead.

References: RFC 7748 (Curve25519 parameters); the CM trace norm equation (e.g. MIT 18.783 Lecture 21); standard 2-isogeny formulas (Vélu).
