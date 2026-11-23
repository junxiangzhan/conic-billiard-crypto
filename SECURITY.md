# Security scope and threat model

This repository is a research implementation of dynamical systems and cryptography. It is not a DH protocol, and it has no constant-time implementation or protocol-level security engineering.

## What is protected

Public: the field $\mathbb F_p$, the boundary conic, the reflection rule, the starting point $(a,0)$ and the second impact $P_1$. Secret: the number of bounces $n$ (modulo the orbit period).

## Hardness is the ECDLP, in both directions

[elliptic_correspondence.md](docs/elliptic_correspondence.md) proves two things:

- **Billiard → curve:** the reflection orbit is $P_n=\pi([n]G)$, where $G$ is a point on the genus-1 curve $J\cong\{y^2=x(x-\epsilon)(x-m)\}$. So counting bounces is no harder than the ECDLP on $J$ with base $G$.
- **Curve → billiard:** conversely, for a point on a curve with full rational 2-torsion that satisfies a square condition, there is a billiard whose bounce-counting problem is that ECDLP. Details and limits of the reduction are in section 4 of that document.

Consequently:

1. **The billiard adds no new hardness assumption.** Its security is exactly the ECDLP of the associated curve; good parameters reach at most the Pollard rho bound $\sqrt r$.
2. **The group order is always a multiple of 4.** All three points of order 2 are rational, so the best case is $N=4r$ with $r$ prime. Standard prime-order curves (secp256k1, P-256) cannot be realized in this model.
3. **Parameter selection follows standard ECC practice:** check for a large prime subgroup, the embedding degree, twist security and so on; see [group_order_search.md](docs/group_order_search.md). Since $4\mid N$, the anomalous case ($N=p$) cannot occur.
4. **All of these conditions have billiard meanings** ([security_dictionary.md](docs/security_dictionary.md)): circles, chords through a focus and parabolas degenerate to multiplicative groups; the group order is the number of states sharing a caustic, $2\cdot\#\text{outside}+\#\text{on}$, while the twist counts the impacts inside; the curve is supersingular, and open to MOV, exactly when as many impacts lie inside the caustic as outside. Embedding degree and CM have only arithmetic meaning and are listed in the comparison table of that document.

## No fixed small-matrix shortcut

One reflection step is a $4\times4$ matrix $M(X)$, but it depends on the state. A natural question: if some fixed matrix $M$ satisfied $X_n=M^nX_0$, would the discrete logarithm drop to a finite field and become weaker?

For good parameters, no such matrix exists. Suppose $\langle G\rangle$ contains an element of prime order $r\neq p$. If it is faithfully represented by a matrix in $GL_d(\mathbb F_p)$, the eigenvalues of that matrix lie in some $\mathbb F_{p^j}$ with $j\le d$, so $r\mid p^j-1$; that is, $d$ is at least the embedding degree $k$.

- Small example, $p=1009$, $r=67$: $k=33$, so at least a $33\times33$ matrix is needed.
- The Curve25519 subgroup: the order of $p$ modulo $\ell$ is exactly $k=(\ell-1)/6$, about $2^{249.4}$. This is computed from the complete factorization $\ell-1=2^2\cdot3\cdot11\cdot q_1\cdot q_2$ (both large factors verified prime with sympy) and agrees with SafeCurves.

So "a fixed small matrix exists" **implies** "the embedding degree is small", and such curves are already broken by MOV. The parabola can be written as a fixed $2\times2$ matrix, which is why that control case is weak.

## Consecutive bounces do not give DH

$$
T^a(T^b(X_0))=T^{a+b}(X_0)
$$

This adds the secrets, while DH needs the product $[a]\,([b]G)=[ab]G$. If $U_A=U_0+[a]G$ and $U_B=U_0+[b]G$ were published and $U_0+[a+b]G$ used as the shared value, an observer would compute it as $U_A+U_B-U_0$ without solving any DLP.

Multiplication needs "add an arbitrary state". There are two implementations:

- **Geometric ladder** ([geometric_group_law.md](docs/geometric_group_law.md)): computed directly in the billiard plane, proved correct, no square roots needed.
- **Jacobi formulas**: used by `python -m conic_billiard dh`, which **immediately recovers the demo private key** with BSGS. The demo subgroup has order 67.

Both compute the same group, so the choice does not affect security.

## Protocol-level concerns

- **Direction bit:** publishing only the impact point loses the direction (the choice of $d$ or $-d$). Without it, the receiver gets the mirrored impact when $a$ is even. When $m-\epsilon$ is a non-square and the peer's point lies in an odd-order subgroup, the direction can be recovered from the conic type; otherwise it must be sent with the public key. See [compressed_points.md](docs/compressed_points.md).
- **Ladder degeneracy:** the ladder degenerates and raises an exception for some $n$; reduce $n$ below the group order first.
- **Subgroup validation:** received points must be checked to lie in the correct subgroup. In billiard terms this is three checks: the state's caustic is the agreed one (otherwise an invalid-curve attack); the impact lies outside the caustic (impacts inside belong to the twist, and compressed decoding rejects them); the orbit is not a short period (small subgroup). See Theorem 5 of [security_dictionary.md](docs/security_dictionary.md).

## Not provided

Key generation, secure group-order selection, trusted point validation, KDF, authentication, CCA-secure encryption, error handling, side-channel protection and a secure random source are not provided. **Do not use this code for real confidential communication.**

## Implementation limits and reporting

Modular square roots, affine group operations and the geometric ladder are not constant-time; some projective exceptional points are handled only in the research tools; there has been no formal security analysis, protocol certification or side-channel testing. Some examples use deliberately small groups.

Do not post secret data in public issues. For research discussion, share only public parameters and reproducible, non-secret test vectors.
