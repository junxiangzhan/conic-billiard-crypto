# Compact Jacobi state format (experimental)

This encoding stores a **full affine lifted fixed-caustic state**, not merely the impact position `(x,y)`, in one field element and one branch bit. It is mathematically related to standard elliptic-curve point compression but not wire-compatible with X25519 or any other established ECC standard.

## Model

Fix odd prime `p`, epsilon `±1` (ellipse / generalized hyperbola), and nonsingular `m` not `0` or `epsilon`:

\[
 c^2+\epsilon s^2=1,\qquad d^2 + m s^2=1.
\]

Set `t=s/(1+c)` when `c != -1`. In the main rational chart:

\[
 s=\frac{2t}{1+\epsilon t^2},\quad
 c=\frac{1-\epsilon t^2}{1+\epsilon t^2},\quad
 w=d(1+\epsilon t^2),
\]

and

\[
 \boxed{w^2=(1+\epsilon t^2)^2-4mt^2=t^4+(2\epsilon-4m)t^2+1}.
\]

Therefore choose a canonical root `w` using one parity bit. Over an odd prime field, `w` and `-w` have opposite parity as canonical integers in `[0,p-1]` unless `w=0`; a zero root **must** encode with sign 0.

## Exact byte layout (version 1)

Let `k=p.bit_length()` and `L=ceil((k+1)/8)` bytes. Represent `t` using `k` low bits of a little-endian integer; bit `k` is a sign selector. Reject any nonzero padding bits above bit `k`.

- `0 <= t < p`: normal affine chart; sign is the parity of integer `w` in `[0,p-1]`.
- `t == p`: reserved **infinity chart** for `(s,c,d)=(0,-1,1)` with sign=0 and `(0,-1,-1)` with sign=1.
- `t > p`: invalid.
- Normal chart requires `1+epsilon*t² !=0`. Reject poles and quadratic nonresidues.

For `p=2^255-19`: `k=255`, `L=32 bytes`, high bit of the last little-endian byte is the branch bit. This is **not** the RFC 7748 X25519 high bit and must not be interpreted as such.

## Checkpoint: p=1009

Given elliptic model `m=635`, `U=(s,c,d)=(647,146,528)` on the Jacobi intersection, the compressed `t` field equals `176`; decoding reconstructs `(647,146,528)` and impact `(x,y)=(730,932)` when the public ellipse axes are `a=5,b=3`.

## Subgroup and security caveats

- An encoded point lies on a **fixed parameterized genus-one state curve**; decoding does **not** prove the point is in the DH subgroup. Verify `[r]P=O` with the correct prime order `r`, and reject the identity.
- The codec, affine `JacobiGroup` methods, and demo are variable-time and incompletely handle exceptional scalar-multiplication cases. This is research software, **not a cryptographic implementation**.
- This encoding has **public model parameters** `p,m,epsilon`; if they are not agreed in advance, transmitting them costs additional bytes.
- Hyperbolic encoding handles `epsilon=-1` algebraic Jacobi states, but **does not** automatically establish a reversible hyperbola mirror-orbit lift or a secure prime-order group.
- The `t==p` reserved code handles two excluded rational-chart points without ambiguity. Invalid padding, invalid square roots, and out-of-field payloads are rejected.
- This is an **affine chart codec**, not a universally complete encoding of the projective Jacobi intersection. In some finite fields the complete genus-one model has rational points at projective infinity that this version does **not** encode. For the chosen Curve25519-related example, the tested group elements lie in the covered affine chart. A fully general group encoder needs a separate treatment of all projective exceptional points.

## Try it

```bash
python -m conic_billiard.compressed_dh_demo
python -m unittest discover -s tests -v
```

The demo uses the known Curve25519-related group structure but **not** the X25519 protocol, its scalar clamping, or key derivation. It intentionally uses fixed demonstration secrets.
