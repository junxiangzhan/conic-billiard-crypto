"""Canonical compressed encoding for generalized Jacobi intersection points.

Research/education only: variable-time Tonelli-Shanks and affine coordinates.
Not suitable for production cryptography or direct use with untrusted keys.

Model over odd prime F_p:
    c^2 + epsilon*s^2 = 1
    d^2 + m*s^2 = 1
where epsilon in {+1, -1} and m not in {0, epsilon}.

One field element t=s/(1+c) plus one bit selecting w=d(1+epsilon*t^2)
recovers (s,c,d), except c=-1, which has a reserved value t=p.

Packed little-endian integer: payload=t in low p.bit_length() bits;
root parity in bit p.bit_length(). Total size ceil((bits+1)/8).

For p=2^255-19 this is exactly 32 bytes. This is *not* the RFC 7748
encoding and is not interoperable with X25519.
"""
from __future__ import annotations
from dataclasses import dataclass

JacobiPoint = tuple[int, int, int]


class InvalidEncoding(ValueError):
    """Non-canonical or invalid Jacobi point encoding."""


def sqrt_mod(a: int, p: int) -> int | None:
    """Tonelli-Shanks: one square root modulo odd prime p, if it exists.

    Returns None for quadratic nonresidues. Variable-time; no side-channel
    resistance, primality proof, or trusted-input guarantee.
    """
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q = p - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    c = pow(z, q, p)
    x = pow(a, (q + 1) // 2, p)
    t = pow(a, q, p)
    m = s
    while t != 1:
        i = 1
        t2 = t * t % p
        while i < m and t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        if i == m:
            raise ArithmeticError("Tonelli-Shanks internal error")
        b = pow(c, 1 << (m - i - 1), p)
        x = x * b % p
        t = t * b * b % p
        c = b * b % p
        m = i
    return x


@dataclass(frozen=True)
class JacobiCodec:
    """Serialize full *affine lifted* Jacobi states, not positions alone.

    For epsilon=-1, the hyperbolic Jacobi lift is handled algebraically;
    constructing it from a hyperbola reflection orbit is a separate step.
    Rational points at projective infinity, when present in the complete
    intersection, are NOT encoded by this affine codec.

    Decoding guarantees membership on the intersection only, NOT subgroup
    membership. Use external validated group/subgroup logic before any DH.
    """
    p: int
    m: int
    epsilon: int = 1

    def __post_init__(self) -> None:
        from .core import is_probable_prime
        if self.p < 3 or self.p % 2 != 1 or not is_probable_prime(self.p):
            raise ValueError("p must be an odd prime")
        if self.epsilon not in (1, -1):
            raise ValueError("epsilon must be +1 or -1")
        m = self.m % self.p
        if m == 0 or m == self.epsilon % self.p:
            raise ValueError("singular generalized Jacobi intersection")
        object.__setattr__(self, 'm', m)

    @property
    def field_bits(self) -> int:
        return self.p.bit_length()

    @property
    def byte_length(self) -> int:
        return (self.field_bits + 1 + 7) // 8

    def valid_point(self, U: JacobiPoint) -> bool:
        s, c, d = U
        p = self.p
        if any(not isinstance(v, int) or not (0 <= v < p) for v in U):
            return False
        return ((c*c + self.epsilon*s*s - 1) % p == 0
                and (d*d + self.m*s*s - 1) % p == 0)

    def encode(self, U: JacobiPoint) -> bytes:
        if not self.valid_point(U):
            raise ValueError("not a canonical on-curve Jacobi state")
        s,c,d = U
        p = self.p
        if (1+c) % p == 0:
            # Necessarily s=0 and d=+1 or -1; c=-1 chart at infinity.
            assert s == 0 and d in (1, p-1)
            t = p  # reserved (out-of-field) payload for infinity chart
            sign = 0 if d == 1 else 1
        else:
            t = s * pow((1+c) % p, -1, p) % p
            den = (1+self.epsilon*t*t) % p
            if den == 0:
                raise ValueError("point maps to invalid rational chart")
            w = d * den % p
            sign = w & 1
        code = t | (sign << self.field_bits)
        return code.to_bytes(self.byte_length, 'little')

    def decode(self, data: bytes) -> JacobiPoint:
        if not isinstance(data, (bytes, bytearray)) or len(data) != self.byte_length:
            raise InvalidEncoding(f"expected exactly {self.byte_length} bytes")
        p = self.p
        bits = self.field_bits
        value = int.from_bytes(data, 'little')
        if value >> (bits + 1):
            raise InvalidEncoding("nonzero padding bits")
        sign = (value >> bits) & 1
        t = value & ((1 << bits) - 1)
        if t > p:
            raise InvalidEncoding("payload outside field or reserved infinity")
        if t == p:
            # one special chart value + sign 0/1 for (0,-1, +1/-1)
            U = (0, p-1, p-1 if sign else 1)
            assert self.valid_point(U)
            return U
        den = (1 + self.epsilon * t*t) % p
        if den == 0:
            raise InvalidEncoding("rational chart pole")
        # square root of the quartic.  Epsilon^2 == 1.
        rhs = ((1 + self.epsilon*t*t)**2 - 4*self.m*t*t) % p
        w = sqrt_mod(rhs, p)
        if w is None:
            raise InvalidEncoding("quartic value is not a quadratic residue")
        if (w & 1) != sign:
            w = (-w) % p
        if (w & 1) != sign:
            # w=0 can only be encoded with sign=0.
            raise InvalidEncoding("noncanonical sign of a zero root")
        inv = pow(den, -1, p)
        s = 2*t*inv % p
        c = (1 - self.epsilon*t*t)*inv % p
        d = w*inv % p
        U = s,c,d
        if not self.valid_point(U):
            raise InvalidEncoding("reconstructed state not on Jacobi intersection")
        return U

    def encode_hex(self, U: JacobiPoint) -> str:
        return self.encode(U).hex()

    def decode_hex(self, value: str) -> JacobiPoint:
        try:
            data = bytes.fromhex(value)
        except (TypeError, ValueError) as exc:
            raise InvalidEncoding("invalid hex string") from exc
        return self.decode(data)
