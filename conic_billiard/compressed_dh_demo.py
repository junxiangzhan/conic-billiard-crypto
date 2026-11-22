"""32-byte public-key encoding demo on a 255-bit elliptic billiard group.

Reuses the standard Curve25519-related prime subgroup, but NOT the X25519
protocol or wire format. No real secret keys or KDF; demonstration only.

Run: python -m conic_billiard.compressed_dh_demo
"""
from __future__ import annotations
from .core import Conic
from .jacobi import ellipse_lift
from .compression import JacobiCodec

P = 2**255 - 19
R = 2**252 + 27742317777372353535851937790883648493
M = 20800338683988658368647408995589388737092878452977063003340006470870624536394
T = 4604094445982880512620096533988940857801848594535175639138443919056596161201


def verify_public(group, point, order=R):
    """Demo-only subgroup check. Uses variable-time affine point addition."""
    if point == group.identity():
        raise ValueError('identity public key is forbidden')
    if group.mul(order, point) != group.identity():
        raise ValueError('public key is not in expected subgroup')
    return point


def demo() -> dict:
    a,b=5,3
    den=pow(1+T*T, -1, P)
    B=(a*(1-T*T)*den%P, 2*b*T*den%P)
    _, group, step = ellipse_lift(Conic.ellipse(P,a,b),a,b,(a,0),B)
    assert group.m == M
    codec=JacobiCodec(P,M)
    base=group.mul(8,step)
    assert base != group.identity() and group.mul(R,base)==group.identity()

    alice_secret, bob_secret=123456789, 987654321 # deterministic examples, NOT secrets
    A=group.mul(alice_secret,base)
    B=group.mul(bob_secret,base)
    Alice_wire=codec.encode(A)
    Bob_wire=codec.encode(B)

    # Parsing is separate from public-key subgroup validation.
    Alice_public=verify_public(group,codec.decode(Alice_wire))
    Bob_public=verify_public(group,codec.decode(Bob_wire))
    shared_a=group.mul(alice_secret,Bob_public)
    shared_b=group.mul(bob_secret,Alice_public)
    assert shared_a==shared_b
    return {
        'field_bits':P.bit_length(),
        'subgroup_bits':R.bit_length(),
        'encoded_public_key_bytes':len(Alice_wire),
        'alice_public_hex':Alice_wire.hex(),
        'bob_public_hex':Bob_wire.hex(),
        'shared_point_agrees':shared_a==shared_b,
        'warning':'Toy, variable-time, unauthenticated; do not use for actual key exchange',
    }


if __name__=='__main__':
    import json
    print(json.dumps(demo(),indent=2))
