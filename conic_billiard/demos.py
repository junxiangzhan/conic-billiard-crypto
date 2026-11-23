"""Reproducible demonstrations: two DH demos and a speed comparison.

All exchanges use the Jacobi lift of an elliptic billiard orbit, i.e. a
known elliptic-curve group law. Variable-time, unauthenticated, no KDF:
NOT for real keys. The hyperbola lift (c^2 - s^2 = 1, d^2 + m s^2 = 1) is
proved and implemented in tools/ (poncelet_pencil.py, dead_ends.py) but has
no DH demo yet; here the hyperbola only gets a small-field period check.

Run via the CLI: python -m conic_billiard dh | bench
"""
from __future__ import annotations
from time import perf_counter
from .core import Conic, matvec4
from .compression import JacobiCodec
from .jacobi import ellipse_lift, project_ellipse, bsgs
from .search import analyze

# Curve25519-related prime subgroup (NOT the X25519 protocol or wire format).
P = 2**255 - 19
R = 2**252 + 27742317777372353535851937790883648493
M = 20800338683988658368647408995589388737092878452977063003340006470870624536394
T = 4604094445982880512620096533988940857801848594535175639138443919056596161201


def small_field_dh() -> dict:
    """F_1009, subgroup of order 67: an exchange that BSGS breaks instantly."""
    p = 1009
    curve = Conic.ellipse(p, 5, 3)
    A, B = (5, 0), (730, 932)
    lam, group, G = ellipse_lift(curve, 5, 3, A, B)
    H = group.mul(2, G)
    assert lam == 964 and group.m == 635
    assert group.mul(134, G) == group.identity()
    assert group.mul(67, H) == group.identity() and H != group.identity()

    # The lifted group orbit matches the original reflection steps.
    state = (A, B)
    for n in range(134):
        assert project_ellipse(group.mul(n, G), 5, 3, p) == state[0]
        state = curve.step(state)
    assert state == (A, B)

    a, b = 17, 29
    public_A, public_B = group.mul(a, H), group.mul(b, H)
    shared = group.mul(a, public_B)
    assert shared == group.mul(b, public_A) == group.mul(a * b, H)
    recovered_a = bsgs(group, H, public_A, 67)
    assert recovered_a == a
    # A naive additive exchange is computable by any observer.
    assert group.add(public_A, public_B) == group.mul(a + b, H)
    return {
        'p': p, 'A': A, 'B': B, 'lambda': lam, 'm': group.m,
        'base_order': 134, 'subgroup_order': 67,
        'Alice_public': public_A, 'Bob_public': public_B,
        'shared_Jacobi': shared,
        'shared_collision_point': project_ellipse(shared, 5, 3, p),
        'BSGS_recovers_Alice_secret': recovered_a,
    }


def verify_public(group, point, order=R):
    """Demo-only subgroup check. Uses variable-time affine point addition."""
    if point == group.identity():
        raise ValueError('identity public key is forbidden')
    if group.mul(order, point) != group.identity():
        raise ValueError('public key is not in expected subgroup')
    return point


def compressed_dh() -> dict:
    """255-bit field, 32-byte encoded public keys, subgroup-validated."""
    a, b = 5, 3
    den = pow(1 + T * T, -1, P)
    B = (a * (1 - T * T) * den % P, 2 * b * T * den % P)
    _, group, step = ellipse_lift(Conic.ellipse(P, a, b), a, b, (a, 0), B)
    assert group.m == M
    codec = JacobiCodec(P, M)
    base = group.mul(8, step)
    assert base != group.identity() and group.mul(R, base) == group.identity()

    alice_secret, bob_secret = 123456789, 987654321  # fixed examples, NOT secrets
    alice_wire = codec.encode(group.mul(alice_secret, base))
    bob_wire = codec.encode(group.mul(bob_secret, base))
    # Parsing is separate from public-key subgroup validation.
    alice_public = verify_public(group, codec.decode(alice_wire))
    bob_public = verify_public(group, codec.decode(bob_wire))
    agrees = group.mul(alice_secret, bob_public) == group.mul(bob_secret, alice_public)
    assert agrees
    return {
        'field_bits': P.bit_length(),
        'subgroup_bits': R.bit_length(),
        'encoded_public_key_bytes': len(alice_wire),
        'alice_public_hex': alice_wire.hex(),
        'bob_public_hex': bob_wire.hex(),
        'shared_point_agrees': agrees,
    }


def _best_ms(fn, repeats):
    best, value = float('inf'), None
    for _ in range(repeats):
        start = perf_counter()
        value = fn()
        best = min(best, perf_counter() - start)
    return value, best * 1000


def matrix_trace(curve, A, B, n):
    state = A, B
    for _ in range(n):
        state = matvec4(curve.transfer_matrix(*state), state, curve.p)
    return state


def benchmark() -> None:
    """Sequential reflection vs 4x4 matrices vs Jacobi fast jump (same result)."""
    p, a, b, n = 1009, 5, 3, 10000
    curve = Conic.ellipse(p, a, b)
    A, B = (a, 0), (730, 932)
    lam, J, G = ellipse_lift(curve, a, b, A, B)
    simple, ts = _best_ms(lambda: curve.trace(A, B, n), 2)
    matrix, tm = _best_ms(lambda: matrix_trace(curve, A, B, n), 2)
    lifted, tf = _best_ms(lambda: (project_ellipse(J.mul(n, G), a, b, p),
                                   project_ellipse(J.mul(n + 1, G), a, b, p)), 7)
    assert simple == matrix == lifted
    print(f'Ellipse F_{p}, n={n}, lambda={lam}')
    print(f'Sequential reflection : {ts:.3f} ms')
    print(f'4x4 state matrices    : {tm:.3f} ms')
    print(f'Jacobi group fastjump : {tf:.3f} ms  (known group law, NOT matrix doubling)')
    print('All return the same state:', simple)

    h = Conic.hyperbola(103, 5, 3)
    stats = analyze(h)
    state, length = stats.max_cycle_initial_state, stats.max_period
    assert h.trace(*state, length) == state == matrix_trace(h, *state, length)
    print(f'Hyperbola F_103: period {length} verified by reflection and 4x4 matrices')
