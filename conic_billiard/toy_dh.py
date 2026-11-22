"""ECDH-like exchange on a Jacobi lift of one elliptic billiard orbit.

Small field and subgroup: intentionally breakable. NOT for real keys.
The corresponding fast-jump implementation for *hyperbola* is not supplied:
that is a separate research task needing proof and validation.
"""
from __future__ import annotations
from .core import Conic
from .jacobi import ellipse_lift, project_ellipse, bsgs


def demo() -> dict:
    p = 1009
    curve = Conic.ellipse(p,5,3)
    A, B = (5,0), (730,932)
    lam, group, G = ellipse_lift(curve, 5,3,A,B)
    H = group.mul(2,G)
    assert lam == 964 and group.m == 635
    assert group.mul(134,G) == group.identity()
    assert group.mul(67,H) == group.identity() and H != group.identity()

    # Ensure actual original reflection steps match the lifted group orbit.
    state = (A,B)
    for n in range(134):
        assert project_ellipse(group.mul(n,G),5,3,p) == state[0]
        state = curve.step(state)
    assert state == (A,B)

    a, b = 17,29
    public_A = group.mul(a,H)
    public_B = group.mul(b,H)
    shared_A = group.mul(a,public_B)
    shared_B = group.mul(b,public_A)
    assert shared_A == shared_B == group.mul(a*b,H)
    recovered_a = bsgs(group,H,public_A,67)
    assert recovered_a == a
    assert group.mul(recovered_a,public_B) == shared_A
    # Incorrect additive exchange is transparently computable by observers.
    naive_shared = group.mul(a+b,H)
    assert group.add(public_A, public_B) == naive_shared
    return {
        'curve': 'ellipse', 'p':p, 'A':A, 'B':B,
        'lambda':lam, 'm':group.m, 'base_order':134, 'subgroup_order':67,
        'Alice_public':public_A, 'Bob_public':public_B,
        'shared_Jacobi':shared_A,
        'shared_collision_point':project_ellipse(shared_A,5,3,p),
        'BSGS_recovers_Alice_secret':recovered_a,
        'note':'toy: ECDH-type group construction, NOT a secure protocol',
    }


if __name__ == '__main__':
    import json
    print(json.dumps(demo(), indent=2, ensure_ascii=False))
