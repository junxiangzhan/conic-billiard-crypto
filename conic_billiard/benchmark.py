"""Reproducible functional comparison, not a cross-platform speed claim."""
from __future__ import annotations
from time import perf_counter
from .core import Conic, matvec4
from .jacobi import ellipse_lift, project_ellipse
from .search import analyze


def elapsed(fn, repeats=3):
    best = float('inf')
    value = None
    for _ in range(repeats):
        start = perf_counter()
        value = fn()
        best = min(best,perf_counter()-start)
    return value,best*1000


def matrix_trace(curve, A, B, n):
    state = A,B
    for _ in range(n):
        M = curve.transfer_matrix(*state)
        state = matvec4(M,state,curve.p)
    return state


def main():
    p=1009; a=5;b=3
    curve=Conic.ellipse(p,a,b)
    A,B=(a,0),(730,932)
    lam,J,G=ellipse_lift(curve,a,b,A,B)
    n=10000
    simple,ts = elapsed(lambda:curve.trace(A,B,n),repeats=2)
    matrix,tm = elapsed(lambda:matrix_trace(curve,A,B,n),repeats=2)
    lifted,tf = elapsed(lambda:(project_ellipse(J.mul(n,G),a,b,p),
                                 project_ellipse(J.mul(n+1,G),a,b,p)),repeats=7)
    assert simple == matrix == lifted
    print('Ellipse F_1009, n=',n, 'lambda=',lam)
    print('Sequential reflection :',round(ts,3),'ms')
    print('4x4 state matrices    :',round(tm,3),'ms')
    print('Jacobi group fastjump :',round(tf,3),'ms')
    print('All return SAME complete consecutive-point state:',simple)
    print('Warning: Jacobi method uses known group law, NOT raw 4x4 matrix doubling.')

    h=Conic.hyperbola(103,5,3)
    stats=analyze(h)
    state=stats.max_cycle_initial_state
    assert state is not None
    length=stats.max_period
    assert h.trace(*state,length)==state
    print('Hyperbola F_103: period',length,'verified by direct reflection and reversible 4x4 matrices')
    assert matrix_trace(h,*state,length)==state
    print('Warning: NO hyperbola fastjump/group lifting is claimed in this prototype.')


if __name__=='__main__':
    main()
