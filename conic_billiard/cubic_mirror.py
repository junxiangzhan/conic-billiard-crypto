"""Finite-field *specular* billiards on a smooth Weierstrass cubic.

Not elliptic-curve chord-and-tangent group addition.  In characteristic != 2,3
this implementation uses the standard bilinear form dx^2+dy^2 and
H_B=I-2 N N^T/(N^T N) with N=grad(y^2-x^3-a*x-b).

A reflected line can have 0, 1, or 2 distinct rational next intersections,
and the algebraic relation is not a globally-defined single-valued map.
"""
from __future__ import annotations
from collections import Counter
from dataclasses import dataclass
from typing import Optional, Tuple, List

Point = Tuple[int, int]
State = Tuple[Point, Point]


def prime_check(n: int) -> bool:
    if n < 2:
        return False
    for q in range(2, int(n**0.5) + 1):
        if n % q == 0:
            return n == q
    return True


def square_roots_table(p: int) -> dict[int, tuple[int, ...]]:
    table: dict[int, list[int]] = {}
    for y in range(p):
        table.setdefault((y*y) % p, []).append(y)
    return {k: tuple(v) for k, v in table.items()}


@dataclass(frozen=True)
class HitResult:
    next_points: tuple[Point, ...]
    status: str
    reflected_direction: Optional[Point]
    polynomial: Optional[tuple[int, int, int]]  # c2*t^2+c1*t+c0
    discriminant: Optional[int]


class CubicMirror:
    def __init__(self, p: int = 101, a: int = -1, b: int = 0):
        if not prime_check(p) or p <= 3:
            raise ValueError('p must be a prime > 3 (educational small-field implementation)')
        self.p = p
        self.a = a % p
        self.b = b % p
        if (4*self.a**3+27*self.b**2) % p == 0:
            raise ValueError('singular cubic')
        self.squares = square_roots_table(p)
        self.points = tuple((x, y) for x in range(p)
                            for y in self.squares.get((x*x*x+self.a*x+self.b) % p, ()))
        self.pointset = set(self.points)

    def valid(self, P: Point) -> bool:
        x,y = P
        return (y*y-x*x*x-self.a*x-self.b) % self.p == 0

    def legendre(self, v: int) -> int:
        v %= self.p
        if v == 0:
            return 0
        return 1 if pow(v,(self.p-1)//2,self.p)==1 else -1

    def chord_norm_class(self, A: Point, B: Point) -> int:
        x=(B[0]-A[0])%self.p; y=(B[1]-A[1])%self.p
        return self.legendre(x*x+y*y)

    def successors(self, A: Point, B: Point) -> HitResult:
        """All DISTINCT affine next hits C != B, including possible C==A.

        If the reflected line is vertical, its third intersection is the
        Weierstrass point O at infinity; we return only its finite root.
        When a root is t=0 (contact at B), this is a tangent degeneracy,
        not a positive-length next hit. Finite fields have no ordering on t.
        """
        p=self.p
        if A not in self.pointset or B not in self.pointset or A==B:
            return HitResult((), 'invalid-state', None, None, None)
        bx,by=B
        nx=(-3*bx*bx-self.a)%p; ny=(2*by)%p
        nn=(nx*nx+ny*ny)%p
        if nn==0:
            return HitResult((), 'isotropic-normal', None, None, None)
        vx=(bx-A[0])%p; vy=(by-A[1])%p
        fac=(2*(vx*nx+vy*ny)*pow(nn,-1,p))%p
        u=(vx-fac*nx)%p; v=(vy-fac*ny)%p
        # F(B+t*w) = t*(c2*t^2+c1*t+c0)
        c2=(-u*u*u)%p
        c1=(v*v-3*bx*u*u)%p
        c0=(2*by*v-(3*bx*bx+self.a)*u)%p
        roots: list[int]
        discr=None
        if c2:
            discr=(c1*c1-4*c2*c0)%p
            roots=[]
            denom=pow((2*c2)%p,-1,p)
            for sr in self.squares.get(discr,()):
                roots.append(((-c1+sr)*denom)%p)
        elif c1:
            roots=[(-c0*pow(c1,-1,p))%p]
        else:
            # Should not be possible for nonsingular cubic with w != 0
            return HitResult((), 'degenerate-intersection', (u,v), (c2,c1,c0), None)
        hitzero=0 in roots
        nextpoints=tuple(sorted({((bx+t*u)%p,(by+t*v)%p) for t in roots if t != 0}))
        assert all(self.valid(C) and C!=B for C in nextpoints)
        status = ('tangent-at-current' if hitzero else
                  'no-affine-next' if not nextpoints else
                  'one-next' if len(nextpoints)==1 else 'two-next')
        return HitResult(nextpoints,status,(u,v),(c2,c1,c0),discr)

    def state_graph_counts(self, sample_states: Optional[int]=None, seed: int=20261118):
        """For sample_states=None examine all ordered pairs A!=B."""
        import random
        if sample_states is None:
            pairs=((A,B) for A in self.points for B in self.points if A!=B)
        else:
            rng=random.Random(seed)
            pairs=[]
            while len(pairs)<sample_states:
                A,B = rng.sample(self.points, 2)
                pairs.append((A,B))
        statuses=Counter()
        degrees=Counter()
        reversals=0
        normclass_ok=0
        normal_edges=0
        witnesses={}
        total=0
        for A,B in pairs:
            total+=1
            result=self.successors(A,B)
            statuses[result.status]+=1
            degrees[len(result.next_points)]+=1
            witnesses.setdefault(result.status,(A,B,result.next_points))
            for C in result.next_points:
                normal_edges+=1
                reverse=self.successors(C,B)
                if A in reverse.next_points: reversals+=1
                if self.chord_norm_class(A,B) == self.chord_norm_class(B,C):
                    normclass_ok+=1
        return {'p':self.p,'a':self.a,'b':self.b,'affine_points':len(self.points),
                'states_tested':total, 'status_counts':dict(statuses),
                'outdegree_counts':dict(degrees), 'edges_tested':normal_edges,
                'reverse_edges_verified':reversals,'chord_norm_squareclass_verified':normclass_ok,
                'witnesses':witnesses}

    def breadth_reachable(self, A: Point, B: Point, steps: int = 10):
        states={(A,B)}
        counts=[1]
        for _ in range(steps):
            next_states=set()
            for X,Y in states:
                for Z in self.successors(X,Y).next_points:
                    next_states.add((Y,Z))
            states=next_states
            counts.append(len(states))
        return counts

    def deterministic_cycle_scan(self, policy: str='min'):
        """Arbitrary finite-field choice of a successor; not geometrically canonical.
        Undefined/terminal states remain terminal (not forcibly connected).
        Returns maximum eventual cycle length for all non-degenerate states.
        """
        if policy not in ('min','max'):
            raise ValueError(policy)
        states=[(A,B) for A in self.points for B in self.points if A!=B]
        nxt={}
        for st in states:
            opts=self.successors(*st).next_points
            nxt[st] = (st[1], min(opts) if policy=='min' else max(opts)) if opts else None
        seen={}
        cycles=Counter()
        examples={}
        terminal=0
        for seed in states:
            if seed in seen: continue
            local={}
            path=[]
            cur=seed
            while cur is not None and cur not in seen and cur not in local:
                local[cur]=len(path)
                path.append(cur)
                cur=nxt[cur]
            if cur is None:
                terminal+=1
            elif cur in local:
                period=len(path)-local[cur]
                cycles[period]+=1
                examples.setdefault(period,cur)
            for st in path:
                seen[st]=True
        return {'policy':policy, 'states':len(states), 'cycles_by_period':dict(sorted(cycles.items())),
                'max_cycle':max(cycles,default=0),
                'max_cycle_witness':examples[max(cycles)] if cycles else None,
                'ending_paths':terminal}


def monomials_upto(deg: int, variables: int=4):
    """Exponent tuples, total degree <= deg."""
    from itertools import product
    return [e for e in product(range(deg+1), repeat=variables) if sum(e)<=deg]


def polynomial_invariant_scan(mirror: CubicMirror, max_degree: int=4,
                              max_edges: int=220):
    """Solve I(A,B)=I(B,C) for low-total-degree polynomial first integrals.

    This uses sampled transition relations to compute a lower bound on rank,
    hence an UPPER bound on the space of *global* polynomial integrals.
    If this bound matches the dimension of obvious curve equations, then
    no other low-degree polynomial integrals exist in this model.
    Not a test for rational, algebraic, or larger-degree invariants.
    """
    import random
    rng=random.Random(20261118)
    p=mirror.p
    edges=[]
    for _ in range(max_edges*20):
        A,B=rng.sample(mirror.points,2)
        opts=mirror.successors(A,B).next_points
        if opts:
            edges.append((A,B,rng.choice(opts)))
        if len(edges)==max_edges:
            break
    if len(edges)<max_edges:
        raise ValueError('Not enough edges to test this sample size')
    results={}
    for degree in range(1,max_degree+1):
        exps=monomials_upto(degree)
        dim=len(exps)
        pivotrows={}
        for A,B,C in edges:
            uv=A+B; vw=B+C
            row=[(pow(vw[0],e[0],p)*pow(vw[1],e[1],p)*pow(vw[2],e[2],p)*pow(vw[3],e[3],p)-
                  pow(uv[0],e[0],p)*pow(uv[1],e[1],p)*pow(uv[2],e[2],p)*pow(uv[3],e[3],p))%p
                 for e in exps]
            for j,old in sorted(pivotrows.items()):
                if row[j]:
                    factor=row[j]
                    row=[(v-factor*vv)%p for v,vv in zip(row,old)]
            pivot=next((i for i,v in enumerate(row) if v), None)
            if pivot is not None:
                inv=pow(row[pivot],-1,p)
                pivotrows[pivot]=[(v*inv)%p for v in row]
        trivial=1 if degree<3 else 1 + 2*(1 if degree==3 else 5)
        results[degree]={
            'monomials':dim, 'rank_on_sample':len(pivotrows),
            'nullity_upper_bound':dim-len(pivotrows),
            'trivial_constant_and_curve_relations_dimension':trivial,
            'no_nontrivial_integral_in_this_polynomial_class':dim-len(pivotrows)==trivial,
        }
    return {'edges_used':len(edges), 'degree_results':results}
