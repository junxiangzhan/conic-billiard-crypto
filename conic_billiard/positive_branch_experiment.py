"""Canonical 'positive' branch conventions for cubic mirror billiards over F_p.

This investigates branch selection; it does NOT produce a secure cryptosystem.
"""
from __future__ import annotations
from collections import Counter
from math import isqrt
from .cubic_mirror import CubicMirror


def options_with_t(m, A, B):
    h=m.successors(A,B)
    if h.reflected_direction is None or h.polynomial is None:
        return (),h
    c2,c1,c0=h.polynomial
    p=m.p
    if c2:
        roots=m.squares.get(h.discriminant, ())
        ts={((-c1+s)*pow(2*c2,-1,p))%p for s in roots}
    elif c1:
        ts={(-c0*pow(c1,-1,p))%p}
    else:
        ts=set()
    u,v=h.reflected_direction
    return tuple(sorted((t, ((B[0]+t*u)%p,(B[1]+t*v)%p)) for t in ts if t)),h


def choose(m, A, B, policy):
    opts,h=options_with_t(m,A,B)
    if not opts: return None
    p=m.p
    if policy=='positive_t':
        allowed=[(t,C) for t,C in opts if 1<=t<=p//2]
        return min(allowed)[1] if allowed else None
    if policy=='min_point':
        return min(C for t,C in opts)
    if policy=='plus_sqrt':
        c2,c1,c0=h.polynomial
        if not c2:
            return opts[0][1]
        sr=min(m.squares[h.discriminant])
        t=(-c1+sr)*pow(2*c2,-1,p)%p
        return next((C for t2,C in opts if t==t2),None)
    raise ValueError(policy)


def analyze(p, policy):
    m=CubicMirror(p,-1,0)
    states=[(A,B) for A in m.points for B in m.points if A!=B]
    nxt={st: ((st[1],C) if (C:=choose(m,*st,policy)) is not None else None) for st in states}
    indeg=Counter(v for v in nxt.values() if v is not None)
    reversed_ok=sum(1 for (A,B),st in nxt.items() if st is not None and nxt.get((st[1],st[0]))==(B,A))
    # state (A,B) -> (B,C), reverse should from (C,B) map -> (B,A).
    seen=set(); lengths=[]; n_terminal=0
    for st in states:
        if st in seen:continue
        path=[]; local={}; cur=st
        while cur is not None and cur not in seen and cur not in local:
            local[cur]=len(path);path.append(cur);cur=nxt.get(cur)
        if cur in local: lengths.append(len(path)-local[cur])
        elif cur is None: n_terminal+=1
        seen.update(path)
    multiple=sum(v>1 for v in indeg.values())
    full=all(nxt[st] is not None for st in states)
    return dict(p=p, policy=policy, points=len(m.points), states=len(states), defined=sum(x is not None for x in nxt.values()), terminal=sum(x is None for x in nxt.values()), indegree0=len(states)-len(indeg), indegree2plus=multiple, max_indegree=max(indeg.values()) if indeg else 0, reverse_branch_success=reversed_ok, cycles=len(lengths), longest_cycle=max(lengths,default=0), longest_cycle_factors=factors(max(lengths,default=0)), permutation=full and multiple==0 and len(indeg)==len(states), sample_noninverse=next(((A,B,st) for (A,B),st in nxt.items() if st is not None and nxt.get((st[1],st[0]))!=(B,A)),None))

def factors(n):
    ans=[]; f=2
    while f*f<=n:
        if n%f==0:
            t=0
            while n%f==0:n//=f;t+=1
            ans.append((f,t))
        f+=1
    if n>1:ans.append((n,1))
    return ans


def counterexample_rescaling(m):
    p=m.p
    for A in m.points:
      for B in m.points:
        if A==B:continue
        opts,h=options_with_t(m,A,B)
        if len(opts)!=2:continue
        u,v=h.reflected_direction
        c2,c1,c0=h.polynomial
        s=min(m.squares[h.discriminant]); t=((-c1+s)*pow(2*c2,-1,p))%p
        original=( (B[0]+t*u)%p,(B[1]+t*v)%p)
        for alpha in range(2, min(p,15)):
          U,V=(alpha*u)%p,(alpha*v)%p
          C2=(-U*U*U)%p; C1=(V*V-3*B[0]*U*U)%p; C0=(2*B[1]*V-(3*B[0]*B[0]+m.a)*U)%p
          disc=(C1*C1-4*C2*C0)%p
          if disc not in m.squares:continue
          s2=min(m.squares[disc]); t2=((-C1+s2)*pow(2*C2,-1,p))%p
          chosen=((B[0]+t2*U)%p,(B[1]+t2*V)%p)
          if chosen!=original:
            return dict(A=A,B=B,original_w=(u,v),alpha=alpha,original_choice=original,rescaled_choice=chosen,options=[C for t,C in opts])
    return None

if __name__=='__main__':
    for p in [101,103]:
        for policy in ['plus_sqrt','positive_t','min_point']:
            print(analyze(p,policy))
        print('rescaling witness',p,counterexample_rescaling(CubicMirror(p,-1,0)))
