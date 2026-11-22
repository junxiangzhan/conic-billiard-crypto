"""Exhaustive state graph exploration for tiny prime fields only."""
from __future__ import annotations
from collections import Counter
from dataclasses import asdict, dataclass
from .core import Conic, State, UndefinedReflection


@dataclass
class CycleStats:
    curve: str
    p: int
    affine_points: int
    candidate_states: int
    defined_states: int
    undefined_by_reason: dict[str, int]
    periodic_states: int
    cycles: int
    max_period: int
    max_cycle_initial_state: State | None
    largest_prime_factor_of_max_period: int | None
    largest_prime_factor_any_period: int | None
    period_with_largest_prime_factor: int | None
    top_periods: list[int]

    def to_dict(self):
        return asdict(self)


def largest_prime_factor(n: int) -> int | None:
    if n < 2:
        return None
    last = None
    d = 2
    while d*d <= n:
        while n % d == 0:
            last = d
            n //= d
        d = 3 if d == 2 else d+2
    if n > 1:
        last = n
    return last


def analyze(conic: Conic) -> CycleStats:
    points = conic.boundary_points()
    transition: dict[State, State] = {}
    bad: Counter[str] = Counter()
    for A in points:
        for B in points:
            if A == B:
                continue
            state = (A, B)
            try:
                transition[state] = conic.step(state)
            except UndefinedReflection as exc:
                bad[exc.reason] += 1
    # This follows from time reversibility wherever the map is defined.
    assert len(set(transition.values())) == len(transition), 'unexpected non-injectivity'

    visited: set[State] = set()
    lengths: list[int] = []
    cyclic_states = 0
    max_state: State | None = None
    max_period = 0
    for start in transition:
        if start in visited:
            continue
        path: list[State] = []
        offset: dict[State, int] = {}
        curr = start
        while curr in transition and curr not in visited and curr not in offset:
            offset[curr] = len(path)
            path.append(curr)
            curr = transition[curr]
        if curr in offset:
            length = len(path) - offset[curr]
            lengths.append(length)
            cyclic_states += length
            if length > max_period:
                max_period = length
                max_state = curr
        visited.update(path)
    period_best_factor = max(lengths, key=lambda r: (largest_prime_factor(r) or 0), default=None)
    return CycleStats(
        curve=conic.label, p=conic.p, affine_points=len(points),
        candidate_states=len(points)*(len(points)-1),
        defined_states=len(transition), undefined_by_reason=dict(sorted(bad.items())),
        periodic_states=cyclic_states, cycles=len(lengths), max_period=max_period,
        max_cycle_initial_state=max_state,
        largest_prime_factor_of_max_period=largest_prime_factor(max_period),
        largest_prime_factor_any_period=largest_prime_factor(period_best_factor) if period_best_factor else None,
        period_with_largest_prime_factor=period_best_factor,
        top_periods=sorted(set(lengths), reverse=True)[:10],
    )
