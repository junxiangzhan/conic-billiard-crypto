"""CLI for a small, non-cryptographic finite-field conic billiard laboratory."""
from __future__ import annotations
import argparse
import json
from .core import Conic, UndefinedReflection
from .search import analyze


def choose(kind: str, p: int, a: int, b: int) -> Conic:
    if kind == 'ellipse':
        return Conic.ellipse(p, a, b)
    if kind == 'hyperbola':
        return Conic.hyperbola(p, a, b)
    if kind == 'parabola':
        return Conic.parabola(p)
    if kind == 'circle':
        return Conic.circle(p, a)
    raise ValueError(kind)


def main():
    parser = argparse.ArgumentParser(description='Toy F_p conic reflection lab; NOT secure cryptography')
    parser.add_argument('command', choices=['scan', 'trace'])
    parser.add_argument('--p', type=int, default=103)
    parser.add_argument('--kind', choices=['ellipse', 'hyperbola', 'parabola', 'circle', 'all'], default='all')
    parser.add_argument('--a', type=int, default=5)
    parser.add_argument('--b', type=int, default=3)
    parser.add_argument('--A', type=int, nargs=2, metavar=('X', 'Y'))
    parser.add_argument('--B', type=int, nargs=2, metavar=('X', 'Y'))
    parser.add_argument('--steps', type=int, default=12)
    args = parser.parse_args()
    if args.command == 'scan':
        kinds = ('ellipse', 'hyperbola', 'parabola') if args.kind == 'all' else (args.kind,)
        for kind in kinds:
            print(json.dumps(analyze(choose(kind, args.p, args.a, args.b)).to_dict(), ensure_ascii=False))
        return
    if args.kind == 'all' or args.A is None or args.B is None:
        parser.error('trace requires --kind, --A x y, --B x y')
    curve = choose(args.kind, args.p, args.a, args.b)
    state = ((args.A[0] % args.p, args.A[1] % args.p),
             (args.B[0] % args.p, args.B[1] % args.p))
    print(f'curve={curve.label} F_{curve.p}: initial state {state}')
    for i in range(args.steps):
        try:
            state = curve.step(state)
        except UndefinedReflection as exc:
            print(f'undefined at transition {i}: {exc.reason}')
            break
        print(f'transition {i+1}: {state}')


if __name__ == '__main__':
    main()
