"""Reproducible cubic *mirror* reflection experiment (not EC group law).

Run: python -m conic_billiard.cubic_experiment [--out results/cubic_mirror.json]
"""
from __future__ import annotations
import argparse, json, time
from .cubic_mirror import CubicMirror, polynomial_invariant_scan


def run_experiment():
    m101=CubicMirror(101,-1,0)
    m1009=CubicMirror(1009,-1,0)
    r101=m101.state_graph_counts()
    r1009=m1009.state_graph_counts()  # COMPLETE enumeration, ~1.08 million states
    pol101=polynomial_invariant_scan(m101)
    pol1009=polynomial_invariant_scan(m1009)
    A,B=(0,0),(2,39)
    A2,B2=(1,0),(2,39)
    r0=m101.successors(A,B)
    r1=m101.successors(A2,B2)
    reach=m1009.breadth_reachable((499,269),(91,322),12)
    return {
        'curve':'y^2=x^3-x mod p',
        'model':'specular reflection with bilinear form dx^2+dy^2',
        '101':r101, '1009':r1009,
        '101_deterministic_min_policy':m101.deterministic_cycle_scan('min'),
        '101_deterministic_max_policy':m101.deterministic_cycle_scan('max'),
        'polynomial_invariant_scan_101':pol101,
        'polynomial_invariant_scan_1009':pol1009,
        '101_examples':{
            'no_rational_next':{'A':A,'B':B,'polynomial':r0.polynomial,
                                'discriminant':r0.discriminant,'result':r0.next_points},
            'two_rational_next':{'A':A2,'B':B2,'polynomial':r1.polynomial,
                                 'discriminant':r1.discriminant,'result':r1.next_points},
        },
        '1009_branch_reach_counts_for_seed_499_269__91_322':reach,
        'caveats':[
            'Finite fields have no physical ordering of ray-intersection parameters.',
            't=0 contact at B excluded as a new hit; finite roots only.',
            'Degree at most 4 polynomial invariant scan excludes only this class; not rational or algebraic invariants.',
            'Min/max successor choice is an arbitrary rule; not a canonical or reversible billiard action.',
            'No direct O(log n) fast jump or secure DH for the cubic mirror model is demonstrated.',
        ],
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',default=None,help='Write JSON report to this path')
    args=parser.parse_args()
    t0=time.perf_counter()
    r=run_experiment()
    print('Full scans:',[(k,r[k]['states_tested'],r[k]['status_counts']) for k in ('101','1009')])
    print('Reversibility:',[(k,r[k]['reverse_edges_verified'],r[k]['edges_tested']) for k in ('101','1009')])
    print('Policy cycle max:',r['101_deterministic_min_policy']['max_cycle'],r['101_deterministic_max_policy']['max_cycle'])
    print('Polynomial invariant nullities, degree <=4:',[(k,r[k]['degree_results'][4]['nullity_upper_bound']) for k in ('polynomial_invariant_scan_101','polynomial_invariant_scan_1009')])
    print('Elapsed seconds:',round(time.perf_counter()-t0,3))
    if args.out:
        from pathlib import Path
        p=Path(args.out)
        p.parent.mkdir(parents=True,exist_ok=True)
        p.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        print('Written:',p)

if __name__=='__main__':
    main()
