"""Run real-oval billiard tests and compare with the finite-field relation."""
from __future__ import annotations
import json
import math
import statistics
import random
from collections import Counter
from pathlib import Path
from .real_oval import oval, next_hit, forward_chain, distance
from .cubic_mirror import CubicMirror


def run(sample_count: int=10_000, chain_count: int=200, chain_length: int=80, seed: int=20261118):
    rng=random.Random(seed)
    stats=Counter()
    max_reverse_err=0.0
    max_cubic_residual=0.0
    max_chain_reverse_err=0.0
    reverse_error_by_length={10:[],20:[],40:[],80:[]}
    initial_states=[]
    examples=[]
    for _ in range(sample_count):
        # Avoid near-tangent impacts to keep double-precision experiments stable.
        while True:
            ta,tb=(rng.random()*math.pi,rng.random()*math.pi)
            A,B=oval(ta),oval(tb)
            if distance(A,B)>.055:
                break
        hit=next_hit(A,B)
        stats['success']+=1
        if hit.other is not None:
            stats['finite_third_intersection']+=1
            if hit.other_t is not None and hit.other_t>0:
                stats['third_point_in_front_of_B']+=1
            else:
                stats['third_point_behind_B']+=1
        else:
            stats['third_point_at_infinity_or_tangency']+=1
        rev=next_hit(hit.point,B)
        err=distance(rev.point,A)
        max_reverse_err=max(err,max_reverse_err)
        from .real_oval import cubic_residual
        max_cubic_residual=max(max_cubic_residual,abs(cubic_residual(hit.point)))
        if len(examples)<3 and hit.other is not None and hit.other_t>hit.t:
            examples.append({'A':A,'B':B,'physical_C':hit.point,
                             'extra_D':hit.other,'t_C':hit.t,'t_D':hit.other_t})
        if len(initial_states)<chain_count:
            initial_states.append((A,B))
    for A,B in initial_states:
        chain=forward_chain(A,B,chain_length)
        for C in chain:
            from .real_oval import cubic_residual
            max_cubic_residual=max(max_cubic_residual,abs(cubic_residual(C)))
        # Reverse at different depths; the map is exactly reversible but
        # chaotic-like numerical sensitivity can amplify rounding error.
        for depth in reverse_error_by_length:
            prefix=chain[:depth+2]
            backwards=forward_chain(prefix[-1],prefix[-2],depth)
            expected=prefix[::-1]
            err=max(distance(P,Q) for P,Q in zip(backwards,expected))
            reverse_error_by_length[depth].append(err)
            if depth==chain_length:
                max_chain_reverse_err=max(max_chain_reverse_err,err)
    mirror=CubicMirror(101,-1,0)
    counts=mirror.state_graph_counts()
    arithmetic_example=mirror.successors((-1%101,0),(0,0))
    # Deterministic boundary angles used for visual illustration.
    visual_chain=forward_chain(oval(.37),oval(1.18),14)
    return {
      'model':'real Euclidean billiard inside bounded oval of y^2=x^3-x',
      'random_initial_pairs':sample_count,
      'root_statistics':dict(stats),
      'maximum_one_step_reverse_error':max_reverse_err,
      'chain_count':chain_count,'additional_bounces_per_chain':chain_length,
      'maximum_full_chain_reverse_error_in_double_precision':max_chain_reverse_err,
      'reverse_roundoff_by_depth':{str(depth):{'median':statistics.median(errors),
                       'max':max(errors),'fraction_above_1e_minus_4':sum(err>1e-4 for err in errors)/len(errors)}
                       for depth,errors in reverse_error_by_length.items()},
      'roundoff_warning':'large reverse errors in 80-step double-precision runs reflect numerical instability, not mathematical non-invertibility; compare optional 100-digit verification',
      'maximum_cubic_residual':max_cubic_residual,
      'two_positive_intersection_examples':examples,
      'finite_field_101':{'affine_points':len(mirror.points),
                          'ordered_states':counts['states_tested'],
                          'status_counts':counts['status_counts'],
                          'number_of_reflected_intersections_for_A_minus1_B_zero':len(arithmetic_example.next_points),
                          'reflected_intersections_A_minus1_B_zero':arithmetic_example.next_points,
                          'warning':'finite field lacks positive-ray parameter and nearest-first ordering'},
      'plot_trajectory':visual_chain
    }


def plot(output: dict, path: Path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(1,2,figsize=(11,4.3))
    x=[-1+i/1199 for i in range(1200)]
    y=[math.sqrt(max(0,z**3-z)) for z in x]
    for ax in axs:
        ax.plot(x,y,linewidth=1.5)
        ax.plot(x,[-z for z in y],linewidth=1.5)
        ax.set_aspect('equal',adjustable='box')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.grid(alpha=.2)
    ex=output['two_positive_intersection_examples'][0]
    A,B,C,D=(ex[k] for k in ['A','B','physical_C','extra_D'])
    axs[0].plot([A[0],B[0],C[0]],[A[1],B[1],C[1]],marker='o',linewidth=1.8)
    for name,P in [('A',A),('B',B),('C',C)]:
        axs[0].annotate(name,P,xytext=(4,4),textcoords='offset points')
    axs[0].set_title('Real first-hit reflection on bounded oval')
    chain=output['plot_trajectory']
    axs[1].plot([P[0] for P in chain],[P[1] for P in chain],marker='.',linewidth=1)
    axs[1].set_title('First 14 reflections (true Euclidean rays)')
    fig.tight_layout()
    fig.savefig(path,dpi=170)
    plt.close(fig)


def main():
    out=run()
    base=Path(__file__).resolve().parent.parent / 'results'
    base.mkdir(parents=True,exist_ok=True)
    f=base/'real_oval_first_hit.json'
    f.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding='utf-8')
    try:
        plot(out,base/'real_oval_first_hit.png')
    except ImportError:
        print("matplotlib not installed: skipped optional plot")
    print(json.dumps({k:v for k,v in out.items() if k not in ('plot_trajectory','two_positive_intersection_examples')}, indent=2))
    print('saved:',f)

if __name__=='__main__': main()
