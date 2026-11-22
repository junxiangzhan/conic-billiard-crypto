"""Physics-compatible real first-hit reflection tests."""
import math
import unittest
from conic_billiard.real_oval import oval, next_hit, forward_chain, distance, cubic_residual
from conic_billiard.cubic_mirror import CubicMirror


class TestRealOval(unittest.TestCase):
    def test_oval_parametrization(self):
        for i in range(1,100):
            P=oval(i*math.pi/100)
            self.assertTrue(-1 <= P[0] <= 0)
            self.assertLess(abs(cubic_residual(P)),1e-14)

    def test_first_intersection_is_oval(self):
        A=(-1.,0.)
        B=(-.5,math.sqrt(3/8))
        hit=next_hit(A,B)
        self.assertLess(abs(hit.point[0]-(-.196619340)),1e-7)
        self.assertGreater(hit.other[0],1)
        self.assertGreater(hit.other_t,hit.t)

    def test_reverse_first_hit(self):
        A,B=oval(.37),oval(1.18)
        C=next_hit(A,B).point
        self.assertLess(distance(next_hit(C,B).point,A),1e-12)

    def test_short_chain(self):
        chain=forward_chain(oval(.37),oval(1.18),12)
        backwards=forward_chain(chain[-1],chain[-2],12)
        self.assertLess(max(distance(P,Q) for P,Q in zip(backwards,chain[::-1])),1e-8)

    def test_finite_field_two_candidates_not_ray(self):
        model=CubicMirror(101,-1,0)
        result=model.successors((100,0),(0,0))
        self.assertEqual(set(result.next_points),{(100,0),(1,0)})
        real_hit=next_hit((-1.,0.),(0.,0.))
        self.assertLess(distance(real_hit.point,(-1.,0.)),1e-12)

if __name__=='__main__': unittest.main()
