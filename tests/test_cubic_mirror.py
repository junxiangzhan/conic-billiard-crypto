import unittest
from conic_billiard.cubic_mirror import CubicMirror, polynomial_invariant_scan


class CubicMirrorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c101=CubicMirror(101)

    def test_affine_points(self):
        self.assertEqual(len(self.c101.points),103)
        self.assertTrue(self.c101.valid((2,39)))

    def test_no_rational_intersection(self):
        r=self.c101.successors((0,0),(2,39))
        self.assertEqual(r.reflected_direction,(98,1))
        self.assertEqual(r.polynomial,(27,48,10))
        self.assertEqual(r.discriminant,12)
        self.assertEqual(r.next_points,())

    def test_two_intersections(self):
        r=self.c101.successors((1,0),(2,39))
        self.assertEqual(r.reflected_direction,(52,63))
        self.assertEqual(r.next_points,((50,35),(96,48)))

    def test_entire_101_relation_reversible(self):
        r=self.c101.state_graph_counts()
        self.assertEqual(r['states_tested'],10506)
        self.assertEqual(r['status_counts']['two-next'],5164)
        self.assertEqual(r['status_counts']['no-affine-next'],5050)
        self.assertEqual(r['edges_tested'],10620)
        self.assertEqual(r['reverse_edges_verified'],r['edges_tested'])
        self.assertEqual(r['chord_norm_squareclass_verified'],r['edges_tested'])

    def test_polynomial_invariants(self):
        for deg, nullity in ((1,1),(2,1),(3,3),(4,11)):
            r=polynomial_invariant_scan(self.c101)['degree_results'][deg]
            self.assertEqual(r['nullity_upper_bound'],nullity)

    def test_isotropic_normal_in_1009(self):
        m=CubicMirror(1009)
        r=m.successors((2,835),(605,465))
        self.assertEqual(r.status,'isotropic-normal')

    def test_min_and_max_policies_differ(self):
        self.assertEqual(self.c101.deterministic_cycle_scan('min')['max_cycle'],4)
        self.assertEqual(self.c101.deterministic_cycle_scan('max')['max_cycle'],14)

if __name__=='__main__': unittest.main()
