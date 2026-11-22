import unittest
from conic_billiard import Conic, UndefinedReflection


class HyperbolaTests(unittest.TestCase):
    def test_hyperbola_equation_and_points(self):
        c = Conic.hyperbola(103, 5, 3)
        self.assertTrue(c.contains((5, 0)))
        self.assertTrue(c.contains((-5 % 103, 0)))
        self.assertFalse(c.contains((0, 0)))
        self.assertEqual(len(c.boundary_points()), 102)

    def test_a_specific_periodic_orbit(self):
        c = Conic.hyperbola(103, 5, 3)
        A, B = (5, 0), (98, 0)
        self.assertEqual(c.step((A, B)), (B, A))
        self.assertEqual(c.step((B, A)), (A, B))

    def test_nonextension_at_affine_infinity(self):
        c = Conic.hyperbola(103, 5, 3)
        pts = c.boundary_points()
        witness = None
        for A in pts:
            for B in pts:
                if A == B:
                    continue
                try:
                    c.next_point(A, B)
                except UndefinedReflection as exc:
                    if exc.reason == 'next_intersection_at_infinity':
                        witness = (A, B)
                        break
            if witness:
                break
        self.assertIsNotNone(witness)
        # We deliberately do not define a made-up continuation.


if __name__ == '__main__':
    unittest.main()
