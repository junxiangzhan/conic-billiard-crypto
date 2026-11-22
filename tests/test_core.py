import unittest
from conic_billiard import Conic, UndefinedReflection
from conic_billiard.core import matvec4, determinant4


class ReflectionTests(unittest.TestCase):
    def test_stored_small_field_counts(self):
        expected = {
            103: {'ellipse': (104, 60, 0), 'hyperbola': (102, 60, 196), 'parabola': (103, 104, 102)},
            127: {'ellipse': (128, 126, 0), 'hyperbola': (126, 74, 244), 'parabola': (127, 128, 126)},
        }
        from conic_billiard import analyze
        for p, cases in expected.items():
            for name, (points, period, undefined) in cases.items():
                with self.subTest(p=p, curve=name):
                    c = getattr(Conic, name)(p)
                    result = analyze(c)
                    self.assertEqual(result.affine_points, points)
                    self.assertEqual(result.max_period, period)
                    self.assertEqual(sum(result.undefined_by_reason.values()), undefined)
                    self.assertIsNotNone(result.largest_prime_factor_of_max_period)

    def test_matrix_matches_step_and_is_invertible(self):
        for c in (Conic.ellipse(103), Conic.hyperbola(103), Conic.parabola(103)):
            pts = c.boundary_points()
            count = 0
            for A in pts[::7]:
                for B in pts[::11]:
                    if A == B:
                        continue
                    try:
                        direct = c.step((A, B))
                        mat = c.transfer_matrix(A, B)
                    except UndefinedReflection:
                        continue
                    self.assertEqual(direct, matvec4(mat, (A, B), c.p))
                    self.assertNotEqual(determinant4(mat, c.p), 0)
                    # Time reversibility on the same collision point B.
                    self.assertEqual(c.next_point(direct[1], B), A)
                    count += 1
            self.assertGreater(count, 10)

    def test_isotropic_cases_and_coincident_points(self):
        c = Conic.hyperbola(103)
        A = c.boundary_points()[0]
        with self.assertRaises(UndefinedReflection) as ctx:
            c.next_point(A, A)
        self.assertEqual(ctx.exception.reason, 'coincident_input_points')
        # p == 1 (mod 4) allows nonzero isotropic normal directions.
        c2 = Conic.ellipse(101)
        points = c2.boundary_points()
        witness = next(B for B in points if sum(x*x for x in c2.normal(B)) % 101 == 0)
        A2 = next(A for A in points if A != witness)
        with self.assertRaises(UndefinedReflection) as ctx2:
            c2.next_point(A2, witness)
        self.assertEqual(ctx2.exception.reason, 'isotropic_normal')

    def test_large_field_single_step_not_exhaustive_search(self):
        p = 2**256-189
        curve = Conic.ellipse(p,5,3)
        t = 42
        z = pow(1+t*t, -1, p)
        A = (5,0)
        B = (5*(1-t*t)*z % p, 3*2*t*z % p)
        C = curve.next_point(A,B)
        self.assertTrue(curve.contains(C))
        self.assertEqual(matvec4(curve.transfer_matrix(A,B), (A,B),p), (B,C))

    def test_new_state_depends_on_previous_state(self):
        c = Conic.ellipse(103)
        pts = c.boundary_points()
        B = pts[19]
        mats = []
        for A in pts:
            if A == B:
                continue
            try:
                mats.append(c.transfer_matrix(A, B))
            except UndefinedReflection:
                pass
            if len(mats) > 3:
                break
        # Dependence on incoming chord is not mathematically forced for all choices,
        # but this particular sample does vary.
        self.assertGreater(len(set(mats)), 1)


if __name__ == '__main__':
    unittest.main()
