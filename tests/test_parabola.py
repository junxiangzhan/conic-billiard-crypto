import unittest
from conic_billiard import Conic, UndefinedReflection


class ParabolaTests(unittest.TestCase):
    def test_known_fixed_matrix_linearization(self):
        p = 103
        c = Conic.parabola(p)
        u, v = 0, 2
        k = (u-v)**2 * pow(1+4*u*v, -1, p) % p
        self.assertEqual(k, 4)
        h = (2+4*k) % p
        first = (u, v)
        state = first
        for _ in range(104):
            A = (state[0], state[0]**2 % p)
            B = (state[1], state[1]**2 % p)
            C = c.next_point(A, B)
            step = (state[1], (h*state[1]-state[0]) % p)
            self.assertEqual(C, (step[1], step[1]*step[1] % p))
            state = step
        self.assertEqual(state, first)

    def test_parabola_not_assumed_secure(self):
        # This is intentionally a regression test for a reducible 2x2 model.
        p = 127
        c = Conic.parabola(p)
        A, B = (0, 0), (2, 4)
        state = (A, B)
        for _ in range(128):
            state = c.step(state)
        self.assertEqual(state, (A, B))


if __name__ == '__main__':
    unittest.main()
