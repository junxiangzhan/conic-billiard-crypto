import unittest
from conic_billiard.cubic_mirror import CubicMirror
from conic_billiard.positive_branch_experiment import analyze, choose

class PositiveRootTests(unittest.TestCase):
    def test_rule_is_deterministic_but_not_permutation(self):
        result=analyze(101, 'plus_sqrt')
        self.assertEqual(result['states'],10506)
        self.assertEqual(result['defined'],5406)
        self.assertEqual(result['indegree2plus'],1339)
        self.assertEqual(result['longest_cycle'],10)
        self.assertFalse(result['permutation'])

    def test_positive_half_rejects_many(self):
        result=analyze(101, 'positive_t')
        self.assertEqual(result['defined'],4080)
        self.assertEqual(result['longest_cycle'],4)

    def test_plus_selected_point(self):
        m=CubicMirror(101,-1,0)
        self.assertEqual(choose(m,(100,0),(0,0),'plus_sqrt'),(100,0))

if __name__=='__main__':
    unittest.main()
