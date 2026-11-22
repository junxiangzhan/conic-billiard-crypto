import unittest
from conic_billiard.toy_dh import demo


class ToyECDHTests(unittest.TestCase):
    def test_reference_ecc_construction_and_attack(self):
        result = demo()
        self.assertEqual(result['subgroup_order'],67)
        self.assertEqual(result['BSGS_recovers_Alice_secret'],17)
        self.assertEqual(tuple(result['shared_collision_point']), (647,425))


if __name__ == '__main__':
    unittest.main()
