import unittest
from conic_billiard.demos import small_field_dh, compressed_dh, verify_public, P, M
from conic_billiard.compression import JacobiCodec
from conic_billiard.jacobi import JacobiGroup


class SmallFieldDHTests(unittest.TestCase):
    def test_reference_construction_and_attack(self):
        result = small_field_dh()
        self.assertEqual(result['subgroup_order'], 67)
        self.assertEqual(result['BSGS_recovers_Alice_secret'], 17)
        self.assertEqual(tuple(result['shared_collision_point']), (647, 425))


class CompressedDHTests(unittest.TestCase):
    def test_exchange_and_32_byte_wire(self):
        result = compressed_dh()
        self.assertTrue(result['shared_point_agrees'])
        self.assertEqual(result['encoded_public_key_bytes'], 32)
        self.assertEqual(len(result['alice_public_hex']), 64)

    def test_reject_identity_and_torsion(self):
        J = JacobiGroup(P, M)
        codec = JacobiCodec(P, M)
        identity = (0, 1, 1)
        self.assertEqual(codec.decode(codec.encode(identity)), identity)
        with self.assertRaises(ValueError):
            verify_public(J, identity)
        # A valid Jacobi state outside the prime-order subgroup.
        torsion = (0, P - 1, 1)
        self.assertTrue(codec.valid_point(torsion))
        self.assertEqual(codec.decode(codec.encode(torsion)), torsion)
        with self.assertRaises(ValueError):
            verify_public(J, torsion)


if __name__ == '__main__':
    unittest.main()
