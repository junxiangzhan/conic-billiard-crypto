import unittest
from conic_billiard.compressed_dh_demo import demo,verify_public,P,M,R
from conic_billiard.compression import JacobiCodec
from conic_billiard.jacobi import JacobiGroup


class TestCompressedDH(unittest.TestCase):
    def test_exchange_and_32_byte_wire(self):
        result=demo()
        self.assertTrue(result['shared_point_agrees'])
        self.assertEqual(result['encoded_public_key_bytes'],32)
        self.assertEqual(len(result['alice_public_hex']),64)

    def test_reject_identity_and_torsion(self):
        J=JacobiGroup(P,M)
        codec=JacobiCodec(P,M)
        identity=(0,1,1)
        self.assertEqual(codec.decode(codec.encode(identity)),identity)
        with self.assertRaises(ValueError):verify_public(J,identity)
        # A valid Jacobi state outside the prime-order subgroup.
        torsion=(0,P-1,1)
        self.assertTrue(codec.valid_point(torsion))
        assert codec.decode(codec.encode(torsion))==torsion
        with self.assertRaises(ValueError): verify_public(J,torsion)


if __name__=='__main__':unittest.main()
