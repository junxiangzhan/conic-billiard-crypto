import unittest
from conic_billiard.compression import JacobiCodec, InvalidEncoding, sqrt_mod
from conic_billiard.jacobi import JacobiGroup, ellipse_lift
from conic_billiard.core import Conic


class TestCompression(unittest.TestCase):
    def test_complete_exhaustive_roundtrip_small_fields(self):
        # Every affine rational Jacobi point for elliptic and hyperbolic models.
        for p in (101,103):
            for eps in (1,-1):
                for m in (2,3,17):
                    if m in (0,eps%p):continue
                    codec=JacobiCodec(p,m,eps)
                    count=0
                    for s in range(p):
                        for c in range(p):
                            if (c*c+eps*s*s-1)%p:continue
                            for d in range(p):
                                U=(s,c,d)
                                if codec.valid_point(U):
                                    self.assertEqual(codec.decode(codec.encode(U)),U)
                                    self.assertEqual(codec.encode(codec.decode(codec.encode(U))), codec.encode(U))
                                    count+=1
                    self.assertGreater(count,30)

    def test_exceptional_c_negative_one_both_signs(self):
        p=103
        for eps in (1,-1):
            codec=JacobiCodec(p,2,eps)
            for d in (1,p-1):
                U=(0,p-1,d)
                data=codec.encode(U)
                self.assertEqual(codec.decode(data),U)
                self.assertEqual(int.from_bytes(data,'little')&((1<<p.bit_length())-1),p)

    def test_nonsquare_invalid_and_pole_and_padding(self):
        codec=JacobiCodec(103,17)
        bits=codec.field_bits
        # Payload strictly bigger than p must be rejected.
        self.assertRaises(InvalidEncoding,codec.decode,(104).to_bytes(codec.byte_length,'little'))
        # Zero square root with nonzero sign is noncanonical if one exists.
        found=False
        for t in range(103):
            if ((1+t*t)**2-4*17*t*t)%103==0:
                self.assertRaises(InvalidEncoding,codec.decode,(t|(1<<bits)).to_bytes(codec.byte_length,'little'))
                found=True;break
        # Non-residues are rejected.
        for t in range(103):
            rhs=((1+t*t)**2-4*17*t*t)%103
            if sqrt_mod(rhs,103) is None:
                self.assertRaises(InvalidEncoding,codec.decode,t.to_bytes(codec.byte_length,'little'))
                break
        self.assertRaises(InvalidEncoding,codec.decode,b'')
        # For p=1009, k=10 and payload+sign occupies only 11 of 16 bits.
        padded=JacobiCodec(1009,635)
        self.assertRaises(InvalidEncoding,padded.decode,(1<<12).to_bytes(2,'little'))
        self.assertRaises(InvalidEncoding,codec.decode,(1<<8).to_bytes(codec.byte_length+1,'little'))

    def test_quartic_zero_root_sign_and_hyperbola_poles(self):
        p=101
        t=2
        m=((1+t*t)**2 * pow(4*t*t,-1,p))%p
        codec=JacobiCodec(p,m,1)
        self.assertEqual(sqrt_mod(((1+t*t)**2-4*m*t*t)%p,p),0)
        bits=codec.field_bits
        self.assertEqual(codec.decode(t.to_bytes(codec.byte_length,'little'))[2],0)
        self.assertRaises(InvalidEncoding, codec.decode,
                          (t|(1<<bits)).to_bytes(codec.byte_length,'little'))
        hyper=JacobiCodec(101,3,-1)
        for pole in (1,100):
            self.assertRaises(InvalidEncoding,hyper.decode,
                              pole.to_bytes(hyper.byte_length,'little'))

    def test_large_curve25519_32_byte_roundtrip(self):
        p=2**255-19
        m=20800338683988658368647408995589388737092878452977063003340006470870624536394
        a,b=5,3
        t=4604094445982880512620096533988940857801848594535175639138443919056596161201
        den=pow((1+t*t)%p,-1,p)
        A=(a,0)
        B=(a*(1-t*t)*den%p,2*b*t*den%p)
        lam,J,G=ellipse_lift(Conic.ellipse(p,a,b),a,b,A,B)
        codec=JacobiCodec(p,m)
        self.assertEqual(J.m, m)
        self.assertEqual(codec.byte_length,32)
        for k in (0,1,2,3,8,42,255):
            U=J.mul(k,G)
            self.assertEqual(codec.decode(codec.encode(U)),U)
        r=2**252+27742317777372353535851937790883648493
        H=J.mul(8,G)
        encoded=codec.encode(H)
        self.assertEqual(len(encoded),32)
        self.assertEqual(codec.decode(encoded),H)
        self.assertEqual(J.mul(r,codec.decode(encoded)),J.identity())

    def test_little_endian_sign_and_chart(self):
        codec=JacobiCodec(1009,635)
        # Known Jacobi point corresponding to collision (730,932).
        U=(647,146,528)
        c=codec.encode(U)
        self.assertEqual(codec.decode(c),U)
        t=int.from_bytes(c,'little')&((1<<codec.field_bits)-1)
        self.assertEqual(t,176)
        self.assertEqual(codec.decode_hex(c.hex()),U)
        self.assertRaises(ValueError,codec.encode,(647,146,0))


if __name__ == '__main__':
    unittest.main()
