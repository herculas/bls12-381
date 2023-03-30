#include <gtest/gtest.h>

#include "group/g2_affine.h"
#include "group/g2_projective.h"
#include "group/g2_prepared.h"
#include "group/constant.h"
#include "scalar/scalar.h"

using bls12_381::field::Fp;
using bls12_381::field::Fp2;
using bls12_381::group::G2Affine;
using bls12_381::group::G2Projective;
using bls12_381::group::G2Prepared;
using bls12_381::scalar::Scalar;

TEST(TestG2, OnCurve) {
    EXPECT_TRUE(G2Affine::identity().is_on_curve());
    EXPECT_TRUE(G2Affine::generator().is_on_curve());
    EXPECT_TRUE(G2Projective::identity().is_on_curve());
    EXPECT_TRUE(G2Projective::generator().is_on_curve());

    Fp2 const z{
            Fp(
                    {
                            0xba7afa1f9a6fe250, 0xfa0f5b595eafe731, 0x3bdc477694c306e7,
                            0x2149be4b3949fa24, 0x64aa6e0649b2078c, 0x12b108ac33643c3e,
                    }
            ),
            Fp(
                    {
                            0x125325df3d35b5a8, 0xdc469ef5555d7fe3, 0x02d716d2443106a9,
                            0x05a1db59a6ff37d0, 0x7cf7784e5300bb8f, 0x16a88922c7a5e844,
                    }
            ),
    };
    G2Affine const gen = G2Affine::generator();
    G2Projective const test1{gen.get_x() * z, gen.get_y() * z, z};
    G2Projective const test2{z, gen.get_y() * z, z};

    EXPECT_TRUE(test1.is_on_curve());
    EXPECT_FALSE(test2.is_on_curve());
}

TEST(TestG2, AffineEquality) {
    G2Affine const a = G2Affine::generator();
    G2Affine const b = G2Affine::identity();

    EXPECT_TRUE(a == a);
    EXPECT_TRUE(b == b);
    EXPECT_TRUE(a != b);
    EXPECT_TRUE(b != a);
}

TEST(TestG2, ProjectiveEquality) {
    G2Projective const a = G2Projective::generator();
    G2Projective const b = G2Projective::identity();

    EXPECT_TRUE(a == a);
    EXPECT_TRUE(b == b);
    EXPECT_TRUE(a != b);
    EXPECT_TRUE(b != a);

    Fp2 const z{
            Fp(
                    {
                            0xba7afa1f9a6fe250, 0xfa0f5b595eafe731, 0x3bdc477694c306e7,
                            0x2149be4b3949fa24, 0x64aa6e0649b2078c, 0x12b108ac33643c3e,
                    }
            ),
            Fp(
                    {
                            0x125325df3d35b5a8, 0xdc469ef5555d7fe3, 0x02d716d2443106a9,
                            0x05a1db59a6ff37d0, 0x7cf7784e5300bb8f, 0x16a88922c7a5e844,
                    }
            ),
    };

    G2Projective const p1{a.get_x() * z, a.get_y() * z, z};
    G2Projective const p2{a.get_x() * z, -a.get_y() * z, z};
    G2Projective const p3{z, a.get_y() * z, z};

    EXPECT_TRUE(p1.is_on_curve());
    EXPECT_TRUE(p2.is_on_curve());
    EXPECT_FALSE(p3.is_on_curve());

    EXPECT_TRUE(a == p1);
    EXPECT_TRUE(p1 == a);
    EXPECT_TRUE(b != p1);
    EXPECT_TRUE(p1 != b);

    EXPECT_TRUE(a != p2);
    EXPECT_TRUE(p2 != a);
    EXPECT_TRUE(b != p2);
    EXPECT_TRUE(p2 != b);

    EXPECT_TRUE(a != p3);
    EXPECT_TRUE(p3 != a);
    EXPECT_TRUE(b != p3);
    EXPECT_TRUE(p3 != b);
}

TEST(TestG2, ProjectiveToAffine) {
    G2Projective const a = G2Projective::generator();
    G2Projective const b = G2Projective::identity();

    EXPECT_TRUE(G2Affine(a).is_on_curve());
    EXPECT_FALSE(G2Affine(a).is_identity());
    EXPECT_TRUE(G2Affine(b).is_on_curve());
    EXPECT_TRUE(G2Affine(b).is_identity());

    Fp2 const z{
            Fp(
                    {
                            0xba7afa1f9a6fe250, 0xfa0f5b595eafe731, 0x3bdc477694c306e7,
                            0x2149be4b3949fa24, 0x64aa6e0649b2078c, 0x12b108ac33643c3e,
                    }
            ),
            Fp(
                    {
                            0x125325df3d35b5a8, 0xdc469ef5555d7fe3, 0x02d716d2443106a9,
                            0x05a1db59a6ff37d0, 0x7cf7784e5300bb8f, 0x16a88922c7a5e844,
                    }
            ),
    };

    G2Projective const c{a.get_x() * z, a.get_y() * z, z};
    EXPECT_EQ(G2Affine(c), G2Affine::generator());
}

TEST(TestG2, Doubleing) {
    G2Projective const temp1 = G2Projective::identity().doubles();
    G2Projective const temp2 = G2Projective::generator().doubles();

    EXPECT_TRUE(temp1.is_identity());
    EXPECT_TRUE(temp1.is_on_curve());
    EXPECT_FALSE(temp2.is_identity());
    EXPECT_TRUE(temp2.is_on_curve());

    G2Affine const temp3{
            Fp2{
                    Fp(
                            {
                                    0xe9d9e2da9620f98b, 0x54f1199346b97f36, 0x3db3b820376bed27,
                                    0xcfdb31c9b0b64f4c, 0x41d7c12786354493, 0x05710794c255c064,
                            }
                    ),
                    Fp(
                            {
                                    0xd6c1d3ca6ea0d06e, 0xda0cbd905595489f, 0x4f5352d43479221d,
                                    0x8ade5d736f8c97e0, 0x48cc8433925ef70e, 0x08d7ea71ea91ef81,
                            }
                    ),
            },
            Fp2{
                    Fp(
                            {
                                    0x15ba26eb4b0d186f, 0x0d086d64b7e9e01e, 0xc8b848dd652f4c78,
                                    0xeecf46a6123bae4f, 0x255e8dd8b6dc812a, 0x164142af21dcf93f,
                            }
                    ),
                    Fp(
                            {
                                    0xf9b4a1a895984db4, 0xd417b114cccff748, 0x6856301fc89f086e,
                                    0x41c777878931e3da, 0x3556b155066a2105, 0x00acf7d325cb89cf,
                            }
                    ),
            },
            false
    };

    EXPECT_EQ(G2Affine(temp2), temp3);
}

TEST(TestG2, Add1) {
    G2Projective const a = G2Projective::identity();
    G2Projective const b = G2Projective::identity();
    G2Projective const c = a + b;

    EXPECT_TRUE(c.is_identity());
    EXPECT_TRUE(c.is_on_curve());
}

TEST(TestG2, Add2) {
    G2Projective const a = G2Projective::identity();
    G2Projective b = G2Projective::generator();

    Fp2 const z{
            Fp({
                       0xba7afa1f9a6fe250, 0xfa0f5b595eafe731, 0x3bdc477694c306e7,
                       0x2149be4b3949fa24, 0x64aa6e0649b2078c, 0x12b108ac33643c3e,
               }),
            Fp({
                       0x125325df3d35b5a8, 0xdc469ef5555d7fe3, 0x02d716d2443106a9,
                       0x05a1db59a6ff37d0, 0x7cf7784e5300bb8f, 0x16a88922c7a5e844,
               }),

    };
    b = G2Projective{b.get_x() * z, b.get_y() * z, z};
    G2Projective const c = a + b;

    EXPECT_FALSE(c.is_identity());
    EXPECT_TRUE(c.is_on_curve());
    EXPECT_TRUE(c == G2Projective::generator());
}

TEST(TestG2, Add3) {
    G2Projective const a = G2Projective::identity();
    G2Projective b = G2Projective::generator();

    Fp2 const z{
            Fp({
                       0xba7afa1f9a6fe250, 0xfa0f5b595eafe731, 0x3bdc477694c306e7,
                       0x2149be4b3949fa24, 0x64aa6e0649b2078c, 0x12b108ac33643c3e,
               }),
            Fp({
                       0x125325df3d35b5a8, 0xdc469ef5555d7fe3, 0x02d716d2443106a9,
                       0x05a1db59a6ff37d0, 0x7cf7784e5300bb8f, 0x16a88922c7a5e844,
               }),

    };
    b = G2Projective{b.get_x() * z, b.get_y() * z, z};
    G2Projective const c = b + a;

    EXPECT_FALSE(c.is_identity());
    EXPECT_TRUE(c.is_on_curve());
    EXPECT_TRUE(c == G2Projective::generator());
}

TEST(TestG2, Add4) {
    G2Projective const a = G2Projective::generator().doubles().doubles();
    G2Projective const b = G2Projective::generator().doubles();
    G2Projective const c = a + b;
    G2Projective d = G2Projective::generator();
    for (int i = 0; i < 5; ++i) d += G2Projective::generator();

    EXPECT_FALSE(c.is_identity());
    EXPECT_TRUE(c.is_on_curve());
    EXPECT_FALSE(d.is_identity());
    EXPECT_TRUE(d.is_on_curve());
    EXPECT_EQ(c, d);
}

TEST(TestG2, Add5) {
    Fp2 beta{
            Fp({
                       0xcd03c9e48671f071, 0x5dab22461fcda5d2, 0x587042afd3851b95,
                       0x8eb60ebe01bacb9e, 0x03f97d6e83d050d2, 0x18f0206554638741,
               }),
            Fp::zero()
    };
    beta = beta.square();

    G2Projective const a = G2Projective::generator().doubles().doubles();
    G2Projective const b{a.get_x() * beta, -a.get_y(), a.get_z()};

    EXPECT_TRUE(a.is_on_curve());
    EXPECT_TRUE(b.is_on_curve());

    G2Projective const c = a + b;
    G2Projective const d{
            Fp2{
                    Fp({
                               0x705abc799ca773d3, 0xfe132292c1d4bf08, 0xf37ece3e07b2b466,
                               0x887e1c43f447e301, 0x1e0970d033bc77e8, 0x1985c81e20a693f2,
                       }),
                    Fp({
                               0x1d79b25db36ab924, 0x23948e4d529639d3, 0x471ba7fb0d006297,
                               0x2c36d4b4465dc4c0, 0x82bbc3cfec67f538, 0x051d2728b67bf952,
                       }),
            },
            Fp2{
                    Fp({
                               0x41b1bbf6576c0abf, 0xb6cc93713f7a0f9a, 0x6b65b43e48f3f01f,
                               0xfb7a4cfcaf81be4f, 0x3e32dadc6ec22cb6, 0x0bb0fc49d79807e3,
                       }),
                    Fp({
                               0x7d1397788f5f2ddf, 0xab2907144ff0d8e8, 0x5b7573e0cdb91f92,
                               0x4cb8932dd31daf28, 0x62bbfac6db052a54, 0x11f95c16d14c3bbe,
                       }),
            },
            Fp2::one(),
    };

    EXPECT_EQ(G2Affine(c), G2Affine(d));
    EXPECT_FALSE(c.is_identity());
    EXPECT_TRUE(c.is_on_curve());
}

TEST(TestG2, MixAdd1) {
    G2Affine const a = G2Affine::identity();
    G2Projective const b = G2Projective::identity();
    G2Projective const c = a + b;

    EXPECT_TRUE(c.is_identity());
    EXPECT_TRUE(c.is_on_curve());
}

TEST(TestG2, MixAdd2) {
    G2Affine const a = G2Affine::identity();
    G2Projective b = G2Projective::generator();

    Fp2 const z{
            Fp({
                       0xba7afa1f9a6fe250, 0xfa0f5b595eafe731, 0x3bdc477694c306e7,
                       0x2149be4b3949fa24, 0x64aa6e0649b2078c, 0x12b108ac33643c3e,
               }),
            Fp({
                       0x125325df3d35b5a8, 0xdc469ef5555d7fe3, 0x02d716d2443106a9,
                       0x05a1db59a6ff37d0, 0x7cf7784e5300bb8f, 0x16a88922c7a5e844,
               }),

    };
    b = G2Projective{b.get_x() * z, b.get_y() * z, z};
    G2Projective const c = a + b;

    EXPECT_FALSE(c.is_identity());
    EXPECT_TRUE(c.is_on_curve());
    EXPECT_TRUE(c == G2Projective::generator());
}

TEST(TestG2, MixAdd3) {
    G2Affine const a = G2Affine::identity();
    G2Projective b = G2Projective::generator();

    Fp2 const z{
            Fp({
                       0xba7afa1f9a6fe250, 0xfa0f5b595eafe731, 0x3bdc477694c306e7,
                       0x2149be4b3949fa24, 0x64aa6e0649b2078c, 0x12b108ac33643c3e,
               }),
            Fp({
                       0x125325df3d35b5a8, 0xdc469ef5555d7fe3, 0x02d716d2443106a9,
                       0x05a1db59a6ff37d0, 0x7cf7784e5300bb8f, 0x16a88922c7a5e844,
               }),

    };
    b = G2Projective{b.get_x() * z, b.get_y() * z, z};
    G2Projective const c = b + a;

    EXPECT_FALSE(c.is_identity());
    EXPECT_TRUE(c.is_on_curve());
    EXPECT_TRUE(c == G2Projective::generator());
}

TEST(TestG2, MixAdd4) {
    G2Projective const a = G2Projective::generator().doubles().doubles();
    G2Projective const b = G2Projective::generator().doubles();
    G2Projective const c = a + b;
    G2Projective d = G2Projective::generator();
    for (int i = 0; i < 5; ++i) d += G2Affine::generator();

    EXPECT_FALSE(c.is_identity());
    EXPECT_TRUE(c.is_on_curve());
    EXPECT_FALSE(d.is_identity());
    EXPECT_TRUE(d.is_on_curve());
    EXPECT_EQ(c, d);
}

TEST(TestG2, MixAdd5) {
    Fp2 beta{
            Fp({
                       0xcd03c9e48671f071, 0x5dab22461fcda5d2, 0x587042afd3851b95,
                       0x8eb60ebe01bacb9e, 0x03f97d6e83d050d2, 0x18f0206554638741,
               }),
            Fp::zero()
    };
    beta = beta.square();

    G2Projective const a = G2Projective::generator().doubles().doubles();
    G2Projective const b{a.get_x() * beta, -a.get_y(), a.get_z()};

    G2Affine const aa(a);

    EXPECT_TRUE(aa.is_on_curve());
    EXPECT_TRUE(b.is_on_curve());

    G2Projective const c = aa + b;
    G2Projective const d{
            Fp2{
                    Fp({
                               0x705abc799ca773d3, 0xfe132292c1d4bf08, 0xf37ece3e07b2b466,
                               0x887e1c43f447e301, 0x1e0970d033bc77e8, 0x1985c81e20a693f2,
                       }),
                    Fp({
                               0x1d79b25db36ab924, 0x23948e4d529639d3, 0x471ba7fb0d006297,
                               0x2c36d4b4465dc4c0, 0x82bbc3cfec67f538, 0x051d2728b67bf952,
                       }),
            },
            Fp2{
                    Fp({
                               0x41b1bbf6576c0abf, 0xb6cc93713f7a0f9a, 0x6b65b43e48f3f01f,
                               0xfb7a4cfcaf81be4f, 0x3e32dadc6ec22cb6, 0x0bb0fc49d79807e3,
                       }),
                    Fp({
                               0x7d1397788f5f2ddf, 0xab2907144ff0d8e8, 0x5b7573e0cdb91f92,
                               0x4cb8932dd31daf28, 0x62bbfac6db052a54, 0x11f95c16d14c3bbe,
                       }),
            },
            Fp2::one(),
    };

    EXPECT_EQ(G2Affine(c), G2Affine(d));
    EXPECT_FALSE(c.is_identity());
    EXPECT_TRUE(c.is_on_curve());
}

TEST(TestG2, ProjectiveNegSub) {
    G2Projective const a = G2Projective::generator().doubles();

    EXPECT_EQ(a + (-a), G2Projective::identity());
    EXPECT_EQ(a + (-a), a - a);
}

TEST(TestG2, AffineNegSub) {
    G2Affine const a = G2Affine::generator();

    EXPECT_EQ(G2Projective(a) +(-a), G2Projective::identity());
    EXPECT_EQ(G2Projective(a) +(-a), G2Projective(a) -a);
}

TEST(TestG2, AffineScalarMul) {
    G2Affine const gen = G2Affine::generator();
    Scalar const a = Scalar::from_raw(
            {0x2b568297a56da71c, 0xd8c39ecb0ef375d1, 0x435c38da67bfbf96, 0x8088a05026b659b2});
    Scalar const b = Scalar::from_raw(
            {0x785fdd9b26ef8b85, 0xc997f25837695c18, 0x4c8dbc39e7b756c1, 0x70d9b6cc6d87df20});
    Scalar const c = a * b;

    EXPECT_EQ(G2Affine(gen * a) * b, gen * c);
}

TEST(TestG2, TorsionFree) {
    G2Affine const a{
            Fp2{
                    Fp({
                               0x89f550c813db6431, 0xa50be8c456cd8a1a, 0xa45b374114cae851,
                               0xbb6190f5bf7fff63, 0x970ca02c3ba80bc7, 0x02b85d24e840fbac,
                       }),
                    Fp({
                               0x6888bc53d70716dc, 0x3dea6b4117682d70, 0xd8f5f930500ca354,
                               0x6b5ecb6556f5c155, 0xc96bef0434778ab0, 0x05081505515006ad,
                       }),
            },
            Fp2{
                    Fp({
                               0x3cf1ea0d434b0f40, 0x1a0dc610e603e333, 0x7f89956160c72fa0,
                               0x25ee03decf6431c5, 0xeee8e206ec0fe137, 0x097592b226dfef28,
                       }),
                    Fp({
                               0x71e8bb5f29247367, 0xa5fe049e211831ce, 0x0ce6b354502a3896,
                               0x93b012000997314e, 0x6759f3b6aa5b42ac, 0x156944c4dfe92bbb,
                       }),
            },
            false,
    };

    EXPECT_FALSE(a.is_torsion_free());
    EXPECT_TRUE(G2Affine::identity().is_torsion_free());
    EXPECT_TRUE(G2Affine::generator().is_torsion_free());
}

TEST(TestG2, MulByX) {
    G2Projective const gen = G2Projective::generator();
    Scalar const x = -Scalar(bls12_381::group::constant::BLS_X);

    auto ta = gen * x;
    auto tb = gen.mul_by_x();

    EXPECT_EQ(gen.mul_by_x(), gen * x);
    G2Projective const point = G2Projective::generator() * Scalar(42);
    EXPECT_EQ(point.mul_by_x(), point * x);
}

TEST(TestG2, Psi) {
    G2Projective const gen = G2Projective::generator();

    Fp2 const z{
            Fp({
                       0x0ef2ddffab187c0a, 0x2424522b7d5ecbfc, 0xc6f341a3398054f4,
                       0x5523ddf409502df0, 0xd55c0b5a88e0dd97, 0x066428d704923e52,
               }),
            Fp({
                       0x538bbe0c95b4878d, 0xad04a50379522881, 0x6d5c05bf5c12fb64,
                       0x4ce4a069a2d34787, 0x59ea6c8d0dffaeaf, 0x0d42a083a75bd6f3,
               }),
    };

    G2Projective const point{
            Fp2{
                    Fp({
                               0xee4c8cb7c047eaf2, 0x44ca22eee036b604, 0x33b3affb2aefe101,
                               0x15d3e45bbafaeb02, 0x7bfc2154cd7419a4, 0x0a2d0c2b756e5edc,
                       }),
                    Fp({
                               0xfc224361029a8777, 0x4cbf2baab8740924, 0xc5008c6ec6592c89,
                               0xecc2c57b472a9c2d, 0x8613eafd9d81ffb1, 0x10fe54daa2d3d495,
                       }),
            } * z,
            Fp2{
                    Fp({
                               0x7de7edc43953b75c, 0x58be1d2de35e87dc, 0x5731d30b0e337b40,
                               0xbe93b60cfeaae4c9, 0x8b22c203764bedca, 0x01616c8d1033b771,
                       }),
                    Fp({
                               0xea126fe476b5733b, 0x85cee68b5dae1652, 0x98247779f7272b04,
                               0xa649c8b468c6e808, 0xb5b9a62dff0c4e45, 0x1555b67fc7bbe73d,
                       }),
            },
            z.square() * z,
    };

    EXPECT_TRUE(point.is_on_curve());

    // psi2(P) = psi(psi(P))
    EXPECT_EQ(gen.psi2(), gen.psi().psi());
    EXPECT_EQ(point.psi2(), point.psi().psi());

    // psi(P) is a morphism
    EXPECT_EQ(gen.doubles().psi(), gen.psi().doubles());
    EXPECT_EQ(point.psi() + gen.psi(), (point + gen).psi());

    // psi(P) behaves in the same way on the same projective point
    auto normalized_point = G2Projective(point);
    EXPECT_EQ(point.psi(), normalized_point.psi());
    EXPECT_EQ(point.psi2(), normalized_point.psi2());
}

TEST(TestG2, CofactorClearance) {
    G2Projective const gen = G2Projective::generator();
    G2Projective const id = G2Projective::identity();

    EXPECT_TRUE(gen.clear_cofactor().is_on_curve());
    EXPECT_TRUE(id.clear_cofactor().is_on_curve());

    Fp2 const z{
            Fp({
                       0x0ef2ddffab187c0a, 0x2424522b7d5ecbfc, 0xc6f341a3398054f4,
                       0x5523ddf409502df0, 0xd55c0b5a88e0dd97, 0x066428d704923e52,
               }),
            Fp({
                       0x538bbe0c95b4878d, 0xad04a50379522881, 0x6d5c05bf5c12fb64,
                       0x4ce4a069a2d34787, 0x59ea6c8d0dffaeaf, 0x0d42a083a75bd6f3,
               }),
    };

    G2Projective const point{
            Fp2{
                    Fp({
                               0xee4c8cb7c047eaf2, 0x44ca22eee036b604, 0x33b3affb2aefe101,
                               0x15d3e45bbafaeb02, 0x7bfc2154cd7419a4, 0x0a2d0c2b756e5edc,
                       }),
                    Fp({
                               0xfc224361029a8777, 0x4cbf2baab8740924, 0xc5008c6ec6592c89,
                               0xecc2c57b472a9c2d, 0x8613eafd9d81ffb1, 0x10fe54daa2d3d495,
                       }),
            } * z,
            Fp2{
                    Fp({
                               0x7de7edc43953b75c, 0x58be1d2de35e87dc, 0x5731d30b0e337b40,
                               0xbe93b60cfeaae4c9, 0x8b22c203764bedca, 0x01616c8d1033b771,
                       }),
                    Fp({
                               0xea126fe476b5733b, 0x85cee68b5dae1652, 0x98247779f7272b04,
                               0xa649c8b468c6e808, 0xb5b9a62dff0c4e45, 0x1555b67fc7bbe73d,
                       }),
            },
            z.square() * z,
    };

    EXPECT_TRUE(point.is_on_curve());
    EXPECT_FALSE(G2Affine(point).is_torsion_free());

    G2Projective const cleared = point.clear_cofactor();

    EXPECT_TRUE(cleared.is_on_curve());
    EXPECT_TRUE(G2Affine(cleared).is_torsion_free());

    std::array<uint8_t, 32> const h_eff_mod_q = {
            0xff, 0xff, 0x01, 0x00, 0x04, 0x00, 0x02, 0xa4,
            0x09, 0x90, 0x06, 0x00, 0x04, 0x90, 0x16, 0xb1,
            0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    };
//    EXPECT_EQ(gen.clear_cofactor(), gen.multiply(h_eff_mod_q));
//    EXPECT_EQ(cleared.clear_cofactor(), cleared.multiply(h_eff_mod_q));
}

TEST(TestG2, BatchNormalize) {
    G2Projective const a = G2Projective::generator().doubles();
    G2Projective const b = a.doubles();
    G2Projective const c = b.doubles();

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k) {
                std::vector<G2Projective> v = {a, b, c};
                if (i == 1) v[0] = G2Projective::identity();
                if (j == 0) v[1] = G2Projective::identity();
                if (k == 0) v[2] = G2Projective::identity();

                std::vector<G2Affine> res = G2Projective::batch_normalize(v);
                G2Affine const expected[3] = {G2Affine(v[0]),
                                        G2Affine(v[1]),
                                        G2Affine(v[2])};

                EXPECT_EQ(res[0], expected[0]);
                EXPECT_EQ(res[1], expected[1]);
                EXPECT_EQ(res[2], expected[2]);
            }
        }
    }
}

TEST(TestG2, CommutativeScalarSubgroupMul) {
    Scalar const a = Scalar::from_raw({
                                        0x1fff3231233ffffd, 0x4884b7fa00034802,
                                        0x998c4fefecbc4ff3, 0x1824b159acc50562
                                });
    G2Affine const g1_a = G2Affine::generator();
    G2Projective const g1_p = G2Projective::generator();

    // commutative
    EXPECT_EQ(g1_a * a, a * g1_a);
    EXPECT_EQ(g1_p * a, a * g1_p);

    // mixed
    EXPECT_EQ(g1_a * a, a * g1_p);
    EXPECT_EQ(g1_p * a, a * g1_a);
}

TEST(G2, AffineBytesUnchecked) {
    const G2Affine gen = G2Affine::generator();
    const G2Affine id = G2Affine::identity();

    const std::array<uint8_t, G2Affine::RAW_SIZE> gen_bytes_array = gen.to_raw_bytes();

    std::vector<uint8_t> gen_bytes_vec{};
    gen_bytes_vec.reserve(G2Affine::RAW_SIZE);
    gen_bytes_vec.insert(gen_bytes_vec.end(), gen_bytes_array.begin(), gen_bytes_array.end());

    const G2Affine gen_recovered = G2Affine::from_slice_unchecked(gen_bytes_vec);

    const std::array<uint8_t, G2Affine::RAW_SIZE> id_bytes_array = id.to_raw_bytes();

    std::vector<uint8_t> id_bytes_vec{};
    id_bytes_vec.reserve(G2Affine::RAW_SIZE);
    id_bytes_vec.insert(id_bytes_vec.end(), id_bytes_array.begin(), id_bytes_array.end());

    const G2Affine id_recovered = G2Affine::from_slice_unchecked(id_bytes_vec);

    EXPECT_EQ(gen, gen_recovered);
    EXPECT_EQ(id, id_recovered);
}

TEST(G2, PreparedSerialization) {
    const G2Prepared g2_prepared{G2Affine::generator()};
    const auto bytes = g2_prepared.to_raw_bytes();
    const auto g2_prepared_recovered = G2Prepared::from_slice_unchecked(bytes);
    const auto bytes_2 = g2_prepared_recovered.to_raw_bytes();
    EXPECT_EQ(bytes, bytes_2);
}