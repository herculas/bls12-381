#include <gtest/gtest.h>

#include "group/gt.h"
#include "pairing/pairing.h"
#include "scalar/scalar.h"

using bls12_381::field::Fp12;
using bls12_381::group::G1Affine;
using bls12_381::group::G2Affine;
using bls12_381::group::G2Prepared;
using bls12_381::group::Gt;
using bls12_381::scalar::Scalar;
using bls12_381::pairing::MillerLoopResult;
using bls12_381::pairing::pairings;
using bls12_381::pairing::multi_miller_loop;

TEST(TestPairing, GtGenerator) {
    Gt p1 = Gt::generator();
    Gt p2 = pairings(G1Affine::generator(), G2Affine::generator());

    EXPECT_EQ(p1, p2);
}

TEST(TestPairing, Bilinearity) {
    Scalar a = Scalar::from_raw({1, 2, 3, 4}).invert().value().square();
    Scalar b = Scalar::from_raw({5, 6, 7, 8}).invert().value().square();
    Scalar c = a * b;

    G1Affine g = G1Affine(G1Affine::generator() * a);
    G2Affine h = G2Affine(G2Affine::generator() * b);
    Gt p = pairings(g, h);

    EXPECT_TRUE(p != Gt::identity());

    G1Affine expected = G1Affine(G1Affine::generator() * c);

    EXPECT_EQ(p, pairings(expected, G2Affine::generator()));
    EXPECT_EQ(p, pairings(G1Affine::generator(), G2Affine::generator()) * c);
}

TEST(TestPairing, Unitary) {
    G1Affine g = G1Affine::generator();
    G2Affine h = G2Affine::generator();
    Gt p = -pairings(g, h);
    Gt q = pairings(g, -h);
    Gt r = pairings(-g, h);

    EXPECT_EQ(p, q);
    EXPECT_EQ(q, r);
}

TEST(TestPairing, MultiMillerLoop) {
    G1Affine a1 = G1Affine::generator();
    G2Affine b1 = G2Affine::generator();

    G1Affine a2 = G1Affine(G1Affine::generator() * Scalar::from_raw({1, 2, 3, 4}).invert().value().square());
    G2Affine b2 = G2Affine(G2Affine::generator() * Scalar::from_raw({4, 2, 2, 4}).invert().value().square());

    G1Affine a3 = G1Affine::identity();
    G2Affine b3 = G2Affine(G2Affine::generator() * Scalar::from_raw({9, 2, 2, 4}).invert().value().square());

    G1Affine a4 = G1Affine(G1Affine::generator() * Scalar::from_raw({5, 5, 5, 5}).invert().value().square());
    G2Affine b4 = G2Affine::identity();

    G1Affine a5 = G1Affine(G1Affine::generator() * Scalar::from_raw({323, 32, 3, 1}).invert().value().square());
    G2Affine b5 = G2Affine(G2Affine::generator() * Scalar::from_raw({4, 2, 2, 9099}).invert().value().square());

    G2Prepared b1_prepared = G2Prepared{b1};
    G2Prepared b2_prepared = G2Prepared{b2};
    G2Prepared b3_prepared = G2Prepared{b3};
    G2Prepared b4_prepared = G2Prepared{b4};
    G2Prepared b5_prepared = G2Prepared{b5};

    Gt expected = pairings(a1, b1) + pairings(a2, b2) + pairings(a3, b3) + pairings(a4, b4) + pairings(a5, b5);

    Gt test = multi_miller_loop({
                                        {a1, b1_prepared},
                                        {a2, b2_prepared},
                                        {a3, b3_prepared},
                                        {a4, b4_prepared},
                                        {a5, b5_prepared},
                                })
            .final_exponentiation();

    EXPECT_EQ(expected, test);
}

TEST(TestPairing, MillerLoopResult) {
    EXPECT_EQ(MillerLoopResult{}.final_exponentiation(), Gt::identity());
}

TEST(TestPairing, TrickingMillerLoopResult) {
    Fp12 data0 = multi_miller_loop({{G1Affine::identity(), G2Prepared(G2Affine::generator())}}).get_data();
    Fp12 data1 = multi_miller_loop({{G1Affine::generator(), G2Prepared(G2Affine::identity())}}).get_data();
    Fp12 data2 = multi_miller_loop({{G1Affine::generator(),  G2Prepared(G2Affine::generator())},
                                    {-G1Affine::generator(), G2Prepared(G2Affine::generator())}})
            .get_data();

    EXPECT_EQ(data0, Fp12::one());
    EXPECT_EQ(data1, Fp12::one());
    EXPECT_NE(data2, Fp12::one());

    Gt point = multi_miller_loop({{G1Affine::generator(),  G2Prepared(G2Affine::generator())},
                                  {-G1Affine::generator(), G2Prepared(G2Affine::generator())}})
            .final_exponentiation();

    EXPECT_EQ(point, Gt::identity());
}