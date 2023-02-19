#include <gtest/gtest.h>

#include "group/gt.h"
#include "pairing/pairing.h"
#include "scalar/scalar.h"

TEST(TestPairing, GtGenerator) {
    bls12_381::group::Gt p1 = bls12_381::group::Gt::generator();
    bls12_381::group::Gt p2 = bls12_381::pairing::pairings(bls12_381::group::G1Affine::generator(),
                                                           bls12_381::group::G2Affine::generator());

    EXPECT_EQ(p1, p2);
}

TEST(TestPairing, Bilinearity) {
    bls12_381::scalar::Scalar a = bls12_381::scalar::Scalar::from_raw({1, 2, 3, 4}).invert().value().square();
    bls12_381::scalar::Scalar b = bls12_381::scalar::Scalar::from_raw({5, 6, 7, 8}).invert().value().square();
    bls12_381::scalar::Scalar c = a * b;

    bls12_381::group::G1Affine g = bls12_381::group::G1Affine(bls12_381::group::G1Affine::generator() * a);
    bls12_381::group::G2Affine h = bls12_381::group::G2Affine(bls12_381::group::G2Affine::generator() * b);
    bls12_381::group::Gt p = bls12_381::pairing::pairings(g, h);

    EXPECT_TRUE(p != bls12_381::group::Gt::identity());

    bls12_381::group::G1Affine expected = bls12_381::group::G1Affine(bls12_381::group::G1Affine::generator() * c);

    EXPECT_EQ(p, bls12_381::pairing::pairings(expected, bls12_381::group::G2Affine::generator()));
    EXPECT_EQ(p, bls12_381::pairing::pairings(bls12_381::group::G1Affine::generator(),
                                              bls12_381::group::G2Affine::generator()) * c);
}

TEST(TestPairing, Unitary) {
    bls12_381::group::G1Affine g = bls12_381::group::G1Affine::generator();
    bls12_381::group::G2Affine h = bls12_381::group::G2Affine::generator();
    bls12_381::group::Gt p = -bls12_381::pairing::pairings(g, h);
    bls12_381::group::Gt q = bls12_381::pairing::pairings(g, -h);
    bls12_381::group::Gt r = bls12_381::pairing::pairings(-g, h);

    EXPECT_EQ(p, q);
    EXPECT_EQ(q, r);
}

TEST(TestPairing, MultiMillerLoop) {
    bls12_381::group::G1Affine a1 = bls12_381::group::G1Affine::generator();
    bls12_381::group::G2Affine b1 = bls12_381::group::G2Affine::generator();

    bls12_381::group::G1Affine a2 = bls12_381::group::G1Affine(bls12_381::group::G1Affine::generator() *
                                                               bls12_381::scalar::Scalar::from_raw(
                                                                       {1, 2, 3, 4}).invert().value().square());
    bls12_381::group::G2Affine b2 = bls12_381::group::G2Affine(bls12_381::group::G2Affine::generator() *
                                                               bls12_381::scalar::Scalar::from_raw(
                                                                       {4, 2, 2, 4}).invert().value().square());

    bls12_381::group::G1Affine a3 = bls12_381::group::G1Affine::identity();
    bls12_381::group::G2Affine b3 = bls12_381::group::G2Affine(bls12_381::group::G2Affine::generator() *
                                                               bls12_381::scalar::Scalar::from_raw(
                                                                       {9, 2, 2, 4}).invert().value().square());

    bls12_381::group::G1Affine a4 = bls12_381::group::G1Affine(bls12_381::group::G1Affine::generator() *
                                                               bls12_381::scalar::Scalar::from_raw(
                                                                       {5, 5, 5, 5}).invert().value().square());
    bls12_381::group::G2Affine b4 = bls12_381::group::G2Affine::identity();

    bls12_381::group::G1Affine a5 = bls12_381::group::G1Affine(bls12_381::group::G1Affine::generator() *
                                                               bls12_381::scalar::Scalar::from_raw(
                                                                       {323, 32, 3, 1}).invert().value().square());
    bls12_381::group::G2Affine b5 = bls12_381::group::G2Affine(bls12_381::group::G2Affine::generator() *
                                                               bls12_381::scalar::Scalar::from_raw(
                                                                       {4, 2, 2, 9099}).invert().value().square());

    bls12_381::group::G2Prepared b1_prepared = bls12_381::group::G2Prepared{b1};
    bls12_381::group::G2Prepared b2_prepared = bls12_381::group::G2Prepared{b2};
    bls12_381::group::G2Prepared b3_prepared = bls12_381::group::G2Prepared{b3};
    bls12_381::group::G2Prepared b4_prepared = bls12_381::group::G2Prepared{b4};
    bls12_381::group::G2Prepared b5_prepared = bls12_381::group::G2Prepared{b5};

    bls12_381::group::Gt expected = bls12_381::pairing::pairings(a1, b1)
                                    + bls12_381::pairing::pairings(a2, b2)
                                    + bls12_381::pairing::pairings(a3, b3)
                                    + bls12_381::pairing::pairings(a4, b4)
                                    + bls12_381::pairing::pairings(a5, b5);

    bls12_381::group::Gt test = bls12_381::pairing::multi_miller_loop({
                                                                              {a1, b1_prepared},
                                                                              {a2, b2_prepared},
                                                                              {a3, b3_prepared},
                                                                              {a4, b4_prepared},
                                                                              {a5, b5_prepared},
                                                                      }).final_exponentiation();

    EXPECT_EQ(expected, test);
}

TEST(TestPairing, MillerLoopResult) {
    EXPECT_EQ(bls12_381::pairing::MillerLoopResult{}.final_exponentiation(), bls12_381::group::Gt::identity());
}

TEST(TestPairing, TrickingMillerLoopResult) {
    bls12_381::field::Fp12 data0 =
            bls12_381::pairing::multi_miller_loop({
                                                          {bls12_381::group::G1Affine::identity(),
                                                           bls12_381::group::G2Prepared(
                                                                   bls12_381::group::G2Affine::generator())}
                                                  }).get_data();
    bls12_381::field::Fp12 data1 =
            bls12_381::pairing::multi_miller_loop({
                                                          {bls12_381::group::G1Affine::generator(),
                                                           bls12_381::group::G2Prepared(
                                                                   bls12_381::group::G2Affine::identity())}
                                                  }).get_data();
    bls12_381::field::Fp12 data2 =
            bls12_381::pairing::multi_miller_loop({
                                                          {bls12_381::group::G1Affine::generator(),
                                                           bls12_381::group::G2Prepared(
                                                                   bls12_381::group::G2Affine::generator())},
                                                          {-bls12_381::group::G1Affine::generator(),
                                                           bls12_381::group::G2Prepared(
                                                                   bls12_381::group::G2Affine::generator())}
                                                  }).get_data();

    EXPECT_EQ(data0, bls12_381::field::Fp12::one());
    EXPECT_EQ(data1, bls12_381::field::Fp12::one());
    EXPECT_NE(data2, bls12_381::field::Fp12::one());

    bls12_381::group::Gt point =
            bls12_381::pairing::multi_miller_loop({
                                                          {bls12_381::group::G1Affine::generator(),
                                                                  bls12_381::group::G2Prepared(
                                                                          bls12_381::group::G2Affine::generator())},
                                                          {-bls12_381::group::G1Affine::generator(),
                                                                  bls12_381::group::G2Prepared(
                                                                          bls12_381::group::G2Affine::generator())}
                                                  }).final_exponentiation();

    EXPECT_EQ(point, bls12_381::group::Gt::identity());
}