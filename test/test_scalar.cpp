#include <gtest/gtest.h>
#include "scalar/constant.h"
#include "scalar/scalar.h"

using bls12_381::scalar::Scalar;

TEST(TestScalar, Inv) {
    uint64_t inv = 1;
    for (int i = 0; i < 63; ++i) {
        inv = inv * inv;
        inv = inv * 0xffffffff00000001;
    }
    inv = -inv;
    EXPECT_EQ(inv, bls12_381::scalar::constant::INV);
}

TEST(TestScalar, Zero) {
    EXPECT_EQ(Scalar::zero(), -Scalar::zero());
    EXPECT_EQ(Scalar::zero(), Scalar::zero() + Scalar::zero());
    EXPECT_EQ(Scalar::zero(), Scalar::zero() - Scalar::zero());
    EXPECT_EQ(Scalar::zero(), Scalar::zero() * Scalar::zero());
}

TEST(TestScalar, HexStr) {
    EXPECT_TRUE(Scalar::zero().get_hex() == "0x0000000000000000000000000000000000000000000000000000000000000000");
    EXPECT_TRUE(Scalar::one().get_hex() == "0x0000000000000000000000000000000000000000000000000000000000000001");
    EXPECT_TRUE(bls12_381::scalar::constant::R2.get_hex() ==
                "0x1824b159acc5056f998c4fefecbc4ff55884b7fa0003480200000001fffffffe");
}

TEST(TestScalar, Equality) {
    EXPECT_EQ(Scalar::zero(), Scalar::zero());
    EXPECT_EQ(Scalar::one(), Scalar::one());
    EXPECT_EQ(bls12_381::scalar::constant::R2, bls12_381::scalar::constant::R2);

    EXPECT_NE(Scalar::zero(), Scalar::one());
    EXPECT_NE(Scalar::one(), bls12_381::scalar::constant::R2);
}

TEST(TestScalar, ToBytes) {
    std::array<uint8_t, Scalar::WIDTH * sizeof(uint64_t)> expected_0 = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };
    std::array<uint8_t, Scalar::WIDTH * sizeof(uint64_t)> expected_1 = {
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };
    std::array<uint8_t, Scalar::WIDTH * sizeof(uint64_t)> expected_r2 = {
            254, 255, 255, 255, 1, 0, 0, 0, 2, 72, 3, 0, 250, 183, 132, 88,
            245, 79, 188, 236, 239, 79, 140, 153, 111, 5, 197, 172, 89, 177, 36, 24
    };
    std::array<uint8_t, Scalar::WIDTH * sizeof(uint64_t)> expected_neg_1 = {
            0, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83,
            5, 216, 161, 9, 8, 216, 57, 51, 72, 125, 157, 41, 83, 167, 237, 115
    };

    std::array<uint8_t, Scalar::WIDTH * sizeof(uint64_t)> bytes_0 = Scalar::zero().to_bytes();
    std::array<uint8_t, Scalar::WIDTH * sizeof(uint64_t)> bytes_1 = Scalar::one().to_bytes();
    std::array<uint8_t, Scalar::WIDTH * sizeof(uint64_t)> bytes_r2 = bls12_381::scalar::constant::R2.to_bytes();
    std::array<uint8_t, Scalar::WIDTH * sizeof(uint64_t)> bytes_neg_1 = (-Scalar::one()).to_bytes();

    EXPECT_TRUE(bytes_0 == expected_0);
    EXPECT_TRUE(bytes_1 == expected_1);
    EXPECT_TRUE(bytes_r2 == expected_r2);
    EXPECT_TRUE(bytes_neg_1 == expected_neg_1);
}

TEST(TestScalar, FromBytes) {
    std::array<uint8_t, Scalar::WIDTH * sizeof(uint64_t)> bytes_0 = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };
    std::array<uint8_t, Scalar::WIDTH * sizeof(uint64_t)> bytes_1 = {
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };
    std::array<uint8_t, Scalar::WIDTH * sizeof(uint64_t)> bytes_r2 = {
            254, 255, 255, 255, 1, 0, 0, 0, 2, 72, 3, 0, 250, 183, 132, 88,
            245, 79, 188, 236, 239, 79, 140, 153, 111, 5, 197, 172, 89, 177, 36, 24
    };
    std::array<uint8_t, Scalar::WIDTH * sizeof(uint64_t)> bytes_neg_1 = {
            0, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83,
            5, 216, 161, 9, 8, 216, 57, 51, 72, 125, 157, 41, 83, 167, 237, 115
    };
    std::array<uint8_t, Scalar::WIDTH * sizeof(uint64_t)> modulus = {
            1, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83,
            5, 216, 161, 9, 8, 216, 57, 51, 72, 125, 157, 41, 83, 167, 237, 115
    };
    std::array<uint8_t, Scalar::WIDTH * sizeof(uint64_t)> larger_than_modulus_0 = {
            2, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83,
            5, 216, 161, 9, 8, 216, 57, 51, 72, 125, 157, 41, 83, 167, 237, 115
    };
    std::array<uint8_t, Scalar::WIDTH * sizeof(uint64_t)> larger_than_modulus_1 = {
            1, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83,
            5, 216, 161, 9, 8, 216, 58, 51, 72, 125, 157, 41, 83, 167, 237, 115
    };
    std::array<uint8_t, Scalar::WIDTH * sizeof(uint64_t)> larger_than_modulus_2 = {
            1, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83,
            5, 216, 161, 9, 8, 216, 57, 51, 72, 125, 157, 41, 83, 167, 237, 116
    };

    Scalar expected_0 = Scalar::from_bytes(bytes_0).value();
    Scalar expected_1 = Scalar::from_bytes(bytes_1).value();
    Scalar expected_r2 = Scalar::from_bytes(bytes_r2).value();
    Scalar expected_neg_1 = Scalar::from_bytes(bytes_neg_1).value();

    EXPECT_TRUE(expected_0 == Scalar::zero());
    EXPECT_TRUE(expected_1 == Scalar::one());
    EXPECT_TRUE(expected_r2 == bls12_381::scalar::constant::R2);
    EXPECT_TRUE(expected_neg_1 == -Scalar::one());

    EXPECT_FALSE(Scalar::from_bytes(modulus).has_value());
    EXPECT_FALSE(Scalar::from_bytes(larger_than_modulus_0).has_value());
    EXPECT_FALSE(Scalar::from_bytes(larger_than_modulus_1).has_value());
    EXPECT_FALSE(Scalar::from_bytes(larger_than_modulus_2).has_value());
}

TEST(TestScalar, FromWide) {
    std::array<uint8_t, Scalar::WIDTH * sizeof(uint64_t) * 2> bytes_r2 = {
            254, 255, 255, 255, 1, 0, 0, 0, 2, 72, 3, 0, 250, 183, 132, 88,
            245, 79, 188, 236, 239, 79, 140, 153, 111, 5, 197, 172, 89, 177, 36, 24,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    };
    std::array<uint8_t, Scalar::WIDTH * sizeof(uint64_t) * 2> bytes_n1 = {
            0, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83,
            5, 216, 161, 9, 8, 216, 57, 51, 72, 125, 157, 41, 83, 167, 237, 115,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    };
    std::array<uint8_t, Scalar::WIDTH * sizeof(uint64_t) * 2> bytes_max = {
            0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
            0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
            0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
            0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    };

    Scalar max_exp({0xc62c1805439b73b1, 0xc2b9551e8ced218e, 0xda44ec81daf9a422, 0x5605aa601c162e79});

    Scalar from_r2 = Scalar::from_bytes_wide(bytes_r2);
    Scalar from_n1 = Scalar::from_bytes_wide(bytes_n1);
    Scalar from_max = Scalar::from_bytes_wide(bytes_max);

    EXPECT_EQ(from_r2, bls12_381::scalar::constant::R2);
    EXPECT_EQ(from_n1, -Scalar::one());
    EXPECT_EQ(from_max, max_exp);
}

TEST(TestScalar, FromRaw) {
    std::array<uint64_t, Scalar::WIDTH> raw_1 = {
            0x00000001fffffffd, 0x5884b7fa00034802,
            0x998c4fefecbc4ff5, 0x1824b159acc5056f,
    };
    std::array<uint64_t, Scalar::WIDTH> raw_2 = {
            0xffffffffffffffff, 0xffffffffffffffff,
            0xffffffffffffffff, 0xffffffffffffffff,
    };
    std::array<uint64_t, Scalar::WIDTH> modulus = {
            0xffffffff00000001, 0x53bda402fffe5bfe,
            0x3339d80809a1d805, 0x73eda753299d7d48,
    };

    EXPECT_EQ(Scalar::from_raw(raw_1), Scalar::from_raw(raw_2));
    EXPECT_EQ(Scalar::from_raw(modulus), Scalar::zero());
    EXPECT_EQ(Scalar::from_raw({1, 0, 0, 0}), bls12_381::scalar::constant::R1);
}

TEST(TestScalar, Neg) {
    Scalar largest({0xffffffff00000000, 0x53bda402fffe5bfe, 0x3339d80809a1d805, 0x73eda753299d7d48});
    Scalar exp_neg_largest = -largest;
    Scalar neg_largest({1, 0, 0, 0});

    EXPECT_EQ(exp_neg_largest, neg_largest);

    Scalar neg_zero = -Scalar::zero();
    EXPECT_EQ(neg_zero, Scalar::zero());

    Scalar neg_neg_largest = -neg_largest;
    EXPECT_EQ(neg_neg_largest, largest);
}

TEST(TestScalar, Add) {
    Scalar largest({0xffffffff00000000, 0x53bda402fffe5bfe, 0x3339d80809a1d805, 0x73eda753299d7d48});
    Scalar largest2({0xffffffff00000000, 0x53bda402fffe5bfe, 0x3339d80809a1d805, 0x73eda753299d7d48});
    Scalar expected({0xfffffffeffffffff, 0x53bda402fffe5bfe, 0x3339d80809a1d805, 0x73eda753299d7d48});
    largest += largest;
    largest2 += Scalar({1, 0, 0, 0});

    EXPECT_EQ(largest, expected);
    EXPECT_EQ(largest2, Scalar::zero());
}

TEST(TestScalar, Sub) {
    const Scalar largest({0xffffffff00000000, 0x53bda402fffe5bfe, 0x3339d80809a1d805, 0x73eda753299d7d48});
    Scalar temp = largest;
    temp -= largest;
    EXPECT_EQ(temp, Scalar::zero());

    temp = Scalar::zero();
    temp -= largest;

    Scalar temp2 = bls12_381::scalar::constant::MODULUS;
    temp2 -= largest;
    EXPECT_EQ(temp, temp2);
}

TEST(TestScalar, Mul) {
    const Scalar largest({0xffffffff00000000, 0x53bda402fffe5bfe, 0x3339d80809a1d805, 0x73eda753299d7d48});
    Scalar current = largest;

    for (int i = 0; i < 100; ++i) {
        Scalar temp = current;
        temp *= current;
        Scalar temp2 = Scalar::zero();

        std::array<uint8_t, 32> bytes = current.to_bytes();

        for (auto iter = bytes.rbegin(); iter != bytes.rend(); ++iter) { // NOLINT(modernize-loop-convert)
            for (int k = 7; k >= 0; --k) {
                bool b = ((*iter >> k) & static_cast<uint8_t>(1)) == static_cast<uint8_t>(1);
                Scalar temp3 = temp2;
                temp2 += temp3;
                if (b) temp2 += current;
            }
        }
        EXPECT_EQ(temp, temp2);
        current += largest;
    }
}

TEST(TestScalar, Squaring) {
    const Scalar largest({0xffffffff00000000, 0x53bda402fffe5bfe, 0x3339d80809a1d805, 0x73eda753299d7d48});
    Scalar current = largest;

    for (int i = 0; i < 100; ++i) {
        Scalar temp = current;
        temp = temp.square();
        Scalar temp2 = Scalar::zero();

        std::array<uint8_t, 32> bytes = current.to_bytes();

        for (auto iter = bytes.rbegin(); iter != bytes.rend(); ++iter) { // NOLINT(modernize-loop-convert)
            for (int k = 7; k >= 0; --k) {
                bool b = ((*iter >> k) & static_cast<uint8_t>(1)) == static_cast<uint8_t>(1);
                Scalar temp3 = temp2;
                temp2 += temp3;
                if (b) temp2 += current;
            }
        }
        EXPECT_EQ(temp, temp2);
        current += largest;
    }
}

TEST(TestScalar, Inversion) {
    EXPECT_FALSE(Scalar::zero().invert().has_value());
    EXPECT_EQ(Scalar::one().invert().value(), Scalar::one());
    EXPECT_EQ((-Scalar::one()).invert().value(), -Scalar::one());

    Scalar temp = bls12_381::scalar::constant::R2;
    for (int i = 0; i < 100; ++i) {
        Scalar temp2 = temp.invert().value();
        temp2 *= temp;
        EXPECT_EQ(temp2, Scalar::one());
        temp += bls12_381::scalar::constant::R2;
    }
}

TEST(TestScalar, InvertIsPow) {
    std::array<uint64_t, Scalar::WIDTH> q_minus_2 = {
            0xfffffffeffffffff, 0x53bda402fffe5bfe,
            0x3339d80809a1d805, 0x73eda753299d7d48,
    };

    Scalar r1 = bls12_381::scalar::constant::R1;
    Scalar r2 = bls12_381::scalar::constant::R1;

    for (int i = 0; i < 100; ++i) {
        r1 = r1.invert().value();
        r2 = r2.pow(q_minus_2);

        EXPECT_EQ(r1, r2);

        r1 += bls12_381::scalar::constant::R1;
        r2 = r1;
    }
}

TEST(TestScalar, Sqrt) {
    EXPECT_EQ(Scalar::zero().sqrt().value(), Scalar::zero());

    Scalar square({0x46cd85a5f273077e, 0x1d30c47dd68fc735, 0x77f656f60beca0eb, 0x494aa01bdf32468d});
    int32_t none_count = 0;
    for (int i = 0; i < 100; ++i) {
        std::optional<Scalar> square_root = square.sqrt();
        if (!square_root.has_value()) {
            none_count += 1;
        } else {
            EXPECT_EQ(square_root.value() * square_root.value(), square);
        }
        square -= Scalar::one();
    }
    EXPECT_EQ(49, none_count);
}

TEST(TestScalar, Double) {
    Scalar a = Scalar::from_raw({0x1fff3231233ffffd, 0x4884b7fa00034802, 0x998c4fefecbc4ff3, 0x1824b159acc50562});
    EXPECT_EQ(a.doubles(), a + a);
}
