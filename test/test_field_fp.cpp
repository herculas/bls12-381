#include <gtest/gtest.h>
#include "field/fp.h"

using bls12_381::field::Fp;

TEST(TestFp, Mul) {
    Fp a({
                 0x0397a38320170cd4, 0x734c1b2c9e761d30, 0x5ed255ad9a48beb5,
                 0x095a3c6b22a7fcfc, 0x2294ce75d4e26a27, 0x13338bd870011ebb,
         });
    Fp b({
                 0xb9c3c7c5b1196af7, 0x2580e2086ce335c1, 0xf49aed3d8a57ef42,
                 0x41f281e49846e878, 0xe0762346c38452ce, 0x0652e89326e57dc0,
         });
    Fp c({
                 0xf96ef3d711ab5355, 0xe8d459ea00f148dd, 0x53f7354a5f00fa78,
                 0x9e34a4f3125c5f83, 0x3fbe0c47ca74c19e, 0x01b06a8bbd4adfe4,
         });

    EXPECT_EQ(a *= b, c);
}

TEST(TestFp, FromBytes) {
    std::array<uint8_t, Fp::WIDTH * sizeof(uint64_t)> array1 = {
            26, 1, 17, 234, 57, 127, 230, 154,
            75, 27, 167, 182, 67, 75, 172, 215,
            100, 119, 75, 132, 243, 133, 18, 191,
            103, 48, 210, 160, 246, 176, 246, 36,
            30, 171, 255, 254, 177, 83, 255, 255,
            185, 254, 255, 255, 255, 255, 170, 170
    };

    std::array<uint8_t, Fp::WIDTH * sizeof(uint64_t)> array2 = {
            27, 1, 17, 234, 57, 127, 230, 154,
            75, 27, 167, 182, 67, 75, 172, 215,
            100, 119, 75, 132, 243, 133, 18, 191,
            103, 48, 210, 160, 246, 176, 246, 36,
            30, 171, 255, 254, 177, 83, 255, 255,
            185, 254, 255, 255, 255, 255, 170, 170
    };

    Fp a = -Fp::one();
    Fp b = Fp::from_bytes(array1).value();
    EXPECT_EQ(a, b);

    std::optional<Fp> c = Fp::from_bytes(array2);
    EXPECT_FALSE(c.has_value());
}

TEST(TestFp, ToBytes) {
    Fp a({
                 0xdc906d9be3f95dc8, 0x8755caf7459691a1, 0xcff1a7f4e9583ab3,
                 0x9b43821f849e2284, 0xf57554f3a2974f3f, 0x085dbea84ed47f79,
         });

    std::array<uint8_t, Fp::WIDTH * sizeof(uint64_t)> bytes{};
    Fp b;

    for (int i = 0; i < 100; ++i) {
        a = a.square();
        bytes = a.to_bytes();
        b = Fp::from_bytes(bytes).value();
        EXPECT_EQ(a, b);
    }
}

TEST(TestFp, Square) {
    Fp a({
                 0xd215d2768e83191b, 0x5085d80f8fb28261, 0xce9a032ddf393a56,
                 0x3e9c4fff2ca0c4bb, 0x6436b6f7f4d95dfb, 0x10606628ad4a4d90,
         });
    Fp b({
                 0x33d9c42a3cb3e235, 0xdad11a094c4cd455, 0xa2f144bd729aaeba,
                 0xd4150932be9ffeac, 0xe27bc7c47d44ee50, 0x14b6a78d3ec7a560,
         });

    EXPECT_EQ(a.square(), b);
}

TEST(TestFp, Sqrt) {
    // a = 4.
    Fp a({
                 0xaa270000000cfff3, 0x53cc0032fc34000a, 0x478fe97a6b0a807f,
                 0xb1d37ebee6ba24d7, 0x8ec9733bbf78ab2f, 0x09d645513d83de7e,
         });

    // b = 2.
    Fp b({
                 0x321300000006554f, 0xb93c0018d6c40005, 0x57605e0db0ddbb51,
                 0x8b256521ed1f9bcb, 0x6cf28d7901622c03, 0x11ebab9dbb81e28c,
         });

    // sqrt(4) = -2.
    EXPECT_EQ(a.sqrt().value(), -b);
}

TEST(TestFp, Inverse) {
    Fp a({
                 0x43b43a5078ac2076, 0x1ce0763046f8962b, 0x724a5276486d735c,
                 0x6f05c2a6282d48fd, 0x2095bd5bb4ca9331, 0x03b35b3894b0f7da,
         });
    Fp b({
                 0x69ecd7040952148f, 0x985ccc2022190f55, 0xe19bba36a9ad2f41,
                 0x19bb16c95219dbd8, 0x14dcacfdfb478693, 0x115ff58afff9a8e1,
         });
    EXPECT_EQ(a.invert().value(), b);
}

TEST(TestFp, Lexical) {
    Fp a({
                 0xa1fafffffffe5557, 0x995bfff976a3fffe, 0x03f41d24d174ceb4,
                 0xf6547998c1995dbd, 0x778a468f507a6034, 0x020559931f7f8103
         });
    Fp b({
                 0x1804000000015554, 0x855000053ab00001, 0x633cb57c253c276f,
                 0x6e22d1ec31ebb502, 0xd3916126f2d14ca2, 0x17fbb8571a006596,
         });
    Fp c({
                 0x43f5fffffffcaaae, 0x32b7fff2ed47fffd, 0x07e83a49a2e99d69,
                 0xeca8f3318332bb7a, 0xef148d1ea0f4c069, 0x040ab3263eff0206,
         });

    EXPECT_FALSE(Fp::zero().lexicographically_largest());
    EXPECT_FALSE(Fp::one().lexicographically_largest());
    EXPECT_FALSE(a.lexicographically_largest());
    EXPECT_TRUE(b.lexicographically_largest());
    EXPECT_TRUE(c.lexicographically_largest());
}

TEST(TestFp, Str) {
    Fp a({
                 0x5360bb5978678032, 0x7dd275ae799e128e, 0x5c5b5071ce4f4dcf,
                 0xcdb21f93078dbb3e, 0xc32365c5e73f474a, 0x115a2a5489babe5b,
         });
    EXPECT_TRUE(
            a.to_hex_str() ==
            "0x104bf052ad3bc99bcb176c24a06a6c3aad4eaf2308fc4d282e106c84a757d061052630515305e59bdddf8111bfdeb704"
    );
}