#include <gtest/gtest.h>
#include "field/fp12.h"

TEST(TestFp12, Arithmetic) {
    Fp12 a{
            Fp6{
                    Fp2{
                            Fp({
                                       0x47f9cb98b1b82d58, 0x5fe911eba3aa1d9d, 0x96bf1b5f4dd81db3,
                                       0x8100d27cc9259f5b, 0xafa20b9674640eab, 0x09bbcea7d8d9497d,
                               }),
                            Fp({
                                       0x0303cb98b1662daa, 0xd93110aa0a621d5a, 0xbfa9820c5be4a468,
                                       0x0ba3643ecb05a348, 0xdc3534bb1f1c25a6, 0x06c305bb19c0e1c1,
                               }),
                    },
                    Fp2{
                            Fp({
                                       0x46f9cb98b162d858, 0x0be9109cf7aa1d57, 0xc791bc55fece41d2,
                                       0xf84c57704e385ec2, 0xcb49c1d9c010e60f, 0x0acdb8e158bfe3c8,
                               }),
                            Fp({
                                       0x8aefcb98b15f8306, 0x3ea1108fe4f21d54, 0xcf79f69fa1b7df3b,
                                       0xe4f54aa1d16b1a3c, 0xba5e4ef86105a679, 0x0ed86c0797bee5cf,
                               }),
                    },
                    Fp2{
                            Fp({
                                       0xcee5cb98b15c2db4, 0x71591082d23a1d51, 0xd76230e944a17ca4,
                                       0xd19e3dd3549dd5b6, 0xa972dc1701fa66e3, 0x12e31f2dd6bde7d6,
                               }),
                            Fp({
                                       0xad2acb98b1732d9d, 0x2cfd10dd06961d64, 0x07396b86c6ef24e8,
                                       0xbd76e2fdb1bfc820, 0x6afea7f6de94d0d5, 0x10994b0c5744c040,
                               }),
                    },
            },
            Fp6{
                    Fp2{
                            Fp({
                                       0x47f9cb98b1b82d58, 0x5fe911eba3aa1d9d, 0x96bf1b5f4dd81db3,
                                       0x8100d27cc9259f5b, 0xafa20b9674640eab, 0x09bbcea7d8d9497d,
                               }),
                            Fp({
                                       0x0303cb98b1662daa, 0xd93110aa0a621d5a, 0xbfa9820c5be4a468,
                                       0x0ba3643ecb05a348, 0xdc3534bb1f1c25a6, 0x06c305bb19c0e1c1,
                               }),
                    },
                    Fp2{
                            Fp({
                                       0x46f9cb98b162d858, 0x0be9109cf7aa1d57, 0xc791bc55fece41d2,
                                       0xf84c57704e385ec2, 0xcb49c1d9c010e60f, 0x0acdb8e158bfe3c8,
                               }),
                            Fp({
                                       0x8aefcb98b15f8306, 0x3ea1108fe4f21d54, 0xcf79f69fa1b7df3b,
                                       0xe4f54aa1d16b1a3c, 0xba5e4ef86105a679, 0x0ed86c0797bee5cf,
                               }),
                    },
                    Fp2{
                            Fp({
                                       0xcee5cb98b15c2db4, 0x71591082d23a1d51, 0xd76230e944a17ca4,
                                       0xd19e3dd3549dd5b6, 0xa972dc1701fa66e3, 0x12e31f2dd6bde7d6,
                               }),
                            Fp({
                                       0xad2acb98b1732d9d, 0x2cfd10dd06961d64, 0x07396b86c6ef24e8,
                                       0xbd76e2fdb1bfc820, 0x6afea7f6de94d0d5, 0x10994b0c5744c040,
                               }),
                    },
            },
    };

    Fp12 b{
            Fp6{
                    Fp2{
                            Fp({
                                       0x47f9cb98b1b82d58, 0x5fe911eba3aa1d9d, 0x96bf1b5f4dd81db3,
                                       0x8100d272c9259f5b, 0xafa20b9674640eab, 0x09bbcea7d8d9497d,
                               }),
                            Fp({
                                       0x0303cb98b1662daa, 0xd93110aa0a621d5a, 0xbfa9820c5be4a468,
                                       0x0ba3643ecb05a348, 0xdc3534bb1f1c25a6, 0x06c305bb19c0e1c1,
                               }),
                    },
                    Fp2{
                            Fp({
                                       0x46f9cb98b162d858, 0x0be9109cf7aa1d57, 0xc791bc55fece41d2,
                                       0xf84c57704e385ec2, 0xcb49c1d9c010e60f, 0x0acdb8e158bfe348,
                               }),
                            Fp({
                                       0x8aefcb98b15f8306, 0x3ea1108fe4f21d54, 0xcf79f69fa1b7df3b,
                                       0xe4f54aa1d16b1a3c, 0xba5e4ef86105a679, 0x0ed86c0797bee5cf,
                               }),
                    },
                    Fp2{
                            Fp({
                                       0xcee5cb98b15c2db4, 0x71591082d23a1d51, 0xd76230e944a17ca4,
                                       0xd19e3dd3549dd5b6, 0xa972dc1701fa66e3, 0x12e31f2dd6bde7d6,
                               }),
                            Fp({
                                       0xad2acb98b1732d9d, 0x2cfd10dd06961d64, 0x07396b86c6ef24e8,
                                       0xbd76e2fdb1bfc820, 0x6afea7f6de94d0d5, 0x10994b0c5744c040,
                               }),
                    },
            },
            Fp6{
                    Fp2{
                            Fp({
                                       0x47f9cb98b1b82d58, 0x5fe911eba3aa1d9d, 0x96bf1b5f4dd21db3,
                                       0x8100d27cc9259f5b, 0xafa20b9674640eab, 0x09bbcea7d8d9497d,
                               }),
                            Fp({
                                       0x0303cb98b1662daa, 0xd93110aa0a621d5a, 0xbfa9820c5be4a468,
                                       0x0ba3643ecb05a348, 0xdc3534bb1f1c25a6, 0x06c305bb19c0e1c1,
                               }),
                    },
                    Fp2{
                            Fp({
                                       0x46f9cb98b162d858, 0x0be9109cf7aa1d57, 0xc791bc55fece41d2,
                                       0xf84c57704e385ec2, 0xcb49c1d9c010e60f, 0x0acdb8e158bfe3c8,
                               }),
                            Fp({
                                       0x8aefcb98b15f8306, 0x3ea1108fe4f21d54, 0xcf79f69fa117df3b,
                                       0xe4f54aa1d16b1a3c, 0xba5e4ef86105a679, 0x0ed86c0797bee5cf,
                               }),
                    },
                    Fp2{
                            Fp({
                                       0xcee5cb98b15c2db4, 0x71591082d23a1d51, 0xd76230e944a17ca4,
                                       0xd19e3dd3549dd5b6, 0xa972dc1701fa66e3, 0x12e31f2dd6bde7d6,
                               }),
                            Fp({
                                       0xad2acb98b1732d9d, 0x2cfd10dd06961d64, 0x07396b86c6ef24e8,
                                       0xbd76e2fdb1bfc820, 0x6afea7f6de94d0d5, 0x10994b0c5744c040,
                               }),
                    },
            },
    };

    Fp12 c{
            Fp6{
                    Fp2{
                            Fp({
                                       0x47f9cb9871b82d58, 0x5fe911eba3aa1d9d, 0x96bf1b5f4dd81db3,
                                       0x8100d27cc9259f5b, 0xafa20b9674640eab, 0x09bbcea7d8d9497d,
                               }),
                            Fp({
                                       0x0303cb98b1662daa, 0xd93110aa0a621d5a, 0xbfa9820c5be4a468,
                                       0x0ba3643ecb05a348, 0xdc3534bb1f1c25a6, 0x06c305bb19c0e1c1,
                               }),
                    },
                    Fp2{
                            Fp({
                                       0x46f9cb98b162d858, 0x0be9109cf7aa1d57, 0x7791bc55fece41d2,
                                       0xf84c57704e385ec2, 0xcb49c1d9c010e60f, 0x0acdb8e158bfe3c8,
                               }),
                            Fp({
                                       0x8aefcb98b15f8306, 0x3ea1108fe4f21d54, 0xcf79f69fa1b7df3b,
                                       0xe4f54aa1d16b133c, 0xba5e4ef86105a679, 0x0ed86c0797bee5cf,
                               }),
                    },
                    Fp2{
                            Fp({
                                       0xcee5cb98b15c2db4, 0x71591082d23a1d51, 0xd76240e944a17ca4,
                                       0xd19e3dd3549dd5b6, 0xa972dc1701fa66e3, 0x12e31f2dd6bde7d6,
                               }),
                            Fp({
                                       0xad2acb98b1732d9d, 0x2cfd10dd06961d64, 0x07396b86c6ef24e8,
                                       0xbd76e2fdb1bfc820, 0x6afea7f6de94d0d5, 0x10994b0c1744c040,
                               }),
                    },
            },
            Fp6{
                    Fp2{
                            Fp({
                                       0x47f9cb98b1b82d58, 0x5fe911eba3aa1d9d, 0x96bf1b5f4dd81db3,
                                       0x8100d27cc9259f5b, 0xafa20b9674640eab, 0x09bbcea7d8d9497d,
                               }),
                            Fp({
                                       0x0303cb98b1662daa, 0xd93110aa0a621d5a, 0xbfa9820c5be4a468,
                                       0x0ba3643ecb05a348, 0xdc3534bb1f1c25a6, 0x06c305bb19c0e1c1,
                               }),
                    },
                    Fp2{
                            Fp({
                                       0x46f9cb98b162d858, 0x0be9109cf7aa1d57, 0xc791bc55fece41d2,
                                       0xf84c57704e385ec2, 0xcb49c1d3c010e60f, 0x0acdb8e158bfe3c8,
                               }),
                            Fp({
                                       0x8aefcb98b15f8306, 0x3ea1108fe4f21d54, 0xcf79f69fa1b7df3b,
                                       0xe4f54aa1d16b1a3c, 0xba5e4ef86105a679, 0x0ed86c0797bee5cf,
                               }),
                    },
                    Fp2{
                            Fp({
                                       0xcee5cb98b15c2db4, 0x71591082d23a1d51, 0xd76230e944a17ca4,
                                       0xd19e3dd3549dd5b6, 0xa972dc1701fa66e3, 0x12e31f2dd6bde7d6,
                               }),
                            Fp({
                                       0xad2acb98b1732d9d, 0x2cfd10dd06961d64, 0x07396b86c6ef24e8,
                                       0xbd76e2fdb1bfc820, 0x6afea7f6de94d0d5, 0x10994b0c57441040,
                               }),
                    },
            },
    };

    a = a.square().invert().value().square() + c;
    b = b.square().invert().value().square() + a;
    c = c.square().invert().value().square() + b;

    EXPECT_EQ(a.square(), a * a);
    EXPECT_EQ(b.square(), b * b);
    EXPECT_EQ(c.square(), c * c);

    EXPECT_EQ((a + b) * c.square(), (c * c * a) + (c * c * b));

    EXPECT_EQ(a.invert().value() * b.invert().value(), (a * b).invert().value());

    EXPECT_TRUE(a != a.frobenius_map());
    EXPECT_EQ(
            a,
            a.frobenius_map().frobenius_map().frobenius_map().frobenius_map()
                    .frobenius_map().frobenius_map().frobenius_map().frobenius_map()
                    .frobenius_map().frobenius_map().frobenius_map().frobenius_map()
    );
}