#ifndef BLS12_381_SCALAR_CONSTANT_H
#define BLS12_381_SCALAR_CONSTANT_H

#include "scalar/scalar.h"

namespace bls12_381::scalar::constant {

const uint32_t MODULUS_LIMBS_32[Scalar::WIDTH * 2] = {
        0x00000001, 0xffffffff, 0xfffe5bfe, 0x53bda402,
        0x09a1d805, 0x3339d808, 0x299d7d48, 0x73eda753,
};

const uint32_t MODULUS_BITS = 255;
const uint32_t S = 32;
const uint64_t INV = 0xfffffffeffffffff;

const Scalar MODULUS(
        {
                0xffffffff00000001, 0x53bda402fffe5bfe,
                0x3339d80809a1d805, 0x73eda753299d7d48,
        }
);

const Scalar GENERATOR(
        {
                0x0000000efffffff1, 0x17e363d300189c0f,
                0xff9c57876f8457b0, 0x351332208fc5a8c4,
        }
);

const Scalar R1(
        {
                0x00000001fffffffe, 0x5884b7fa00034802,
                0x998c4fefecbc4ff5, 0x1824b159acc5056f,
        }
);

const Scalar R2(
        {
                0xc999e990f3f29c6d, 0x2b6cedcb87925c23,
                0x05d314967254398f, 0x0748d9d99f59ff11,
        }
);

const Scalar R3(
        {
                0xc62c1807439b73af, 0x1b3e0d188cf06990,
                0x73d13c71c7b5f418, 0x6e2a5bb9c8db33e9,
        }
);

const Scalar ROOT_OF_UNITY(
        {
                0xb9b58d8c5f0e466a, 0x5b1b4c801819d7ec,
                0x0af53ae352a31e64, 0x5bf3adda19e9b27b,
        }
);

} // namespace bls12_381::scalar::constant

#endif //BLS12_381_SCALAR_CONSTANT_H