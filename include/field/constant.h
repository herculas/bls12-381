#ifndef BLS12_381_FIELD_CONSTANT_H
#define BLS12_381_FIELD_CONSTANT_H

#include "field/fp.h"

/// The modulus of the finite field.
const uint64_t MODULUS[Fp::WIDTH] = {
        0xb9feffffffffaaab, 0x1eabfffeb153ffff, 0x6730d2a0f6b0f624,
        0x64774b84f38512bf, 0x4b1ba7b6434bacd7, 0x1a0111ea397fe69a,
};

/// INV = -(p ^ {-1} mod 2 ^ 64) mod 2 ^ 64.
const uint64_t INV = 0x89f3fffcfffcfffd;

/// R1 = 2 ^ 384 mod p.
const Fp R1(
        {
                0x760900000002fffd, 0xebf4000bc40c0002, 0x5f48985753c758ba,
                0x77ce585370525745, 0x5c071a97a256ec6d, 0x15f65ec3fa80e493,
        }
);

/// R2 = 2 ^ (384 * 2) mod p.
const Fp R2(
        {
                0xf4df1f341c341746, 0x0a76e6a609d104f1, 0x8de5476c4c95b6d5,
                0x67eb88a9939d83c0, 0x9a793e85b519952d, 0x11988fe592cae3aa,
        }
);

/// R3 = 2 ^ (384 * 3) mod p.
const Fp R3(
        {
                0xed48ac6bd94ca1e0, 0x315f831e03a7adf8, 0x9a53352a615e29dd,
                0x34c04e5e921e1761, 0x2512d43565724728, 0x0aa6346091755d4d,
        }
);

const Fp B(
        {
                0xaa270000000cfff3, 0x53cc0032fc34000a, 0x478fe97a6b0a807f,
                0xb1d37ebee6ba24d7, 0x8ec9733bbf78ab2f, 0x09d645513d83de7e,
        }
);

/// A non-trivial third root of unity in Fp.
const Fp BETA(
        {
                0x30f1361b798a64e8, 0xf3b8ddab7ece5a2a, 0x16a8ca3ac61577f7,
                0xc26a2ff874fd029b, 0x3636b76660701c6e, 0x051ba4ab241b6160,
        }
);

#endif //BLS12_381_FIELD_CONSTANT_H