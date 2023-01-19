#include "field/fp.h"

// P = 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
const uint64_t MODULUS_P[6] = {
        0xb9feffffffffaaab,
        0x1eabfffeb153ffff,
        0x6730d2a0f6b0f624,
        0x64774b84f38512bf,
        0x4b1ba7b6434bacd7,
        0x1a0111ea397fe69a,
};

// R = 2 ^ 384 mod p
const Fp R = Fp(
        {
                0x760900000002fffd,
                0xebf4000bc40c0002,
                0x5f48985753c758ba,
                0x77ce585370525745,
                0x5c071a97a256ec6d,
                0x15f65ec3fa80e493
        }
);

// R2 = 2 ^ (384 * 2) mod p
const Fp R2 = Fp(
        {
                0xf4df1f341c341746,
                0x0a76e6a609d104f1,
                0x8de5476c4c95b6d5,
                0x67eb88a9939d83c0,
                0x9a793e85b519952d,
                0x11988fe592cae3aa,
        }
);

// R3 = 2 ^ (384 * 3) mod p
const Fp R3 = Fp(
        {
                0xed48ac6bd94ca1e0,
                0x315f831e03a7adf8,
                0x9a53352a615e29dd,
                0x34c04e5e921e1761,
                0x2512d43565724728,
                0x0aa6346091755d4d,
        }
);

Fp::Fp() : data{} {}

Fp::Fp(uint64_t val) : data{val} {}

Fp::Fp(const std::vector<uint64_t> &vch) : data{} {
    assert(vch.size() == sizeof(this->data));
    std::memcpy(this->data, vch.data(), sizeof(this->data));
}

Fp Fp::zero() {
    return Fp({0, 0, 0, 0, 0, 0});
}

Fp Fp::one() {
    return R;
}

bool Fp::is_zero() {
    return this->data[0] == 0 && this->data[1] == 0
           && this->data[2] == 0 && this->data[3] == 0
           && this->data[4] == 0 && this->data[5] == 0;
}
