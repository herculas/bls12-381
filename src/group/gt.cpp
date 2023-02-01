#include "group/gt.h"
#include "scalar/scalar.h"
#include "pairing/miller_loop.h"

namespace bls12_381::group {

Gt::Gt() : data{field::Fp12::one()} {}

Gt::Gt(const field::Fp12 &point) : data{point} {}

Gt::Gt(field::Fp12 &&point) : data{point} {}

Gt Gt::identity() {
    return Gt{};
}

Gt Gt::generator() {
    return Gt{
            field::Fp12{
                    field::Fp6{
                            field::Fp2{
                                    field::Fp({
                                                      0x1972e433a01f85c5, 0x97d32b76fd772538, 0xc8ce546fc96bcdf9,
                                                      0xcef63e7366d40614, 0xa611342781843780, 0x13f3448a3fc6d825,
                                              }),
                                    field::Fp({
                                                      0xd26331b02e9d6995, 0x9d68a482f7797e7d, 0x9c9b29248d39ea92,
                                                      0xf4801ca2e13107aa, 0xa16c0732bdbcb066, 0x083ca4afba360478,
                                              }),
                            },
                            field::Fp2{
                                    field::Fp({
                                                      0x59e261db0916b641, 0x2716b6f4b23e960d, 0xc8e55b10a0bd9c45,
                                                      0x0bdb0bd99c4deda8, 0x8cf89ebf57fdaac5, 0x12d6b7929e777a5e,
                                              }),
                                    field::Fp({
                                                      0x5fc85188b0e15f35, 0x34a06e3a8f096365, 0xdb3126a6e02ad62c,
                                                      0xfc6f5aa97d9a990b, 0xa12f55f5eb89c210, 0x1723703a926f8889,
                                              }),
                            },
                            field::Fp2{
                                    field::Fp({
                                                      0x93588f2971828778, 0x43f65b8611ab7585, 0x3183aaf5ec279fdf,
                                                      0xfa73d7e18ac99df6, 0x64e176a6a64c99b0, 0x179fa78c58388f1f,
                                              }),
                                    field::Fp({
                                                      0x672a0a11ca2aef12, 0x0d11b9b52aa3f16b, 0xa44412d0699d056e,
                                                      0xc01d0177221a5ba5, 0x66e0cede6c735529, 0x05f5a71e9fddc339,
                                              }),
                            },
                    },
                    field::Fp6{
                            field::Fp2{
                                    field::Fp({
                                                      0xd30a88a1b062c679, 0x5ac56a5d35fc8304, 0xd0c834a6a81f290d,
                                                      0xcd5430c2da3707c7, 0xf0c27ff780500af0, 0x09245da6e2d72eae,
                                              }),
                                    field::Fp({
                                                      0x9f2e0676791b5156, 0xe2d1c8234918fe13, 0x4c9e459f3c561bf4,
                                                      0xa3e85e53b9d3e3c1, 0x820a121e21a70020, 0x15af618341c59acc,
                                              }),
                            },
                            field::Fp2{
                                    field::Fp({
                                                      0x7c95658c24993ab1, 0x73eb38721ca886b9, 0x5256d749477434bc,
                                                      0x8ba41902ea504a8b, 0x04a3d3f80c86ce6d, 0x18a64a87fb686eaa,
                                              }),
                                    field::Fp({
                                                      0xbb83e71bb920cf26, 0x2a5277ac92a73945, 0xfc0ee59f94f046a0,
                                                      0x7158cdf3786058f7, 0x7cc1061b82f945f6, 0x03f847aa9fdbe567,
                                              }),
                            },
                            field::Fp2{
                                    field::Fp({
                                                      0x8078dba56134e657, 0x1cd7ec9a43998a6e, 0xb1aa599a1a993766,
                                                      0xc9a0f62f0842ee44, 0x8e159be3b605dffa, 0x0c86ba0d4af13fc2,
                                              }),
                                    field::Fp({
                                                      0xe80ff2a06a52ffb1, 0x7694ca48721a906c, 0x7583183e03b08514,
                                                      0xf567afdd40cee4e2, 0x9a6d96d2e526a5fc, 0x197e9f49861f2242,
                                              }),
                            },
                    },
            }
    };
}

Gt Gt::random() {
    while (true) {
        field::Fp12 inner = field::Fp12::random();
        if (!inner.is_zero())
            return pairing::MillerLoopResult{inner}.final_exponentiation();
    }
}

bool Gt::is_identity() const {
    return *this == Gt::identity();
}

Gt Gt::doubles() const {
    return Gt{this->data.square()};
}

Gt Gt::operator-() const {
    return Gt{this->data.conjugate()};
}

Gt &Gt::operator=(const Gt &rhs) {
    if (*this == rhs) return *this;
    this->data = rhs.data;
    return *this;
}

Gt &Gt::operator+=(const Gt &rhs) {
    *this = Gt{this->data * rhs.data};
    return *this;
}

Gt &Gt::operator-=(const Gt &rhs) {
    *this = Gt{this->data * (-rhs.data)};
    return *this;
}

Gt &Gt::operator*=(const scalar::Scalar &rhs) {
    Gt acc = Gt::identity();
    auto bytes = rhs.to_bytes();
    for (auto iter = bytes.rbegin(); iter != bytes.rend(); ++iter) {
        for (int i = 7; i >= 0; --i) {
            if (iter == bytes.rbegin() && i == 7) continue;
            uint8_t bit = (*iter >> i) & static_cast<uint8_t>(1);
            acc = acc.doubles();
            if (bit != 0) acc = acc + *this;
        }
    }
    *this = acc;
    return *this;
}

} // namespace bls12_381::group