#include "scalar/scalar.h"

#include "group/g1_affine.h"
#include "group/g1_projective.h"
#include "scalar/constant.h"
#include "utils/arith.h"
#include "utils/bit.h"
#include "utils/encode.h"
#include "utils/random.h"

namespace bls12_381::scalar {

Scalar::Scalar() : data{0} {}

Scalar::Scalar(uint64_t val) : data{0} {
    *this = Scalar({val, 0, 0, 0}) * constant::R2;
}

Scalar::Scalar(std::array<uint64_t, Scalar::WIDTH> &&data) : data{data} {}

Scalar::Scalar(const std::array<uint64_t, Scalar::WIDTH> &data) : data{data} {}

Scalar Scalar::zero() {
    return Scalar{};
}

Scalar Scalar::one() {
    return constant::R1;
}

Scalar Scalar::random() {
    std::array<uint8_t, Scalar::WIDTH * sizeof(uint64_t) * 2> bytes{};
    for (uint8_t &byte: bytes) byte = bls12_381::util::random::getRandom<uint8_t>();
    return Scalar::from_bytes_wide(bytes);
}

Scalar Scalar::montgomery_reduce(const std::array<uint64_t, Scalar::WIDTH * 2> &ts) {
    uint64_t k;
    uint64_t carry = 0;
    uint64_t carry2 = 0;

    std::array<uint64_t, Scalar::WIDTH * 2> r{};
    for (int i = 0; i < Scalar::WIDTH * 2; ++i)
        r[i] = ts[i];

    for (int i = 0; i < Scalar::WIDTH; ++i) {
        k = r[i] * constant::INV;
        std::tie(std::ignore, carry) = bls12_381::util::arithmetic::mac(r[i], k, constant::MODULUS.data[0], 0);
        for (int j = 1; j < Scalar::WIDTH; ++j)
            std::tie(r[i + j], carry) = bls12_381::util::arithmetic::mac(r[i + j], k, constant::MODULUS.data[j], carry);
        std::tie(r[i + Scalar::WIDTH], carry2) = bls12_381::util::arithmetic::adc(r[i + Scalar::WIDTH], carry2, carry);
    }
    return Scalar({r[4], r[5], r[6], r[7]}).subtract_modulus();
}

Scalar Scalar::from_raw(const std::array<uint64_t, Scalar::WIDTH> &values) {
    return Scalar(values) * constant::R2;
}

Scalar Scalar::from_bytes_wide(const std::array<uint8_t, Scalar::WIDTH * sizeof(uint64_t) * 2> &bytes) {
    std::array<std::array<uint8_t, sizeof(uint64_t)>, Scalar::WIDTH * 2> array{};
    std::array<uint64_t, Scalar::WIDTH * 2> data{};

    for (int i = 0; i < bytes.size(); ++i)
        array[i / sizeof(uint64_t)][i % sizeof(uint64_t)] = bytes[i];
    for (int i = 0; i < data.size(); ++i)
        data[i] = bls12_381::util::bit_operation::le_bytes_to_uint64(array[i]);

    return Scalar::reduce(data);
}

std::optional<Scalar> Scalar::from_bytes(const std::array<uint8_t, Scalar::WIDTH * sizeof(uint64_t)> &bytes) {
    std::array<std::array<uint8_t, sizeof(uint64_t)>, Scalar::WIDTH> array{};
    std::array<uint64_t, Scalar::WIDTH> data{};

    for (int i = 0; i < bytes.size(); ++i)
        array[i / sizeof(uint64_t)][i % sizeof(uint64_t)] = bytes[i];
    for (int i = 0; i < data.size(); ++i)
        data[i] = bls12_381::util::bit_operation::le_bytes_to_uint64(array[i]);

    Scalar temp({data[0], data[1], data[2], data[3]});

    uint64_t borrow = 0;
    for (int i = 0; i < Scalar::WIDTH; ++i)
        std::tie(std::ignore, borrow) = bls12_381::util::arithmetic::sbb(temp.data[i], constant::MODULUS.data[i],
                                                                         borrow);

    uint8_t is_some = static_cast<uint8_t>(borrow) & 1;
    temp *= constant::R2;

    if (is_some) {
        return temp;
    } else {
        return std::nullopt;
    }
}

bool Scalar::is_zero() const {
    return this->data[0] == 0
           && this->data[1] == 0
           && this->data[2] == 0
           && this->data[3] == 0;
}

std::string Scalar::getHex() const {
    std::array<uint8_t, Scalar::WIDTH * sizeof(uint64_t)> bytes = this->to_bytes();
    std::reverse(bytes.begin(), bytes.end());

    std::string res = "0x";
    res += bls12_381::util::encoding::hexStr(bytes);

    return res;
}

std::array<uint8_t, Scalar::WIDTH * sizeof(uint64_t)> Scalar::to_bytes() const {
    std::array<uint64_t, Scalar::WIDTH * 2> contents{0};
    for (int i = 0; i < std::size(contents); ++i)
        if (i < Scalar::WIDTH)
            contents[i] = this->data[i];

    Scalar point = Scalar::montgomery_reduce(contents);

    std::array<uint8_t, sizeof(uint64_t)> temp{};
    std::array<uint8_t, Scalar::WIDTH * sizeof(uint64_t)> bytes{0};

    for (int i = 0; i < Scalar::WIDTH; ++i) {
        temp = bls12_381::util::bit_operation::uint64_to_le_bytes(point.data[i]);
        for (int j = 0; j < sizeof(uint64_t); ++j)
            bytes[i * 8 + j] = temp[j];
    }

    return bytes;
}

Scalar Scalar::doubles() const {
    return *this + *this;
}

Scalar Scalar::square() const {
    uint64_t carry = 0;
    std::array<uint64_t, Scalar::WIDTH * 2> temp{0};

    for (int i = 0; i < Scalar::WIDTH - 1; ++i) {
        carry = 0;
        for (int j = 0; j < Scalar::WIDTH - i - 2; ++j) {
            int32_t anchor = i * 2 + j + 1;
            std::tie(temp[anchor], carry) = bls12_381::util::arithmetic::mac(temp[anchor], this->data[i],
                                                                             this->data[i + j + 1], carry);
        }
        std::tie(temp[i + Scalar::WIDTH - 1], temp[i + Scalar::WIDTH]) =
                bls12_381::util::arithmetic::mac(temp[i + Scalar::WIDTH - 1], this->data[i],
                                                 this->data[Scalar::WIDTH - 1], carry);
    }

    temp[7] = temp[6] >> 63;
    for (int i = 6; i >= 2; --i) temp[i] = (temp[i] << 1) | (temp[i - 1] >> 63);
    temp[1] = temp[1] << 1;

    carry = 0;
    temp[0] = 0;
    for (int i = 0; i < Scalar::WIDTH * 2; ++i) {
        if (i % 2 == 0) {
            std::tie(temp[i], carry) = bls12_381::util::arithmetic::mac(temp[i], this->data[i / 2], this->data[i / 2],
                                                                        carry);
        } else {
            std::tie(temp[i], carry) = bls12_381::util::arithmetic::adc(temp[i], 0, carry);
        }
    }

    return Scalar::montgomery_reduce(temp);
}

Scalar Scalar::subtract_modulus() const {
    return *this - constant::MODULUS;
}

Scalar Scalar::pow_vartime(const std::array<uint64_t, Scalar::WIDTH> &exp) const {
    Scalar res = Scalar::one();
    for (int i = Scalar::WIDTH - 1; i >= 0; --i) {
        for (int j = 63; j >= 0; --j) {
            res = res.square();
            if (((exp[i] >> j) & 0x01) == 0x01) res *= *this;
        }
    }
    return res;
}

std::optional<Scalar> Scalar::sqrt() const {
    std::array<uint64_t, Scalar::WIDTH> exp = {
            0x7fff2dff7fffffff, 0x04d0ec02a9ded201,
            0x94cebea4199cec04, 0x0000000039f6d3a9,
    };
    Scalar w = this->pow_vartime(exp);

    uint32_t v = constant::S;
    Scalar x = *this * w;
    Scalar b = x * w;
    Scalar z = constant::ROOT_OF_UNITY;

    for (int i = constant::S; i >= 1; --i) {
        int32_t k = 1;
        Scalar temp = b.square();
        bool j_less_than_v = true;

        for (int j = 2; j < i; ++j) {
            bool temp_is_one = temp == Scalar::one();
            Scalar squared = (temp_is_one ? z : temp).square();
            temp = temp_is_one ? temp : squared;
            Scalar new_z = temp_is_one ? squared : z;
            j_less_than_v &= (j != v);
            k = temp_is_one ? k : j;
            z = j_less_than_v ? new_z : z;
        }

        Scalar res = x * z;
        x = (b == Scalar::one()) ? x : res;
        z = z.square();
        b *= z;
        v = k;
    }

    if ((x * x) == *this) {
        return x;
    } else {
        return std::nullopt;
    }
}

void multi_square(Scalar &n, int32_t times) {
    for (int i = 0; i < times; ++i)
        n = n.square();
}

std::optional<Scalar> Scalar::invert() const {
    Scalar t0 = this->square();
    Scalar t1 = t0 * *this;
    Scalar t16 = t0.square();
    Scalar t6 = t16.square();
    Scalar t5 = t6 * t0;
    t0 = t6 * t16;
    Scalar t12 = t5 * t16;
    Scalar t2 = t6.square();
    Scalar t7 = t5 * t6;
    Scalar t15 = t0 * t5;
    Scalar t17 = t12.square();
    t1 *= t17;
    Scalar t3 = t7 * t2;
    Scalar t8 = t1 * t17;
    Scalar t4 = t8 * t2;
    Scalar t9 = t8 * t7;
    t7 = t4 * t5;
    Scalar t11 = t4 * t17;
    t5 = t9 * t17;
    Scalar t14 = t7 * t15;
    Scalar t13 = t11 * t12;
    t12 = t11 * t17;
    t15 *= t12;
    t16 *= t15;
    t3 *= t16;
    t17 *= t3;
    t0 *= t17;
    t6 *= t0;
    t2 *= t6;

    multi_square(t0, 8);
    t0 *= t17;
    multi_square(t0, 9);
    t0 *= t16;
    multi_square(t0, 9);
    t0 *= t15;
    multi_square(t0, 9);
    t0 *= t15;
    multi_square(t0, 7);
    t0 *= t14;
    multi_square(t0, 7);
    t0 *= t13;
    multi_square(t0, 10);
    t0 *= t12;
    multi_square(t0, 9);
    t0 *= t11;
    multi_square(t0, 8);
    t0 *= t8;
    multi_square(t0, 8);
    t0 *= *this;
    multi_square(t0, 14);
    t0 *= t9;
    multi_square(t0, 10);
    t0 *= t8;
    multi_square(t0, 15);
    t0 *= t7;
    multi_square(t0, 10);
    t0 *= t6;
    multi_square(t0, 8);
    t0 *= t5;
    multi_square(t0, 16);
    t0 *= t3;
    multi_square(t0, 8);
    t0 *= t2;
    multi_square(t0, 7);
    t0 *= t4;
    multi_square(t0, 9);
    t0 *= t2;
    multi_square(t0, 8);
    t0 *= t3;
    multi_square(t0, 8);
    t0 *= t2;
    multi_square(t0, 8);
    t0 *= t2;
    multi_square(t0, 8);
    t0 *= t2;
    multi_square(t0, 8);
    t0 *= t3;
    multi_square(t0, 8);
    t0 *= t2;
    multi_square(t0, 8);
    t0 *= t2;
    multi_square(t0, 5);
    t0 *= t1;
    multi_square(t0, 5);
    t0 *= t1;

    if (*this == Scalar::zero()) {
        return std::nullopt;
    } else {
        return t0;
    }
}

Scalar Scalar::reduce(const std::array<uint64_t, Scalar::WIDTH * 2> &limbs) {
    Scalar d0({limbs[0], limbs[1], limbs[2], limbs[3]});
    Scalar d1({limbs[4], limbs[5], limbs[6], limbs[7]});
    return d0 * constant::R2 + d1 * constant::R3;
}

Scalar &Scalar::operator=(const Scalar &rhs) {
    if (this == &rhs) return *this;
    for (int i = 0; i < this->data.size(); ++i)
        this->data[i] = rhs.data[i];
    return *this;
}

Scalar &Scalar::operator+=(const Scalar &rhs) {
    uint64_t carry = 0;
    uint64_t d[Scalar::WIDTH];
    for (int i = 0; i < Scalar::WIDTH; ++i)
        std::tie(d[i], carry) = bls12_381::util::arithmetic::adc(this->data[i], rhs.data[i], carry);
    *this = Scalar({d[0], d[1], d[2], d[3]}).subtract_modulus();
    return *this;
}

Scalar &Scalar::operator-=(const Scalar &rhs) {
    uint64_t borrow = 0;
    uint64_t d[Scalar::WIDTH];
    for (int i = 0; i < Scalar::WIDTH; ++i)
        std::tie(d[i], borrow) = bls12_381::util::arithmetic::sbb(this->data[i], rhs.data[i], borrow);

    uint64_t carry = 0;
    for (int i = 0; i < Scalar::WIDTH; ++i)
        std::tie(d[i], carry) = bls12_381::util::arithmetic::adc(d[i], constant::MODULUS.data[i] & borrow, carry);

    *this = Scalar({d[0], d[1], d[2], d[3]});
    return *this;
}

Scalar &Scalar::operator*=(const Scalar &rhs) {
    uint64_t carry = 0;
    std::array<uint64_t, Scalar::WIDTH * 2> temp{0};

    for (int i = 0; i < Scalar::WIDTH; ++i) {
        carry = 0;
        for (int j = 0; j < Scalar::WIDTH - 1; ++j)
            std::tie(temp[i + j], carry) = bls12_381::util::arithmetic::mac(i == 0 ? 0 : temp[i + j], this->data[i],
                                                                            rhs.data[j], carry);
        std::tie(temp[i + Scalar::WIDTH - 1], temp[i + Scalar::WIDTH]) =
                bls12_381::util::arithmetic::mac(i == 0 ? 0 : temp[i + Scalar::WIDTH - 1], this->data[i],
                                                 rhs.data[Scalar::WIDTH - 1], carry);
    }

    *this = Scalar::montgomery_reduce(temp);
    return *this;
}

Scalar Scalar::operator-() const {
    uint64_t borrow = 0;
    uint64_t d[Scalar::WIDTH];

    for (int i = 0; i < Scalar::WIDTH; ++i)
        std::tie(d[i], borrow) = bls12_381::util::arithmetic::sbb(constant::MODULUS.data[i], this->data[i], borrow);

    bool dec = (this->data[0] | this->data[1] | this->data[2] | this->data[3]) == 0;
    uint64_t mask = static_cast<uint64_t>(dec) - 1;

    return Scalar(
            {
                    d[0] & mask, d[1] & mask,
                    d[2] & mask, d[3] & mask,
            }
    );
}

group::G1Projective operator*(const Scalar &a, const group::G1Affine &b) {
    return group::G1Projective(b) *= a;
}

group::G1Projective operator*(const Scalar &a, const group::G1Projective &b) {
    return group::G1Projective(b) *= a;
}
} // namespace bls12_381::scalar