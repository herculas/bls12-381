#include "field/fp.h"

#include <cassert>

#ifdef __cpp_lib_ranges
#include <ranges>
#endif

#include "utils/bit.h"

#include "field/constant.h"
#include "utils/arith.h"
#include "utils/encode.h"

namespace bls12_381::field {

using rng::core::RngCore;
using rng::util::bit::from_be_bytes;
using rng::util::bit::to_be_bytes;

using util::arithmetic::adc;
using util::arithmetic::sbb;
using util::arithmetic::mac;
using util::encoding::hex_str;

Fp::Fp() noexcept: data{0} {}

Fp::Fp(const Fp &fp) noexcept = default;

Fp::Fp(uint64_t val) noexcept: data{val} {}

Fp::Fp(const std::array<uint64_t, Fp::WIDTH> &data) noexcept: data{data} {}

Fp::Fp(Fp &&fp) noexcept = default;

Fp::Fp(std::array<uint64_t, Fp::WIDTH> &&data) noexcept: data{data} {}

Fp::~Fp() noexcept = default;

Fp Fp::zero() noexcept {
    return Fp{};
}

Fp Fp::one() noexcept {
    return constant::R1;
}

Fp Fp::random(RngCore &rng) {
    std::array<uint8_t, Fp::BYTE_SIZE * 2> randoms{};
    rng.fill_bytes(randoms);

    std::array<uint64_t, Fp::WIDTH * 2> bytes{};
    for (int i = 0; i < bytes.size(); ++i) {
        bytes[i] = from_be_bytes<uint64_t>(
                {
                        randoms[8 * i + 0], randoms[8 * i + 1], randoms[8 * i + 2], randoms[8 * i + 3],
                        randoms[8 * i + 4], randoms[8 * i + 5], randoms[8 * i + 6], randoms[8 * i + 7],
                }
        );
    }
    return Fp::reduce(bytes);
}

Fp Fp::montgomery_reduce(const std::array<uint64_t, Fp::WIDTH * 2> &ts) {
    uint64_t carry = 0;

    std::array<uint64_t, Fp::WIDTH * 2> r{};
    for (int i = 0; i < Fp::WIDTH * 2; ++i)
        r[i] = i < Fp::WIDTH ? ts[i] : 0;

    for (int i = 0; i < Fp::WIDTH; ++i) {
        const uint64_t k = r[i] * constant::INV;
        carry = 0;
        r[0] = mac(r[i], k, constant::MODULUS[0], carry);
        for (int j = 1; j < Fp::WIDTH; ++j)
            r[i + j] = mac(r[i + j], k, constant::MODULUS[j], carry);
        r[i + Fp::WIDTH] = adc(ts[i + Fp::WIDTH], r[i + Fp::WIDTH], carry);
        r[(i + Fp::WIDTH + 1) % (Fp::WIDTH * 2)] = carry;
    }

    return Fp({r[6], r[7], r[8], r[9], r[10], r[11]}).subtract_modulus();
}

Fp Fp::sum_of_products(const std::vector<Fp> &a, const std::vector<Fp> &b) {
    assert(a.size() == b.size());
    const auto len = static_cast<int32_t>(a.size());

    uint64_t u[Fp::WIDTH] = {0};
    uint64_t carry = 0;

    for (int j = 0; j < Fp::WIDTH; ++j) {
        uint64_t t[Fp::WIDTH + 1] = {0};
        std::copy(u, u + Fp::WIDTH, t); // NOLINT
        for (int i = 0; i < len; ++i) {
            carry = 0;
            for (int k = 0; k < Fp::WIDTH; ++k)
                t[k] = mac(t[k], a[i].data[j], b[i].data[k], carry);
            t[Fp::WIDTH] = adc(t[Fp::WIDTH], 0, carry);
        }

        uint64_t r[Fp::WIDTH + 1] = {0};
        const uint64_t k = t[0] * constant::INV;
        carry = 0;
        mac(t[0], k, constant::MODULUS[0], carry);

        for (int i = 1; i < Fp::WIDTH; ++i)
            r[i] = mac(t[i], k, constant::MODULUS[i], carry);
        r[Fp::WIDTH] = adc(t[Fp::WIDTH], 0, carry);

        for (int i = 0; i < Fp::WIDTH; ++i)
            u[i] = r[i + 1];
    }

    return Fp({u[0], u[1], u[2], u[3], u[4], u[5]}).subtract_modulus();
}

/// Tries to convert a big-endian byte representation of a scalar into an `Fp`.
std::optional<Fp> Fp::from_bytes(const std::array<uint8_t, Fp::BYTE_SIZE> &bytes) {
    std::array<std::array<uint8_t, sizeof(uint64_t)>, Fp::WIDTH> array{};
    std::array<uint64_t, Fp::WIDTH> data{};

    for (int i = 0; i < bytes.size(); ++i)
        array[i / sizeof(uint64_t)][i % sizeof(uint64_t)] = bytes[i];
    for (int i = 0; i < data.size(); ++i)
        data[i] = from_be_bytes<uint64_t>(array[Fp::WIDTH - i - 1]);

    Fp temp({data[0], data[1], data[2], data[3], data[4], data[5]});

    // try to subtract the modulus
    uint64_t borrow = 0;
    for (int i = 0; i < Fp::WIDTH; ++i)
        sbb(temp.data[i], constant::MODULUS[i], borrow);

    // if the element is smaller than the modulus, the subtraction would underflow, generating a `borrow` of
    // 0xfff...fff. Otherwise, it would be 0.
    if (borrow == 0) {
        return std::nullopt;
    } else {
        // convert to Montgomery form by computing (a.R ^ 0 * R ^ 2) / R = a.R
        return temp *= constant::R2;
    }
}

bool Fp::is_zero() const {
    return this->data[0] == 0 && this->data[1] == 0
           && this->data[2] == 0 && this->data[3] == 0
           && this->data[4] == 0 && this->data[5] == 0;
}

/// Computes if this element is strictly lexicographically larger than its negation. This can be determined by
/// checking if the element is larger than (p - 1) / 2. If we subtract by (p - 1) / 2 + 1 and there is no
/// underflow, then the element must be larger than (p - 1) / 2.
bool Fp::lexicographically_largest() const {
    std::array<uint64_t, Fp::WIDTH * 2> contents{0};
    for (int i = 0; i < std::size(contents); ++i)
        if (i < WIDTH)
            contents[i] = this->data[i];

    const Fp temp = Fp::montgomery_reduce(contents);

    const uint64_t subs[WIDTH] = {
            0xdcff7fffffffd556, 0x0f55ffff58a9ffff, 0xb39869507b587b12,
            0xb23ba5c279c2895f, 0x258dd3db21a5d66b, 0x0d0088f51cbff34d
    };

    uint64_t borrow = 0;
    for (int i = 0; i < WIDTH; ++i)
        sbb(temp.data[i], subs[i], borrow);

    // If the element is smaller, the subtraction would underflow, producing a borrow of 0xfff...fff,
    // otherwise, it will be 0x000...000.
    return !(static_cast<uint8_t>(borrow) & 1);
}

const std::array<uint64_t, Fp::WIDTH> &Fp::get_data() const {
    return this->data;
}

/// Converts an element of `Fp` into a big-endian byte array.
std::array<uint8_t, Fp::BYTE_SIZE> Fp::to_bytes() const {
    std::array<uint64_t, Fp::WIDTH * 2> contents{0};
    for (int i = 0; i < Fp::WIDTH; ++i)
        contents[i] = this->data[i];

    // turn this into canonical form by computing (a.R) / R = a.
    const Fp point = Fp::montgomery_reduce(contents);

    std::array<uint8_t, sizeof(uint64_t)> temp{};
    std::array<uint8_t, Fp::BYTE_SIZE> bytes{0};

    for (int i = 0; i < Fp::WIDTH; ++i) {
        temp = to_be_bytes<uint64_t>(point.data[Fp::WIDTH - 1 - i]);
        for (int j = 0; j < sizeof(uint64_t); ++j)
            bytes[i * 8 + j] = temp[j];
    }
    return bytes;
}

std::string Fp::to_hex_str() const {
    std::array<uint8_t, Fp::BYTE_SIZE> bytes = this->to_bytes();

    std::string res = "0x";
    res += hex_str(bytes);
    return res;
}

Fp Fp::square() const {
    uint64_t carry = 0;
    std::array<uint64_t, Fp::WIDTH * 2> temp{0};

    for (int i = 0; i < Fp::WIDTH - 1; ++i) {
        carry = 0;
        for (int j = 0; j < Fp::WIDTH - i - 2; ++j) {
            const int32_t anchor = i * 2 + j + 1;
            temp[anchor] = mac(temp[anchor], this->data[i], this->data[i + j + 1], carry);
        }
        temp[i + Fp::WIDTH - 1] = mac(temp[i + Fp::WIDTH - 1], this->data[i],
                                      this->data[Fp::WIDTH - 1], carry);
        temp[i + Fp::WIDTH] = carry;
    }

    temp[11] = temp[10] >> 63;
    for (int i = 10; i >= 2; --i) temp[i] = (temp[i] << 1) | (temp[i - 1] >> 63);
    temp[1] = temp[1] << 1;

    carry = 0;
    temp[0] = 0;
    for (int i = 0; i < Fp::WIDTH * 2; ++i) {
        if (i % 2 == 0) {
            temp[i] = mac(temp[i], this->data[i / 2], this->data[i / 2], carry);
        } else {
            temp[i] = adc(temp[i], 0, carry);
        }
    }
    return Fp::montgomery_reduce(temp);
}

Fp Fp::subtract_modulus() const {
    uint64_t borrow = 0;
    uint64_t r[Fp::WIDTH] = {0};
    uint64_t d[Fp::WIDTH] = {0};

    for (int i = 0; i < Fp::WIDTH; ++i)
        r[i] = sbb(this->data[i], constant::MODULUS[i], borrow);

    // If underflow occurs on the final limb, borrow = 0xfff...fff, otherwise borrow = 0x000...000.
    // Thus, we use it as a mask.
    const uint64_t borrow_flip = borrow ^ 0xffffffffffffffff;

    for (int i = 0; i < Fp::WIDTH; ++i) {
        d[i] = (this->data[i] & borrow) | (r[i] & borrow_flip);
    }

    return Fp({d[0], d[1], d[2], d[3], d[4], d[5]});
}

Fp Fp::pow(const std::array<uint64_t, Fp::WIDTH> &exp) const {
    Fp res = Fp::one();
#ifdef __cpp_lib_ranges
    for (const uint64_t &e: exp | std::views::reverse) {
        for (const int32_t &i: std::views::iota(0,64) | std::views::reverse) {
            res = res.square();
            if (((e >> i) & 0x01) == 0x01) res *= *this;
        }
    }
#else
    for (int i = Fp::WIDTH - 1; i >= 0; i--) {
        for (int32_t j = 63; j >= 0; j--) {
            res = res.square();
            if (((exp[i] >> j) & 0x01) == 0x01) res *= *this;
        }
    }
#endif
    return res;
}

std::optional<Fp> Fp::sqrt() const {
    const std::array<uint64_t, Fp::WIDTH> exp = {
            0xee7fbfffffffeaab, 0x07aaffffac54ffff, 0xd9cc34a83dac3d89,
            0xd91dd2e13ce144af, 0x92c6e9ed90d2eb35, 0x0680447a8e5ff9a6,
    };
    const Fp sqrt = this->pow(exp);

    if (sqrt.square() == *this) {
        return sqrt;
    } else {
        return std::nullopt;
    }
}

/// Computes the multiplicative inverse of this field element.
/// Returns `null` when this element is zero.
std::optional<Fp> Fp::invert() const {
    // modulus - 2
    const std::array<uint64_t, Fp::WIDTH> exp = {
            0xb9feffffffffaaa9, 0x1eabfffeb153ffff, 0x6730d2a0f6b0f624,
            0x64774b84f38512bf, 0x4b1ba7b6434bacd7, 0x1a0111ea397fe69a,
    };
    const Fp result = this->pow(exp);

    if (this->is_zero()) {
        return std::nullopt;
    } else {
        return result;
    }
}

/// Reduces a big-endian 64-bit limb representation of a 768-bit number.
Fp Fp::reduce(const std::array<uint64_t, Fp::WIDTH * 2> &limbs) {
    const Fp d1({limbs[11], limbs[10], limbs[9], limbs[8], limbs[7], limbs[6]});
    const Fp d0({limbs[5], limbs[4], limbs[3], limbs[2], limbs[1], limbs[0]});
    return d0 * constant::R2 + d1 * constant::R3;
}

Fp &Fp::operator=(const Fp &rhs) = default;

Fp &Fp::operator=(Fp &&rhs) noexcept = default;

Fp &Fp::operator+=(const Fp &rhs) {
    uint64_t carry = 0;
    uint64_t d[Fp::WIDTH];
    for (int i = 0; i < Fp::WIDTH; ++i)
        d[i] = adc(this->data[i], rhs.data[i], carry);
    *this = Fp({d[0], d[1], d[2], d[3], d[4], d[5]}).subtract_modulus();
    return *this;
}

Fp &Fp::operator-=(const Fp &rhs) {
    *this += -rhs;
    return *this;
}

Fp &Fp::operator*=(const Fp &rhs) {
    uint64_t carry = 0;
    std::array<uint64_t, Fp::WIDTH * 2> temp{0};

    for (int i = 0; i < Fp::WIDTH; ++i) {
        carry = 0;
        for (int j = 0; j < Fp::WIDTH - 1; ++j)
            temp[i + j] = mac(i == 0 ? 0 : temp[i + j], this->data[i], rhs.data[j], carry);
        temp[i + Fp::WIDTH - 1] = mac(i == 0 ? 0 : temp[i + Fp::WIDTH - 1], this->data[i], rhs.data[Fp::WIDTH - 1],
                                      carry);
        temp[i + Fp::WIDTH] = carry;
    }

    *this = Fp::montgomery_reduce(temp);
    return *this;
}

Fp Fp::operator-() const {
    uint64_t borrow = 0;
    uint64_t d[Fp::WIDTH];

    for (int i = 0; i < Fp::WIDTH; ++i)
        d[i] = sbb(constant::MODULUS[i], this->data[i], borrow);

    const bool dec =
            (this->data[0] | this->data[1] | this->data[2] | this->data[3] | this->data[4] | this->data[5]) == 0;
    const uint64_t mask = static_cast<uint64_t>(dec) - 1;

    return Fp(
            {
                    d[0] & mask, d[1] & mask, d[2] & mask,
                    d[3] & mask, d[4] & mask, d[5] & mask
            }
    );
}

Fp Fp::operator+(const Fp &rhs) const {
    uint64_t carry = 0;
    uint64_t d[Fp::WIDTH];
    for (int i = 0; i < Fp::WIDTH; ++i)
        d[i] = adc(this->data[i], rhs.data[i], carry);
    return Fp({d[0], d[1], d[2], d[3], d[4], d[5]}).subtract_modulus();
}

Fp Fp::operator-(const Fp &rhs) const {
    return *this + -rhs;
}

Fp Fp::operator*(const Fp &rhs) const {
    uint64_t carry = 0;
    std::array<uint64_t, Fp::WIDTH * 2> temp{0};

    for (int i = 0; i < Fp::WIDTH; ++i) {
        carry = 0;
        for (int j = 0; j < Fp::WIDTH - 1; ++j)
            temp[i + j] = mac(i == 0 ? 0 : temp[i + j], this->data[i], rhs.data[j], carry);
        temp[i + Fp::WIDTH - 1] = mac(i == 0 ? 0 : temp[i + Fp::WIDTH - 1], this->data[i], rhs.data[Fp::WIDTH - 1],
                                      carry);
        temp[i + Fp::WIDTH] = carry;
    }
    return Fp::montgomery_reduce(temp);
}

} // namespace bls12_381::field