#include "field/fp.h"

#include <vector>
#include "utils/random.h"
#include "utils/arith.h"
#include "utils/encode.h"

/// The modulus of the finite field.
const uint64_t MODULUS[6] = {
        0xb9feffffffffaaab, 0x1eabfffeb153ffff, 0x6730d2a0f6b0f624,
        0x64774b84f38512bf, 0x4b1ba7b6434bacd7, 0x1a0111ea397fe69a
};

/// INV = -(p ^ {-1} mod 2 ^ 64) mod 2 ^ 64.
const uint64_t INV = 0x89f3fffcfffcfffd;

/// R1 = 2 ^ 384 mod p.
const Fp R1 = Fp({
                         0x760900000002fffd, 0xebf4000bc40c0002, 0x5f48985753c758ba,
                         0x77ce585370525745, 0x5c071a97a256ec6d, 0x15f65ec3fa80e493
                 });

/// R2 = 2 ^ (384 * 2) mod p.
const Fp R2 = Fp({
                         0xf4df1f341c341746, 0x0a76e6a609d104f1, 0x8de5476c4c95b6d5,
                         0x67eb88a9939d83c0, 0x9a793e85b519952d, 0x11988fe592cae3aa,
                 });

/// R3 = 2 ^ (384 * 3) mod p.
const Fp R3 = Fp({
                         0xed48ac6bd94ca1e0, 0x315f831e03a7adf8, 0x9a53352a615e29dd,
                         0x34c04e5e921e1761, 0x2512d43565724728, 0x0aa6346091755d4d,
                 });

Fp::Fp() : data{} {}

Fp::Fp(uint64_t val) : data{val} {}

Fp::Fp(const std::vector<uint64_t> &vch) : data{} {
    assert(vch.size() * 8 == sizeof(this->data));
    std::memcpy(this->data, vch.data(), sizeof(this->data));
}

Fp Fp::zero() {
    return {};
}

Fp Fp::one() {
    return R1;
}

Fp Fp::random() {
    uint64_t randoms[WIDTH * 2];
    for (uint64_t &random: randoms) random = getRandom();
    return reduce(randoms);
}

/// Tries to convert a big-endian byte representation of a scalar into an `Fp`.
std::optional<Fp> Fp::from_bytes(const std::span<uint8_t> bytes) {
    assert(bytes.size() == Fp::WIDTH * sizeof(uint64_t));

    uint8_t array[WIDTH][sizeof(uint64_t)];
    for (int i = 0; i < bytes.size(); ++i) {
        array[i / sizeof(uint64_t)][i % sizeof(uint64_t)] = bytes[i];
    }

    Fp temp({
                    be_bytes_to_uint64(array[5]),
                    be_bytes_to_uint64(array[4]),
                    be_bytes_to_uint64(array[3]),
                    be_bytes_to_uint64(array[2]),
                    be_bytes_to_uint64(array[1]),
                    be_bytes_to_uint64(array[0])
            });

    // try to subtract the modulus
    uint64_t borrow = 0;
    for (int i = 0; i < WIDTH; ++i) {
        std::tie(std::ignore, borrow) = sbb(temp.data[i], MODULUS[i], borrow);
    }

    // if the element is smaller than the modulus, the subtraction would underflow, generating a `borrow` of
    // 0xfff...fff. Otherwise, it would be 0.
    if (borrow == 0) {
        return std::nullopt;
    } else {
        // convert to Montgomery form by computing (a.R ^ 0 * R ^ 2) / R = a.R
        return temp *= R2;
    }
}

/// Converts an element of `Fp` into a big-endian byte array.
uint8_t *Fp::to_bytes(std::span<uint8_t> bytes) const {
    assert(bytes.size() == Fp::WIDTH * sizeof(uint64_t));

    uint64_t contents[WIDTH * 2];
    for (int i = 0; i < std::size(contents); ++i) {
        if (i < WIDTH)
            contents[i] = this->data[i];
        else
            contents[i] = 0;
    }

    // turn this into canonical form by computing (a.R) / R = a.
    Fp point = Fp::montgomery_reduce(contents);

    uint8_t temp[sizeof(uint64_t)];
    for (int i = 0; i < WIDTH; ++i) {
        uint64_to_be_bytes(point.data[WIDTH - 1 - i], temp);
        for (int j = 0; j < sizeof(uint64_t); ++j) {
            bytes[i * 8 + j] = temp[j];
        }
    }
    return bytes.data();
}

std::string Fp::getHex() const {
    std::string res;

    uint8_t bytes[48];
    std::ignore = this->to_bytes(bytes);

    res += "0x";
    res += hexStr(bytes);

    return res;
}

Fp Fp::montgomery_reduce(std::span<uint64_t> ts) {
    assert(ts.size() == Fp::WIDTH * 2);

    uint64_t k;
    uint64_t carry = 0;
    uint64_t r0 = 0, r1 = 0, r2 = 0, r3 = 0, r4 = 0, r5 = 0, r6 = 0, r7 = 0, r8 = 0, r9 = 0, r10 = 0, r11 = 0;

    k = ts[0] * INV;
    std::tie(r0, carry) = mac(ts[0], k, MODULUS[0], 0);
    std::tie(r1, carry) = mac(ts[1], k, MODULUS[1], carry);
    std::tie(r2, carry) = mac(ts[2], k, MODULUS[2], carry);
    std::tie(r3, carry) = mac(ts[3], k, MODULUS[3], carry);
    std::tie(r4, carry) = mac(ts[4], k, MODULUS[4], carry);
    std::tie(r5, carry) = mac(ts[5], k, MODULUS[5], carry);
    std::tie(r6, r7) = adc(ts[6], 0, carry);

    k = r1 * INV;
    std::tie(r0, carry) = mac(r1, k, MODULUS[0], 0);
    std::tie(r2, carry) = mac(r2, k, MODULUS[1], carry);
    std::tie(r3, carry) = mac(r3, k, MODULUS[2], carry);
    std::tie(r4, carry) = mac(r4, k, MODULUS[3], carry);
    std::tie(r5, carry) = mac(r5, k, MODULUS[4], carry);
    std::tie(r6, carry) = mac(r6, k, MODULUS[5], carry);
    std::tie(r7, r8) = adc(ts[7], r7, carry);

    k = r2 * INV;
    std::tie(r0, carry) = mac(r2, k, MODULUS[0], 0);
    std::tie(r3, carry) = mac(r3, k, MODULUS[1], carry);
    std::tie(r4, carry) = mac(r4, k, MODULUS[2], carry);
    std::tie(r5, carry) = mac(r5, k, MODULUS[3], carry);
    std::tie(r6, carry) = mac(r6, k, MODULUS[4], carry);
    std::tie(r7, carry) = mac(r7, k, MODULUS[5], carry);
    std::tie(r8, r9) = adc(ts[8], r8, carry);

    k = r3 * INV;
    std::tie(r0, carry) = mac(r3, k, MODULUS[0], 0);
    std::tie(r4, carry) = mac(r4, k, MODULUS[1], carry);
    std::tie(r5, carry) = mac(r5, k, MODULUS[2], carry);
    std::tie(r6, carry) = mac(r6, k, MODULUS[3], carry);
    std::tie(r7, carry) = mac(r7, k, MODULUS[4], carry);
    std::tie(r8, carry) = mac(r8, k, MODULUS[5], carry);
    std::tie(r9, r10) = adc(ts[9], r9, carry);

    k = r4 * INV;
    std::tie(r0, carry) = mac(r4, k, MODULUS[0], 0);
    std::tie(r5, carry) = mac(r5, k, MODULUS[1], carry);
    std::tie(r6, carry) = mac(r6, k, MODULUS[2], carry);
    std::tie(r7, carry) = mac(r7, k, MODULUS[3], carry);
    std::tie(r8, carry) = mac(r8, k, MODULUS[4], carry);
    std::tie(r9, carry) = mac(r9, k, MODULUS[5], carry);
    std::tie(r10, r11) = adc(ts[10], r10, carry);

    k = r5 * INV;
    std::tie(r0, carry) = mac(r5, k, MODULUS[0], 0);
    std::tie(r6, carry) = mac(r6, k, MODULUS[1], carry);
    std::tie(r7, carry) = mac(r7, k, MODULUS[2], carry);
    std::tie(r8, carry) = mac(r8, k, MODULUS[3], carry);
    std::tie(r9, carry) = mac(r9, k, MODULUS[4], carry);
    std::tie(r10, carry) = mac(r10, k, MODULUS[5], carry);
    std::tie(r11, r0) = adc(ts[11], r11, carry);

    return Fp({r6, r7, r8, r9, r10, r11}).subtract_modulus();
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
    uint64_t contents[WIDTH * 2];
    for (int i = 0; i < std::size(contents); ++i) {
        if (i < WIDTH)
            contents[i] = this->data[i];
        else
            contents[i] = 0;
    }
    Fp temp = Fp::montgomery_reduce(contents);

    uint64_t subs[WIDTH] = {
            0xdcff7fffffffd556, 0x0f55ffff58a9ffff, 0xb39869507b587b12,
            0xb23ba5c279c2895f, 0x258dd3db21a5d66b, 0x0d0088f51cbff34d
    };

    uint64_t borrow = 0;
    for (int i = 0; i < WIDTH; ++i) {
        std::tie(std::ignore, borrow) = sbb(temp.data[i], subs[i], borrow);
    }

    // If the element is smaller, the subtraction would underflow, producing a borrow of 0xfff...fff,
    // otherwise, it will be 0x000...000.

    return !(static_cast<uint8_t>(borrow) & 1);
}

Fp Fp::pow_vartime(std::span<uint64_t> exp) const {
    assert(exp.size() == Fp::WIDTH);

    Fp res = Fp::one();
    for (int i = 5; i >= 0; i--) {
        for (int32_t j = 63; j >= 0; j--) {
            res = res.square();

            if (((exp[i] >> j) & 0x01) == 0x01) {
                res *= *this;
            }
        }
    }
    return res;
}

std::optional<Fp> Fp::sqrt() const {
    uint64_t exp[6] = {
            0xee7fbfffffffeaab, 0x07aaffffac54ffff, 0xd9cc34a83dac3d89,
            0xd91dd2e13ce144af, 0x92c6e9ed90d2eb35, 0x0680447a8e5ff9a6,
    };

    Fp sqrt = this->pow_vartime(exp);

    Fp temp = sqrt.square();
    if (temp == *this) {
        return sqrt;
    } else {
        return std::nullopt;
    }
}

Fp Fp::square() const {
    uint64_t carry = 0;
    uint64_t t0 = 0, t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5 = 0, t6 = 0, t7 = 0, t8 = 0, t9 = 0, t10 = 0, t11 = 0;

    std::tie(t1, carry) = mac(0, this->data[0], this->data[1], 0);
    std::tie(t2, carry) = mac(0, this->data[0], this->data[2], carry);
    std::tie(t3, carry) = mac(0, this->data[0], this->data[3], carry);
    std::tie(t4, carry) = mac(0, this->data[0], this->data[4], carry);
    std::tie(t5, t6) = mac(0, this->data[0], this->data[5], carry);

    std::tie(t3, carry) = mac(t3, this->data[1], this->data[2], 0);
    std::tie(t4, carry) = mac(t4, this->data[1], this->data[3], carry);
    std::tie(t5, carry) = mac(t5, this->data[1], this->data[4], carry);
    std::tie(t6, t7) = mac(t6, this->data[1], this->data[5], carry);

    std::tie(t5, carry) = mac(t5, this->data[2], this->data[3], 0);
    std::tie(t6, carry) = mac(t6, this->data[2], this->data[4], carry);
    std::tie(t7, t8) = mac(t7, this->data[2], this->data[5], carry);

    std::tie(t7, carry) = mac(t7, this->data[3], this->data[4], 0);
    std::tie(t8, t9) = mac(t8, this->data[3], this->data[5], carry);

    std::tie(t9, t10) = mac(t9, this->data[4], this->data[5], 0);

    t11 = t10 >> 63;
    t10 = (t10 << 1) | (t9 >> 63);
    t9 = (t9 << 1) | (t8 >> 63);
    t8 = (t8 << 1) | (t7 >> 63);
    t7 = (t7 << 1) | (t6 >> 63);
    t6 = (t6 << 1) | (t5 >> 63);
    t5 = (t5 << 1) | (t4 >> 63);
    t4 = (t4 << 1) | (t3 >> 63);
    t3 = (t3 << 1) | (t2 >> 63);
    t2 = (t2 << 1) | (t1 >> 63);
    t1 = t1 << 1;

    std::tie(t0, carry) = mac(0, this->data[0], this->data[0], 0);
    std::tie(t1, carry) = adc(t1, 0, carry);
    std::tie(t2, carry) = mac(t2, this->data[1], this->data[1], carry);
    std::tie(t3, carry) = adc(t3, 0, carry);
    std::tie(t4, carry) = mac(t4, this->data[2], this->data[2], carry);
    std::tie(t5, carry) = adc(t5, 0, carry);
    std::tie(t6, carry) = mac(t6, this->data[3], this->data[3], carry);
    std::tie(t7, carry) = adc(t7, 0, carry);
    std::tie(t8, carry) = mac(t8, this->data[4], this->data[4], carry);
    std::tie(t9, carry) = adc(t9, 0, carry);
    std::tie(t10, carry) = mac(t10, this->data[5], this->data[5], carry);
    std::tie(t11, carry) = adc(t11, 0, carry);

    uint64_t res[12] = {
            t0, t1, t2, t3,
            t4, t5, t6, t7,
            t8, t9, t10, t11
    };

    return Fp::montgomery_reduce(res);
}

/// Computes the multiplicative inverse of this field element.
/// Returns `null` when this element is zero.
std::optional<Fp> Fp::invert() const {
    // modulus - 2
    uint64_t exp[6] = {
            0xb9feffffffffaaa9, 0x1eabfffeb153ffff, 0x6730d2a0f6b0f624,
            0x64774b84f38512bf, 0x4b1ba7b6434bacd7, 0x1a0111ea397fe69a,
    };
    Fp result = this->pow_vartime(exp);
    if (this->is_zero()) {
        return std::nullopt;
    } else {
        return result;
    }
}

Fp Fp::subtract_modulus() const {
    uint64_t borrow = 0;
    uint64_t r0 = 0, r1 = 0, r2 = 0, r3 = 0, r4 = 0, r5 = 0;

    std::tie(r0, borrow) = sbb(this->data[0], MODULUS[0], 0);
    std::tie(r1, borrow) = sbb(this->data[1], MODULUS[1], borrow);
    std::tie(r2, borrow) = sbb(this->data[2], MODULUS[2], borrow);
    std::tie(r3, borrow) = sbb(this->data[3], MODULUS[3], borrow);
    std::tie(r4, borrow) = sbb(this->data[4], MODULUS[4], borrow);
    std::tie(r5, borrow) = sbb(this->data[5], MODULUS[5], borrow);

    // If underflow occurs on the final limb, borrow = 0xfff...fff, otherwise borrow = 0x000...000.
    // Thus, we use it as a mask.
    uint64_t borrow_flip = borrow ^ 0xffffffffffffffff;

    uint64_t d0 = (this->data[0] & borrow) | (r0 & borrow_flip);
    uint64_t d1 = (this->data[1] & borrow) | (r1 & borrow_flip);
    uint64_t d2 = (this->data[2] & borrow) | (r2 & borrow_flip);
    uint64_t d3 = (this->data[3] & borrow) | (r3 & borrow_flip);
    uint64_t d4 = (this->data[4] & borrow) | (r4 & borrow_flip);
    uint64_t d5 = (this->data[5] & borrow) | (r5 & borrow_flip);

    return Fp({d0, d1, d2, d3, d4, d5});
}

template<uint32_t N>
Fp Fp::sum_of_products(std::span<Fp> a, std::span<Fp> b) {
    assert(a.size() == N);
    assert(b.size() == N);

    uint64_t u0 = 0, u1 = 0, u2 = 0, u3 = 0, u4 = 0, u5 = 0;

    for (int j = 0; j < WIDTH; ++j) {
        uint64_t t0 = u0, t1 = u1, t2 = u2, t3 = u3, t4 = u4, t5 = u5, t6 = 0;

        for (int i = 0; i < N; ++i) {
            uint64_t carry = 0;
            std::tie(t0, carry) = mac(t0, a[i].data[j], b[i].data[0], carry);
            std::tie(t1, carry) = mac(t1, a[i].data[j], b[i].data[1], carry);
            std::tie(t2, carry) = mac(t2, a[i].data[j], b[i].data[2], carry);
            std::tie(t3, carry) = mac(t3, a[i].data[j], b[i].data[3], carry);
            std::tie(t4, carry) = mac(t4, a[i].data[j], b[i].data[4], carry);
            std::tie(t5, carry) = mac(t5, a[i].data[j], b[i].data[5], carry);
            std::tie(t6, std::ignore) = adc(t6, 0, carry);
        }

        uint64_t r1, r2, r3, r4, r5, r6;
        uint64_t carry = 0;

        uint64_t k = t0 * INV;
        std::tie(std::ignore, carry) = mac(t0, k, MODULUS[0], carry);
        std::tie(r1, carry) = mac(t1, k, MODULUS[1], carry);
        std::tie(r2, carry) = mac(t2, k, MODULUS[2], carry);
        std::tie(r3, carry) = mac(t3, k, MODULUS[3], carry);
        std::tie(r4, carry) = mac(t4, k, MODULUS[4], carry);
        std::tie(r5, carry) = mac(t5, k, MODULUS[5], carry);
        std::tie(r6, std::ignore) = adc(t6, 0, carry);

        u0 = r1, u1 = r2, u2 = r3, u3 = r4, u4 = r5, u5 = r6;
    }

    return Fp({u0, u1, u2, u3, u4, u5}).subtract_modulus();
}

Fp &Fp::operator+=(const Fp &rhs) {
    uint64_t carry = 0;
    uint64_t d0 = 0, d1 = 0, d2 = 0, d3 = 0, d4 = 0, d5 = 0;

    std::tie(d0, carry) = adc(this->data[0], rhs.data[0], 0);
    std::tie(d1, carry) = adc(this->data[1], rhs.data[1], carry);
    std::tie(d2, carry) = adc(this->data[2], rhs.data[2], carry);
    std::tie(d3, carry) = adc(this->data[3], rhs.data[3], carry);
    std::tie(d4, carry) = adc(this->data[4], rhs.data[4], carry);
    std::tie(d5, carry) = adc(this->data[5], rhs.data[5], carry);

    *this = Fp({d0, d1, d2, d3, d4, d5}).subtract_modulus();
    return *this;
}

Fp &Fp::operator-=(const Fp &rhs) {
    *this += -rhs;
    return *this;
}

Fp &Fp::operator*=(const Fp &rhs) {
    uint64_t carry = 0;
    uint64_t t0 = 0, t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5 = 0, t6 = 0, t7 = 0, t8 = 0, t9 = 0, t10 = 0, t11 = 0;

    std::tie(t0, carry) = mac(0, this->data[0], rhs.data[0], 0);
    std::tie(t1, carry) = mac(0, this->data[0], rhs.data[1], carry);
    std::tie(t2, carry) = mac(0, this->data[0], rhs.data[2], carry);
    std::tie(t3, carry) = mac(0, this->data[0], rhs.data[3], carry);
    std::tie(t4, carry) = mac(0, this->data[0], rhs.data[4], carry);
    std::tie(t5, t6) = mac(0, this->data[0], rhs.data[5], carry);

    std::tie(t1, carry) = mac(t1, this->data[1], rhs.data[0], 0);
    std::tie(t2, carry) = mac(t2, this->data[1], rhs.data[1], carry);
    std::tie(t3, carry) = mac(t3, this->data[1], rhs.data[2], carry);
    std::tie(t4, carry) = mac(t4, this->data[1], rhs.data[3], carry);
    std::tie(t5, carry) = mac(t5, this->data[1], rhs.data[4], carry);
    std::tie(t6, t7) = mac(t6, this->data[1], rhs.data[5], carry);

    std::tie(t2, carry) = mac(t2, this->data[2], rhs.data[0], 0);
    std::tie(t3, carry) = mac(t3, this->data[2], rhs.data[1], carry);
    std::tie(t4, carry) = mac(t4, this->data[2], rhs.data[2], carry);
    std::tie(t5, carry) = mac(t5, this->data[2], rhs.data[3], carry);
    std::tie(t6, carry) = mac(t6, this->data[2], rhs.data[4], carry);
    std::tie(t7, t8) = mac(t7, this->data[2], rhs.data[5], carry);

    std::tie(t3, carry) = mac(t3, this->data[3], rhs.data[0], 0);
    std::tie(t4, carry) = mac(t4, this->data[3], rhs.data[1], carry);
    std::tie(t5, carry) = mac(t5, this->data[3], rhs.data[2], carry);
    std::tie(t6, carry) = mac(t6, this->data[3], rhs.data[3], carry);
    std::tie(t7, carry) = mac(t7, this->data[3], rhs.data[4], carry);
    std::tie(t8, t9) = mac(t8, this->data[3], rhs.data[5], carry);

    std::tie(t4, carry) = mac(t4, this->data[4], rhs.data[0], 0);
    std::tie(t5, carry) = mac(t5, this->data[4], rhs.data[1], carry);
    std::tie(t6, carry) = mac(t6, this->data[4], rhs.data[2], carry);
    std::tie(t7, carry) = mac(t7, this->data[4], rhs.data[3], carry);
    std::tie(t8, carry) = mac(t8, this->data[4], rhs.data[4], carry);
    std::tie(t9, t10) = mac(t9, this->data[4], rhs.data[5], carry);

    std::tie(t5, carry) = mac(t5, this->data[5], rhs.data[0], 0);
    std::tie(t6, carry) = mac(t6, this->data[5], rhs.data[1], carry);
    std::tie(t7, carry) = mac(t7, this->data[5], rhs.data[2], carry);
    std::tie(t8, carry) = mac(t8, this->data[5], rhs.data[3], carry);
    std::tie(t9, carry) = mac(t9, this->data[5], rhs.data[4], carry);
    std::tie(t10, t11) = mac(t10, this->data[5], rhs.data[5], carry);

    uint64_t ts[12] = {
            t0, t1, t2, t3,
            t4, t5, t6, t7,
            t8, t9, t10, t11
    };

    *this = Fp::montgomery_reduce(ts);
    return *this;
}

Fp &Fp::operator=(const Fp &rhs) {
    if (this == &rhs) return *this;
    for (int i = 0; i < std::size(this->data); ++i) {
        this->data[i] = rhs.data[i];
    }
    return *this;
}

Fp Fp::operator-() const {
    uint64_t borrow = 0;
    uint64_t d0 = 0, d1 = 0, d2 = 0, d3 = 0, d4 = 0, d5 = 0;

    std::tie(d0, borrow) = sbb(MODULUS[0], this->data[0], 0);
    std::tie(d1, borrow) = sbb(MODULUS[1], this->data[1], borrow);
    std::tie(d2, borrow) = sbb(MODULUS[2], this->data[2], borrow);
    std::tie(d3, borrow) = sbb(MODULUS[3], this->data[3], borrow);
    std::tie(d4, borrow) = sbb(MODULUS[4], this->data[4], borrow);
    std::tie(d5, borrow) = sbb(MODULUS[5], this->data[5], borrow);

    bool dec = (this->data[0] | this->data[1] | this->data[2] | this->data[3] | this->data[4] | this->data[5]) == 0;
    uint64_t mask = static_cast<uint64_t>(dec) - 1;

    return Fp({d0 & mask, d1 & mask, d2 & mask, d3 & mask, d4 & mask, d5 & mask});
}

/// Reduces a big-endian 64-bit limb representation of a 768-bit number.
Fp Fp::reduce(std::span<uint64_t> limbs) {
    assert(limbs.size() == Fp::WIDTH * 2);

    Fp d1({limbs[11], limbs[10], limbs[9], limbs[8], limbs[7], limbs[6]});
    Fp d0({limbs[5], limbs[4], limbs[3], limbs[2], limbs[1], limbs[0]});
    return d0 * R2 + d1 * R3;
}
