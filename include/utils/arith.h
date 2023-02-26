#ifndef BLS12_381_ARITH_H
#define BLS12_381_ARITH_H

#include <cstdint>

namespace bls12_381::util::arithmetic {

/// Computes a + b + carry, returning the result and updating the carry.
constexpr inline uint64_t adc(uint64_t a, uint64_t b, uint64_t &carry) {
    auto res = static_cast<__uint128_t>(a) + static_cast<__uint128_t>(b) + static_cast<__uint128_t>(carry);
    carry = static_cast<uint64_t>((res >> 64));
    return static_cast<uint64_t>(res);
}

/// Computes a - (b + borrow), returning the result and updating the borrow.
constexpr inline uint64_t sbb(uint64_t a, uint64_t b, uint64_t &borrow) {
    auto res = static_cast<__uint128_t>(a) - static_cast<__uint128_t>(b) - static_cast<__uint128_t>(borrow >> 63);
    borrow = static_cast<uint64_t>(res >> 64);
    return static_cast<uint64_t>(res);
}

/// Computes a + (b * c) + carry, returning the result and updating the carry.
constexpr inline uint64_t mac(uint64_t a, uint64_t b, uint64_t c, uint64_t &carry) {
    auto res = static_cast<__uint128_t>(a) + (static_cast<__uint128_t>(b) * static_cast<__uint128_t>(c)) +
               static_cast<__uint128_t>(carry);
    carry = static_cast<uint64_t>(res >> 64);
    return static_cast<uint64_t>(res);
}

} // namespace bls12_381::util::arithmetic

#endif //BLS12_381_ARITH_H