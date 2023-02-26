#ifndef BLS12_381_ARITH_H
#define BLS12_381_ARITH_H

#include <cstdint>

#ifdef _MSC_VER
#include "intrin.h"
#endif

namespace bls12_381::util::arithmetic {

/// Computes a + b + carry, returning the result and updating the carry.
constexpr inline uint64_t adc(uint64_t a, uint64_t b, uint64_t &carry) {
#ifdef __SIZEOF_INT128__
    auto res = static_cast<__uint128_t>(a) + static_cast<__uint128_t>(b) + static_cast<__uint128_t>(carry);
    carry = static_cast<uint64_t>((res >> 64));
    return static_cast<uint64_t>(res);
#else

#ifdef _MSC_VER
    _addcarry_u64()
#endif

#endif
}

/// Computes a - (b + borrow), returning the result and updating the borrow.
constexpr inline uint64_t sbb(uint64_t a, uint64_t b, uint64_t &borrow) {
#ifdef __SIZEOF_INT128__
    auto res = static_cast<__uint128_t>(a) - static_cast<__uint128_t>(b) - static_cast<__uint128_t>(borrow >> 63);
    borrow = static_cast<uint64_t>(res >> 64);
    return static_cast<uint64_t>(res);
#else

#ifdef _MSC_VER
    _subborrow_u64()
#endif

#endif
}

/// Computes a + (b * c) + carry, returning the result and updating the carry.
constexpr inline uint64_t mac(uint64_t a, uint64_t b, uint64_t c, uint64_t &carry) {
#ifdef __SIZEOF_INT128__
    auto res = static_cast<__uint128_t>(a) + (static_cast<__uint128_t>(b) * static_cast<__uint128_t>(c)) +
               static_cast<__uint128_t>(carry);
    carry = static_cast<uint64_t>(res >> 64);
    return static_cast<uint64_t>(res);
#else

#ifdef _MSC_VER
    _mul64() -> _addcarry_u64()
#endif

#endif
}

} // namespace bls12_381::util::arithmetic

#endif //BLS12_381_ARITH_H``````````