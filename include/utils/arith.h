#ifndef BLS12_381_ARITH_H
#define BLS12_381_ARITH_H

#include <cstdint>

namespace bls12_381::util::arithmetic {

/**
 * @brief Computes a + b + carry, returning the result and updating the carry.
 * @param a A 64-bit unsigned integer.
 * @param b A 64-bit unsigned integer.
 * @param carry The carry, a 64-bit unsigned integer, could be updated by this function.
 * @return The result of the addition.
 */
constexpr inline uint64_t adc(uint64_t a, uint64_t b, uint64_t &carry) {
    const __uint128_t res = static_cast<__uint128_t>(a) + static_cast<__uint128_t>(b) + static_cast<__uint128_t>(carry);
    carry = static_cast<uint64_t>((res >> 64));
    return static_cast<uint64_t>(res);
}

/**
 * @brief Computes a - (b + borrow), returning the result and updating the borrow.
 * @param a A 64-bit unsigned integer.
 * @param b A 64-bit unsigned integer.
 * @param borrow The borrow, a 64-bit unsigned integer, could be updated by this function.
 * @return The result of the subtraction.
 */
constexpr inline uint64_t sbb(uint64_t a, uint64_t b, uint64_t &borrow) {
    const __uint128_t res =
            static_cast<__uint128_t>(a) - static_cast<__uint128_t>(b) - static_cast<__uint128_t>(borrow >> 63);
    borrow = static_cast<uint64_t>(res >> 64);
    return static_cast<uint64_t>(res);
}

/**
 * @brief Computes a + (b * c) + carry, returning the result and updating the carry.
 * @param a A 64-bit unsigned integer.
 * @param b A 64-bit unsigned integer.
 * @param c A 64-bit unsigned integer.
 * @param carry The carry, a 64-bit unsigned integer, could be updated by this function.
 * @return The result of the calculation.
 */
constexpr inline uint64_t mac(uint64_t a, uint64_t b, uint64_t c, uint64_t &carry) {
    const __uint128_t res = static_cast<__uint128_t>(a) + (static_cast<__uint128_t>(b) * static_cast<__uint128_t>(c)) +
                            static_cast<__uint128_t>(carry);
    carry = static_cast<uint64_t>(res >> 64);
    return static_cast<uint64_t>(res);
}

} // namespace bls12_381::util::arithmetic

#endif //BLS12_381_ARITH_H