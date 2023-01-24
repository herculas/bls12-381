#ifndef BLS12_381_ARITH_H
#define BLS12_381_ARITH_H

#include <cstdint>
#include <tuple>

/// Computes a + b + carry, returning the result and the new carry over.
std::tuple<uint64_t, uint64_t> adc(uint64_t a, uint64_t b, uint64_t carry);

/// Computes a - (b + borrow), returning the result and the new borrow.
std::tuple<uint64_t, uint64_t> sbb(uint64_t a, uint64_t b, uint64_t borrow);

/// Computes a + (b * c) + carry, returning the result and the new carry over.
std::tuple<uint64_t, uint64_t> mac(uint64_t a, uint64_t b, uint64_t c, uint64_t carry);

/// Creates an uint64_t value from a byte array in big endian.
uint64_t be_bytes_to_uint64(const uint8_t *bytes);

uint8_t *uint64_to_be_bytes(uint64_t value, uint8_t *bytes);

#endif //BLS12_381_ARITH_H
