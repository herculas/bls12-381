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

#endif //BLS12_381_ARITH_H
