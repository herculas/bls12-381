#ifndef BLS12_381_BIT_H
#define BLS12_381_BIT_H

#include <cstdint>
#include <span>

/// Converts a byte array (in big endian) to a 64-bit unsigned integer.
uint64_t be_bytes_to_uint64(std::span<uint8_t> bytes);

/// Converts a 64-bit unsigned integer to a byte array.
uint8_t *uint64_to_be_bytes(uint64_t value, std::span<uint8_t> bytes);

#endif //BLS12_381_BIT_H
