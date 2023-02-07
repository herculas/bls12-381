#ifndef BLS12_381_BIT_H
#define BLS12_381_BIT_H

#include <cstdint>
#include <array>

namespace bls12_381::util::bit_operation {

/// Converts a byte array to a 64-bit unsigned integer.
constexpr inline uint64_t be_bytes_to_uint64(const std::array<uint8_t, sizeof(uint64_t)> bytes) {
    uint64_t sum = 0;
    for (int i = 0; i < sizeof(uint64_t); ++i)
        sum += static_cast<uint64_t>(bytes[i] & 0x00ff) << (56 - i * 8);
    return sum;
}

constexpr inline uint64_t le_bytes_to_uint64(const std::array<uint8_t, sizeof(uint64_t)> bytes) {
    uint64_t sum = 0;
    for (int i = 0; i < sizeof(uint64_t); ++i)
        sum += static_cast<uint64_t>(bytes[i] & 0x00ff) << (i * 8);
    return sum;
}

/// Converts a 64-bit unsigned integer to a byte array.
constexpr inline std::array<uint8_t, sizeof(uint64_t)> uint64_to_be_bytes(uint64_t value) {
    std::array<uint8_t, sizeof(uint64_t)> bytes{0};
    for (int i = 0; i < sizeof(uint64_t); ++i)
        bytes[i] = static_cast<uint8_t>(value >> (56 - i * 8)) & 0x00ff;
    return bytes;
}

constexpr inline std::array<uint8_t, sizeof(uint64_t)> uint64_to_le_bytes(uint64_t value) {
    std::array<uint8_t, sizeof(uint64_t)> bytes{0};
    for (int i = 0; i < sizeof(uint64_t); ++i)
        bytes[i] = static_cast<uint8_t>(value >> (i * 8)) & 0x00ff;
    return bytes;
}

} // namespace bls12_381::util::bit_operation

#endif //BLS12_381_BIT_H