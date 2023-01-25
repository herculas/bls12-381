#include "utils/bit.h"

uint64_t be_bytes_to_uint64(std::span<uint8_t> bytes) {
    uint64_t sum = 0;
    for (int i = 0; i < 8; ++i)
        sum += static_cast<uint64_t>(bytes[i] & 0x00ff) << (56 - i * 8);
    return sum;
}

uint8_t *uint64_to_be_bytes(uint64_t value, std::span<uint8_t> bytes) {
    for (int i = 0; i < 8; ++i)
        bytes[i] = static_cast<uint8_t>(value >> (56 - i * 8)) & 0x00ff;
    return bytes.data();
}