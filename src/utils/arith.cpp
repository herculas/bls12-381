#include "utils/arith.h"
#include <tuple>

std::tuple<uint64_t, uint64_t> adc(uint64_t a, uint64_t b, uint64_t carry) {
    auto res = static_cast<__uint128_t>(a) + static_cast<__uint128_t>(b) + static_cast<__uint128_t>(carry);
    return {static_cast<uint64_t >(res), static_cast<uint64_t>(res >> 64)};
}

std::tuple<uint64_t, uint64_t> sbb(uint64_t a, uint64_t b, uint64_t borrow) {
    auto res = static_cast<__uint128_t>(a) - static_cast<__uint128_t>(b) - static_cast<__uint128_t>(borrow >> 63);
    return {static_cast<uint64_t>(res), static_cast<uint64_t>(res >> 64)};
}

std::tuple<uint64_t, uint64_t> mac(uint64_t a, uint64_t b, uint64_t c, uint64_t carry) {
    auto res = static_cast<__uint128_t>(a) + (static_cast<__uint128_t>(b) * static_cast<__uint128_t>(c)) +
               static_cast<__uint128_t>(carry);
    return {static_cast<uint64_t>(res), static_cast<uint64_t>(res >> 64)};
}

uint64_t be_bytes_to_uint64(const uint8_t *bytes) {
    uint64_t sum = 0;
    for (int i = 0; i < 8; ++i)
        sum += static_cast<uint64_t>(bytes[i] & 0x00ff) << (56 - i * 8);
    return sum;
}

uint8_t *uint64_to_be_bytes(uint64_t value, uint8_t *const bytes) {
    for (int i = 0; i < 8; ++i)
        bytes[i] = static_cast<uint8_t>(value >> (56 - i * 8)) & 0x00ff;
    return bytes;
}