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
