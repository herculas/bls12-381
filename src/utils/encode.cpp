#include "utils/encode.h"

std::string hexStr(const std::span<uint8_t> span) {
    std::string result;
    static constexpr char hexMap[16] = {
            '0', '1', '2', '3',
            '4', '5', '6', '7',
            '8', '9', 'a', 'b',
            'c', 'd', 'e', 'f'
    };
    result.reserve(span.size() * 2);
    for (uint8_t v: span) {
        result.push_back(hexMap[v >> 4]);
        result.push_back(hexMap[v & 15]);
    }
    return result;
}