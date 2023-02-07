#ifndef BLS12_381_ENCODE_H
#define BLS12_381_ENCODE_H

#include <string>
#include <span>

namespace bls12_381::util::encoding {

inline std::string hex_str(const std::span<uint8_t> span) {
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

} // namespace bls12_381::util::encoding

#endif //BLS12_381_ENCODE_H