#ifndef BLS12_381_ENCODE_H
#define BLS12_381_ENCODE_H

#include <string>
#include <span>

namespace bls12_381::util::encoding {

constexpr char hex_map[16] = { // NOLINT(cppcoreguidelines-avoid-c-arrays)
        '0', '1', '2', '3',
        '4', '5', '6', '7',
        '8', '9', 'a', 'b',
        'c', 'd', 'e', 'f'
};

inline std::string hex_str(const std::span<uint8_t> digits) {
    std::string result;
    result.reserve(digits.size() * 2);
    for (uint8_t const digit: digits) {
        result.push_back(hex_map[digit >> 4]);
        result.push_back(hex_map[digit & 15]);
    }
    return result;
}

} // namespace bls12_381::util::encoding

#endif //BLS12_381_ENCODE_H