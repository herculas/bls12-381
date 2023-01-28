#ifndef BLS12_381_ENCODE_H
#define BLS12_381_ENCODE_H

#include <string>
#include <span>

namespace bls12_381::util::encoding {

std::string hexStr(std::span<uint8_t> span);

} // namespace bls12_381::util::encoding

#endif //BLS12_381_ENCODE_H