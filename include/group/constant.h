#ifndef BLS12_381_GROUP_CONSTANT_H
#define BLS12_381_GROUP_CONSTANT_H

#include <cstdint>

namespace bls12_381::group::constant {

/// The BLS parameter x for BLS12-381 is -0xd201000000010000.
const uint64_t BLS_X = 0xd201000000010000;
const bool BLS_X_IS_NEGATIVE = true;

} // namespace bls12_381::group::constant

#endif //BLS12_381_GROUP_CONSTANT_H