#ifndef BLS12_381_RANDOM_H
#define BLS12_381_RANDOM_H

#include <cstdint>
#include <random>

uint64_t getRandom(uint64_t max = std::numeric_limits<uint64_t>::max()) noexcept;

#endif //BLS12_381_RANDOM_H
