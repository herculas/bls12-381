#ifndef BLS12_381_RANDOM_H
#define BLS12_381_RANDOM_H

#include <cstdint>
#include <random>

namespace bls12_381::util::random {

template<typename T>
T getRandom(T max = std::numeric_limits<T>::max()) noexcept {
    std::random_device device;
    std::default_random_engine randomEngine(device());
    std::uniform_int_distribution<uint64_t> distribution{0, max};
    return distribution(randomEngine);
}

} // namespace bls12_381::util::random

#endif //BLS12_381_RANDOM_H\