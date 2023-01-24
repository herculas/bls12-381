#include "utils/random.h"

uint64_t getRandom(uint64_t max) noexcept {
    std::random_device device;
    std::default_random_engine randomEngine(device());
    std::uniform_int_distribution<uint64_t> distribution{0, max};
    return distribution(randomEngine);
}