#include <gtest/gtest.h>
#include <array>

#include "utils/bit.h"

using rng::util::bit::from_le_bytes;
using rng::util::bit::from_be_bytes;
using rng::util::bit::to_le_bytes;
using rng::util::bit::to_be_bytes;

TEST(TestUtil, BytesToUint64BE) {
    uint64_t const expected = 0x0001020304050607;
    std::array<uint8_t, sizeof(uint64_t)> const array = {0, 1, 2, 3, 4, 5, 6, 7};
    auto result = from_be_bytes<uint64_t>(array);
    EXPECT_EQ(result, expected);
}

TEST(TestUtil, Uint64BEToBytes) {
    uint64_t const value = 0x0001020304050607;
    std::array<uint8_t, sizeof(uint64_t)> array = to_be_bytes<uint64_t>(value);
    std::array<uint8_t, sizeof(uint64_t)> array_expected = {0, 1, 2, 3, 4, 5, 6, 7};
    for (int i = 0; i < std::size(array); ++i)
        EXPECT_EQ(array[i], array_expected[i]);
}

TEST(TestUtil, BytesToUint64LE) {
    uint64_t const expected = 0x0001020304050607;
    std::array<uint8_t, sizeof(uint64_t)> const array = {7, 6, 5, 4, 3, 2, 1, 0};
    auto result = from_le_bytes<uint64_t>(array);
    EXPECT_EQ(result, expected);
}

TEST(TestUtil, Uint64LEToBytes) {
    uint64_t const value = 0x0001020304050607;
    std::array<uint8_t, sizeof(uint64_t)> array = to_le_bytes<uint64_t>(value);
    std::array<uint8_t, sizeof(uint64_t)> array_expected = {7, 6, 5, 4, 3, 2, 1, 0};
    for (int i = 0; i < std::size(array); ++i) EXPECT_EQ(array[i], array_expected[i]);
}