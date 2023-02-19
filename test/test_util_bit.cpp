#include <gtest/gtest.h>
#include "utils/arith.h"
#include "utils/bit.h"

using bls12_381::util::bit_operation::be_bytes_to_uint64;
using bls12_381::util::bit_operation::le_bytes_to_uint64;
using bls12_381::util::bit_operation::uint64_to_be_bytes;
using bls12_381::util::bit_operation::uint64_to_le_bytes;

TEST(TestUtil, BytesToUint64BE) {
    uint64_t expected = 0x0001020304050607;
    std::array<uint8_t, sizeof(uint64_t)> array = {0, 1, 2, 3, 4, 5, 6, 7};
    uint64_t result = be_bytes_to_uint64(array);
    EXPECT_EQ(result, expected);
}

TEST(TestUtil, Uint64BEToBytes) {
    uint64_t value = 0x0001020304050607;
    std::array<uint8_t, sizeof(uint64_t)> array = uint64_to_be_bytes(value);
    std::array<uint8_t, sizeof(uint64_t)> array_expected = {0, 1, 2, 3, 4, 5, 6, 7};
    for (int i = 0; i < std::size(array); ++i)
        EXPECT_EQ(array[i], array_expected[i]);
}

TEST(TestUtil, BytesToUint64LE) {
    uint64_t expected = 0x0001020304050607;
    std::array<uint8_t, sizeof(uint64_t)> array = {7, 6, 5, 4, 3, 2, 1, 0};
    uint64_t result = le_bytes_to_uint64(array);
    EXPECT_EQ(result, expected);
}

TEST(TestUtil, Uint64LEToBytes) {
    uint64_t value = 0x0001020304050607;
    std::array<uint8_t, sizeof(uint64_t)> array = uint64_to_le_bytes(value);
    std::array<uint8_t, sizeof(uint64_t)> array_expected = {7, 6, 5, 4, 3, 2, 1, 0};

    for (int i = 0; i < std::size(array); ++i)
        EXPECT_EQ(array[i], array_expected[i]);
}