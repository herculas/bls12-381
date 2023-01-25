#include <gtest/gtest.h>
#include <span>
#include "utils/arith.h"
#include "utils/bit.h"

TEST(TestUtil, BytesToUint64BE) {
    uint64_t expected = 0x0001020304050607;
    uint8_t array[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    uint64_t result = be_bytes_to_uint64(array);

    EXPECT_EQ(result, expected);
}

TEST(TestUtil, Uint64BEToBytes) {
    uint64_t value = 0x0001020304050607;
    uint8_t array[8];
    uint64_to_be_bytes(value, array);
    const uint8_t array_expected[8] = {0, 1, 2, 3, 4, 5, 6, 7};

    for (int i = 0; i < std::size(array); ++i) {
        EXPECT_EQ(array[i], array_expected[i]);
    }
}