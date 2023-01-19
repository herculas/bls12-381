#include <gtest/gtest.h>
#include "utils/arith.h"

TEST(HelloTest, BasicAssertions) {
    EXPECT_STRNE("hello", "world");
    EXPECT_EQ(7 * 6, 42);
    auto a = adc(std::numeric_limits<uint64_t>::max() / 2, std::numeric_limits<uint64_t>::max() / 2 + 1, 1);
    std::cout << std::get<0>(a) << std::endl << std::get<1>(a) << std::endl;
}
