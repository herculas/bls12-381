//#include "scalar/multi.h"
//
//#include <cassert>
//#include <cstddef>
//#include <cstdint>
//#include <array>
//
//namespace bls12_381::scalar::multi {
//
//using group::G1Affine;
//using group::G1Projective;
//
//auto to_radix_2w_size_hint(size_t w) -> size_t {
//    assert(w >= 6);
//    assert(w <= 8);
//
//    size_t digits_count;
//    switch (w) {
//        case 6:
//        case 7:
//            digits_count = (256 + w - 1) / w;
//            break;
//        case 8:
//            digits_count = (256 + w - 1) / w + 1;
//            break;
//        default:
//            assert(false);
//    }
//    assert(digits_count <= 43);
//    return digits_count;
//}
//
//auto to_radix_2w(const Scalar &scalar, size_t w) -> std::array<int8_t, 43> {
//    assert(w >= 6);
//    assert(w <= 8);
//
//    const auto scalar64x4 = scalar.to_bytes();
//
//    const uint64_t radix = 1 << w;
//    const uint64_t window_mask = radix - 1;
//
//    uint64_t carry = 0;
//    std::array<int8_t, 43> digits = {0};
//    const size_t digits_count = (256 + w - 1) / w;
//    for (int i = 0; i < digits_count; ++i) {
//        const size_t bit_offset = i * w;
//        const size_t u64_index = bit_offset / 64;
//        const size_t bit_index = bit_offset % 64;
//        uint64_t bit_buf;
//        if (bit_index < 64 - w || u64_index == 3) {
//            bit_buf = scalar64x4[u64_index] >> bit_index;
//        } else {
//            bit_buf = (scalar64x4[u64_index] >> bit_index) | (scalar64x4[u64_index + 1] << (64 - bit_index));
//        }
//        const uint64_t coeff = carry + (bit_buf & window_mask);
//        carry = (coeff + (radix / 2)) >> w;
//        digits[i] = static_cast<int8_t>(static_cast<int64_t>(coeff) - static_cast<int64_t>(carry << w));
//    }
//
//    if (w == 8)
//        digits[digits_count] += static_cast<int8_t>(carry); // NOLINT(cppcoreguidelines-narrowing-conversions)
//    else
//        digits[digits_count - 1] += static_cast<int8_t>(carry << w); // NOLINT(cppcoreguidelines-narrowing-conversions)
//
//    return digits;
//}
//
//auto log2(size_t x) -> uint32_t {
//    if (x <= 1)return 0;
//    const uint32_t n = __builtin_clzll(x);
//    return static_cast<uint32_t>(sizeof(uint64_t)) * 8 - n;
//}
//
//auto ln_without_floats(size_t a) -> size_t {
//    return static_cast<size_t>(log2(a) * 69 / 100);
//}
//
//auto pippenger(const std::vector<G1Projective> &points, const std::vector<Scalar> &scalars) -> G1Projective {
//
//
////    const size_t size = scalars.size();
////    const size_t w = size < 500 ? 6 : (size < 800 ? 7 : 8);
////
////    const size_t max_digit = 1 << w;
////    const size_t digits_count = to_radix_2w_size_hint(w);
////    const size_t buckets_count = max_digit / 2;
////
////    std::vector<std::array<int8_t, 43>> scalar_digits{};
////    scalar_digits.reserve(scalars.size());
////    for (const auto &scalar: scalars)
////        scalar_digits.push_back(to_radix_2w(scalar, w));
////
////    std::vector<G1Projective> buckets{buckets_count, G1Projective::identity()};
////
////    G1Projective columns{};
////
////    for (int64_t digit_index = static_cast<int64_t>(digits_count) - 1; digit_index >= 0; --digit_index) {
////        for (int64_t i = 0; i < buckets_count; i++)
////            buckets[digit_index] = G1Projective::identity();
////
////        const size_t count_scalar = scalars.size();
////        for (int j = 0; j < count_scalar; ++j) {
////            const auto digits = scalar_digits[j];
////            const auto &point = points[j];
////
////            const auto digit = static_cast<int16_t>(digits[digit_index]);
////            if (digit > 0) {
////                const auto b = static_cast<size_t>(digit - 1);
////                buckets[b] = buckets[b] + point;
////            } else if (digit < 0) {
////                const auto b = static_cast<size_t>(-digit - 1);
////                buckets[b] = buckets[b] - point;
////            }
////        }
////
////        G1Projective buckets_intermediate_sum = buckets[buckets_count - 1];
////        G1Projective buckets_sum = buckets[buckets_count - 1];
////        for (int64_t i = static_cast<int64_t>(digits_count) - 2; i >= 0; --i) {
////            buckets_intermediate_sum += buckets[i];
////            buckets_sum += buckets_intermediate_sum;
////        }
////        columns = buckets_sum;
////    }
//
//
//}
//
//auto mul_by_pow_2(const G1Projective &point, uint32_t k) -> G1Projective {
//    return group::G1Projective();
//}
//
//auto msm_variable_base(const std::vector<G1Affine> &points, const std::vector<Scalar> &scalars) -> G1Projective {
//    return group::G1Projective();
//}
//
//} // namespace bls12_381::scalar::multi