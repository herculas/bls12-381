#include "group/g2_prepared.h"

#include <cassert>
#include <utility>
#include <vector>

#include "utils/bit.h"

#include "field/fp.h"
#include "field/fp2.h"
#include "group/g2_affine.h"
#include "group/g2_projective.h"
#include "pairing/pairing.h"

namespace bls12_381::group {

using rng::util::bit::to_le_bytes;
using rng::util::bit::from_le_bytes;

using field::Fp;
using field::Fp2;
using pairing::MillerLoopDriver;

using coeff_vec = std::vector<std::tuple<Fp2, Fp2, Fp2>>;

G2Prepared::G2Prepared(bool infinity, coeff_vec coefficients)
        : infinity{infinity}, coefficients{std::move(coefficients)} {}

struct Helper : MillerLoopDriver<void> {
    G2Projective current;
    G2Affine base;
    coeff_vec coefficients;

    Helper(G2Projective current, G2Affine base, coeff_vec coefficients)
            : current{std::move(current)}, base{std::move(base)}, coefficients{std::move(coefficients)} {}

    void doubling_step() override {
        auto coeffs = pairing::doubling_step(this->current);
        this->coefficients.push_back(coeffs);
    }

    void addition_step() override {
        auto coeffs = pairing::addition_step(this->current, this->base);
        this->coefficients.push_back(coeffs);
    }

    void square_output() override {}

    void conjugate() override {}

    void one() override {}
};

G2Prepared::G2Prepared(const G2Affine &point) : infinity{}, coefficients{} {
    bool is_identity = point.is_identity();
    G2Affine q = (is_identity) ? G2Affine::generator() : point;
    coeff_vec coeffs_temp;
    coeffs_temp.reserve(68);
    Helper helper{G2Projective{q}, q, coeffs_temp};
    pairing::miller_loop(helper);

    assert(helper.coefficients.size() == 68);
    *this = G2Prepared{is_identity, helper.coefficients};
}

G2Prepared::G2Prepared(G2Affine &&point) : infinity{}, coefficients{} {
    bool is_identity = point.is_identity();
    G2Affine q = (is_identity) ? G2Affine::generator() : point;
    coeff_vec coeffs_temp;
    coeffs_temp.reserve(68);
    Helper helper{G2Projective{q}, q, coeffs_temp};
    pairing::miller_loop(helper);

    assert(helper.coefficients.size() == 68);
    *this = G2Prepared{is_identity, helper.coefficients};
}

bool G2Prepared::is_identity() const {
    return this->infinity;
}

const coeff_vec &G2Prepared::get_coeffs() const {
    return this->coefficients;
}

auto G2Prepared::from_slice_unchecked(const std::vector<uint8_t> &bytes) -> G2Prepared {
    const size_t unit_size = Fp::BYTE_SIZE * 2 * 3;
    coeff_vec coeffs{};
    coeffs.reserve(bytes.size() / unit_size);

    for (int i = 0; i < bytes.size(); i += unit_size) {
        std::array<uint64_t, Fp::WIDTH> a_c0_data{};
        std::array<uint64_t, Fp::WIDTH> a_c1_data{};
        std::array<uint64_t, Fp::WIDTH> b_c0_data{};
        std::array<uint64_t, Fp::WIDTH> b_c1_data{};
        std::array<uint64_t, Fp::WIDTH> c_c0_data{};
        std::array<uint64_t, Fp::WIDTH> c_c1_data{};

        for (int j = 0; j < Fp::WIDTH; ++j) {
            std::array<uint8_t, 8> a_c0_bytes{};
            std::array<uint8_t, 8> a_c1_bytes{};
            std::array<uint8_t, 8> b_c0_bytes{};
            std::array<uint8_t, 8> b_c1_bytes{};
            std::array<uint8_t, 8> c_c0_bytes{};
            std::array<uint8_t, 8> c_c1_bytes{};

            std::copy(bytes.begin() + i + j * 48, bytes.begin() + i + j * 48 + 8, a_c0_bytes.begin());
            std::copy(bytes.begin() + i + j * 48 + 8, bytes.begin() + i + j * 48 + 16, a_c1_bytes.begin());
            std::copy(bytes.begin() + i + j * 48 + 16, bytes.begin() + i + j * 48 + 24, b_c0_bytes.begin());
            std::copy(bytes.begin() + i + j * 48 + 24, bytes.begin() + i + j * 48 + 32, b_c1_bytes.begin());
            std::copy(bytes.begin() + i + j * 48 + 32, bytes.begin() + i + j * 48 + 40, c_c0_bytes.begin());
            std::copy(bytes.begin() + i + j * 48 + 40, bytes.begin() + i + j * 48 + 48, c_c1_bytes.begin());

            a_c0_data[j] = from_le_bytes<uint64_t>(a_c0_bytes);
            a_c1_data[j] = from_le_bytes<uint64_t>(a_c1_bytes);
            b_c0_data[j] = from_le_bytes<uint64_t>(b_c0_bytes);
            b_c1_data[j] = from_le_bytes<uint64_t>(b_c1_bytes);
            c_c0_data[j] = from_le_bytes<uint64_t>(c_c0_bytes);
            c_c1_data[j] = from_le_bytes<uint64_t>(c_c1_bytes);
        }

        coeffs.push_back(
                {
                        Fp2{Fp{a_c0_data}, Fp{a_c1_data}},
                        Fp2{Fp{b_c0_data}, Fp{b_c1_data}},
                        Fp2{Fp{c_c0_data}, Fp{c_c1_data}}
                }
        );
    }
    return G2Prepared{false, coeffs};
}

std::vector<uint8_t> G2Prepared::to_raw_bytes() const {
    const size_t unit_size = Fp::BYTE_SIZE * 2 * 3;
    const size_t size = unit_size * this->coefficients.size();
    std::vector<uint8_t> bytes{};
    bytes.reserve(size);

    for (const auto &[a, b, c]: this->coefficients) {
        for (int i = 0; i < Fp::WIDTH; ++i) {
            const std::array<uint8_t, 8> a_c0_bytes = to_le_bytes<uint64_t>(a.get_c0().get_data()[i]);
            const std::array<uint8_t, 8> a_c1_bytes = to_le_bytes<uint64_t>(a.get_c1().get_data()[i]);
            const std::array<uint8_t, 8> b_c0_bytes = to_le_bytes<uint64_t>(b.get_c0().get_data()[i]);
            const std::array<uint8_t, 8> b_c1_bytes = to_le_bytes<uint64_t>(b.get_c1().get_data()[i]);
            const std::array<uint8_t, 8> c_c0_bytes = to_le_bytes<uint64_t>(c.get_c0().get_data()[i]);
            const std::array<uint8_t, 8> c_c1_bytes = to_le_bytes<uint64_t>(c.get_c1().get_data()[i]);

            bytes.insert(bytes.end(), a_c0_bytes.begin(), a_c0_bytes.end());
            bytes.insert(bytes.end(), a_c1_bytes.begin(), a_c1_bytes.end());
            bytes.insert(bytes.end(), b_c0_bytes.begin(), b_c0_bytes.end());
            bytes.insert(bytes.end(), b_c1_bytes.begin(), b_c1_bytes.end());
            bytes.insert(bytes.end(), c_c0_bytes.begin(), c_c0_bytes.end());
            bytes.insert(bytes.end(), c_c1_bytes.begin(), c_c1_bytes.end());
        }
    }
    return bytes;
}

} // namespace bls12_381::group

