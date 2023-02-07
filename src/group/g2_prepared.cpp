#include "group/g2_prepared.h"

#include <cassert>
#include <vector>

#include "field/fp2.h"
#include "group/g2_affine.h"
#include "group/g2_projective.h"
#include "pairing/pairing.h"

namespace bls12_381::group {

G2Prepared::G2Prepared(bool infinity, const std::vector<std::tuple<field::Fp2, field::Fp2, field::Fp2>> &coefficients)
        : infinity{infinity}, coefficients{coefficients} {}

struct Helper : pairing::MillerLoopDriver<void> {
    G2Projective current;
    G2Affine base;
    std::vector<std::tuple<field::Fp2, field::Fp2, field::Fp2>> coefficients;

    Helper(G2Projective current, G2Affine base,
           const std::vector<std::tuple<field::Fp2, field::Fp2, field::Fp2>> &coefficients)
            : current{std::move(current)}, base{std::move(base)}, coefficients{coefficients} {}

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
    std::vector<std::tuple<field::Fp2, field::Fp2, field::Fp2>> coeffs_temp;
    coeffs_temp.reserve(68);
    Helper helper{G2Projective{q}, q, coeffs_temp};
    pairing::miller_loop(helper);

    assert(helper.coefficients.size() == 68);
    *this = G2Prepared{is_identity, helper.coefficients};
}

G2Prepared::G2Prepared(G2Affine &&point) : infinity{}, coefficients{} {
    bool is_identity = point.is_identity();
    G2Affine q = (is_identity) ? G2Affine::generator() : point;
    std::vector<std::tuple<field::Fp2, field::Fp2, field::Fp2>> coeffs_temp;
    coeffs_temp.reserve(68);
    Helper helper{G2Projective{q}, q, coeffs_temp};
    pairing::miller_loop(helper);

    assert(helper.coefficients.size() == 68);
    *this = G2Prepared{is_identity, helper.coefficients};
}

bool G2Prepared::is_identity() const {
    return this->infinity;
}

std::vector<std::tuple<field::Fp2, field::Fp2, field::Fp2>> G2Prepared::get_coeffs() const {
    return this->coefficients;
}

} // namespace bls12_381::group

