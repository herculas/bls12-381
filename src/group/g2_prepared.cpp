#include "group/g2_prepared.h"

#include <vector>

#include "field/fp2.h"
#include "group/g2_affine.h"
#include "group/g2_projective.h"
#include "pairing/miller_loop_driver.h"

namespace bls12_381::group {

struct Helper : pairing::MillerLoopDriver<void> {
    G2Projective current;
    G2Affine base;
    std::vector<std::array<field::Fp2, 3>> coefficients;

    Helper(const G2Projective &current, const G2Affine &base,
           const std::vector<std::array<field::Fp2, 3>> &coefficients) : current{current}, base{base},
                                                                         coefficients{coefficients} {}

    void doubling_step() override {
        auto coeffs = bls12_381::pairing::doubling_step(this->current);
        this->coefficients.push_back(coeffs);
    }

    void addition_step() override {
        auto coeffs = bls12_381::pairing::addition_step(this->current, this->base);
        this->coefficients.push_back(coeffs);
    }

    void square_output() override {}

    void conjugate() override {}

    void one() override {}
};

G2Prepared::G2Prepared(bool infinity,
                       const std::vector<std::array<field::Fp2, 3>> &coefficients) : infinity{infinity},
                                                                                     coefficients{coefficients} {}


G2Prepared::G2Prepared(const G2Affine &point) : infinity{}, coefficients{} {
    bool is_identity = point.is_identity();
    G2Affine q = (is_identity) ? G2Affine::generator() : point;
    Helper helper{G2Projective{q}, q, std::vector<std::array<field::Fp2, 3>>(68)};
    bls12_381::pairing::miller_loop(helper);

    assert(helper.coefficients.size() == 68);
    *this = G2Prepared{is_identity, helper.coefficients};
}

G2Prepared::G2Prepared(G2Affine &&point) : infinity{}, coefficients{} {
    bool is_identity = point.is_identity();
    G2Affine q = (is_identity) ? G2Affine::generator() : point;
    Helper helper{G2Projective{q}, q, std::vector<std::array<field::Fp2, 3>>(68)};
    bls12_381::pairing::miller_loop(helper);

    assert(helper.coefficients.size() == 68);
    *this = G2Prepared{is_identity, helper.coefficients};
}

} // namespace bls12_381::group

