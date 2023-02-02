#include "pairing/pairing.h"

#include <tuple>
#include <vector>

#include "field/fp12.h"

namespace bls12_381::pairing {

template<>
void miller_loop(MillerLoopDriver<void> &driver) {
    bool found_one = false;
    for (int i = 63; i >= 0; --i) {
        bool bit = (((group::constant::BLS_X >> 1) >> i) & 1) == 1;
        if (!found_one) {
            found_one = bit;
            continue;
        }
        driver.doubling_step();
        if (bit) driver.addition_step();
        driver.square_output();
    }
    driver.doubling_step();
    driver.conjugate();
}

field::Fp12 ell(const field::Fp12 &f, const std::tuple<field::Fp2, field::Fp2, field::Fp2> &coefficients,
                const group::G1Affine &p) {
    field::Fp c00 = std::get<0>(coefficients).get_c0() * p.get_y();
    field::Fp c01 = std::get<0>(coefficients).get_c1() * p.get_y();
    field::Fp c10 = std::get<1>(coefficients).get_c0() * p.get_x();
    field::Fp c11 = std::get<1>(coefficients).get_c1() * p.get_x();

    return f.mul_by_fp2(std::get<2>(coefficients), field::Fp2{c10, c11}, field::Fp2{c00, c01});
}

struct HelperMulti : pairing::MillerLoopDriver<field::Fp12> {
    std::vector<std::tuple<group::G1Affine, group::G2Prepared>> terms;
    size_t index;

    HelperMulti(const std::vector<std::tuple<group::G1Affine, group::G2Prepared>> &terms, size_t index) : terms{terms},
                                                                                                          index{index} {}

    field::Fp12 doubling_step(field::Fp12 &f) override {
        size_t idx = this->index;
        for (auto term: this->terms) {
            bool either_identity = std::get<0>(term).is_identity() | std::get<1>(term).is_identity();
            field::Fp12 new_f = pairing::ell(f, std::get<1>(term).get_coeffs()[idx], std::get<0>(term));
            f = either_identity ? f : new_f;
        }
        this->index += 1;
        return f;
    }

    field::Fp12 addition_step(field::Fp12 &f) override {
        size_t idx = this->index;
        for (auto term: this->terms) {
            bool either_identity = std::get<0>(term).is_identity() | std::get<1>(term).is_identity();
            field::Fp12 new_f = pairing::ell(f, std::get<1>(term).get_coeffs()[idx], std::get<0>(term));
            f = either_identity ? f : new_f;
        }
        this->index += 1;
        return f;
    }

    field::Fp12 square_output(field::Fp12 &f) override {
        return f.square();
    }

    field::Fp12 conjugate(field::Fp12 &f) override {
        return f.conjugate();
    }

    field::Fp12 one() override {
        return field::Fp12::one();
    }
};

MillerLoopResult multi_miller_loop(const std::vector<std::tuple<group::G1Affine, group::G2Prepared>> &terms) {
    HelperMulti helper{terms, 0};
    field::Fp12 temp = miller_loop(helper);
    return MillerLoopResult{temp};
}

struct HelperPairing : pairing::MillerLoopDriver<field::Fp12> {
    group::G2Projective current;
    group::G2Affine base;
    group::G1Affine p;

    HelperPairing(const group::G2Projective &current,
                  const group::G2Affine &base, const group::G1Affine &p)
            : current{current}, base{base}, p{p} {}

    field::Fp12 doubling_step(field::Fp12 &f) override {
        auto coeffs = pairing::doubling_step(this->current);
        return pairing::ell(f, coeffs, this->p);
    }

    field::Fp12 addition_step(field::Fp12 &f) override {
        auto coeffs = pairing::addition_step(this->current, this->base);
        return pairing::ell(f, coeffs, this->p);
    }

    field::Fp12 square_output(field::Fp12 &f) override {
        return f.square();
    }

    field::Fp12 conjugate(field::Fp12 &f) override {
        return f.conjugate();
    }

    field::Fp12 one() override {
        return field::Fp12::one();
    }
};

group::Gt pairings(const group::G1Affine &p, const group::G2Affine &q) {
    bool either_identity = p.is_identity() | q.is_identity();
    group::G1Affine temp_p = either_identity ? group::G1Affine::generator() : p;
    group::G2Affine temp_q = either_identity ? group::G2Affine::generator() : q;

    HelperPairing helper{group::G2Projective{temp_q}, temp_q, temp_p};
    field::Fp12 temp = pairing::miller_loop(helper);
    MillerLoopResult temp2{either_identity ? field::Fp12::one() : temp};
    return temp2.final_exponentiation();
}

} // namespace bls12_381::pairing