#include "pairing/pairing.h"

#include <tuple>
#include <vector>

#include "field/fp12.h"

namespace bls12_381::pairing {

using field::Fp;
using field::Fp2;
using field::Fp12;
using group::G1Affine;
using group::G2Affine;
using group::G2Prepared;
using group::G2Projective;
using group::Gt;
using coeff_t = std::tuple<Fp2, Fp2, Fp2>;

template<>
void miller_loop(MillerLoopDriver<void> &driver) {
    bool found_one = false;
    for (int i = 63; i >= 0; --i) {
        bool const bit = (((group::constant::BLS_X >> 1) >> i) & 1) == 1;
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

Fp12 ell(const Fp12 &f, const coeff_t &coefficients, const G1Affine &p) {
    Fp const c00 = std::get<0>(coefficients).get_c0() * p.get_y();
    Fp const c01 = std::get<0>(coefficients).get_c1() * p.get_y();
    Fp const c10 = std::get<1>(coefficients).get_c0() * p.get_x();
    Fp const c11 = std::get<1>(coefficients).get_c1() * p.get_x();

    return f.mul_by_fp2(std::get<2>(coefficients), Fp2{c10, c11}, Fp2{c00, c01});
}

struct HelperMulti : MillerLoopDriver<Fp12> {
private:
    std::vector<std::tuple<G1Affine, G2Prepared>> terms;
    size_t index;
public:
    HelperMulti(const std::vector<std::tuple<G1Affine, G2Prepared>> &terms, size_t index) : terms{terms},
                                                                                            index{index} {}

    Fp12 doubling_step(Fp12 &f) override {
        size_t const idx = this->index;
        for (auto term: this->terms) {
            bool const either_identity = std::get<0>(term).is_identity() | std::get<1>(term).is_identity();
            Fp12 const new_f = ell(f, std::get<1>(term).get_coeffs()[idx], std::get<0>(term));
            f = either_identity ? f : new_f;
        }
        this->index += 1;
        return f;
    }

    Fp12 addition_step(Fp12 &f) override {
        size_t const idx = this->index;
        for (auto term: this->terms) {
            bool const either_identity = std::get<0>(term).is_identity() | std::get<1>(term).is_identity();
            Fp12 const new_f = ell(f, std::get<1>(term).get_coeffs()[idx], std::get<0>(term));
            f = either_identity ? f : new_f;
        }
        this->index += 1;
        return f;
    }

    Fp12 square_output(Fp12 &f) override {
        return f.square();
    }

    Fp12 conjugate(Fp12 &f) override {
        return f.conjugate();
    }

    Fp12 one() override {
        return Fp12::one();
    }
};

MillerLoopResult multi_miller_loop(const std::vector<std::tuple<G1Affine, G2Prepared>> &terms) {
    HelperMulti helper{terms, 0};
    Fp12 const temp = miller_loop(helper);
    return MillerLoopResult{temp};
}

struct HelperPairing : MillerLoopDriver<Fp12> {
private:
    G2Projective current;
    G2Affine base;
    G1Affine p;
public:
    HelperPairing(G2Projective current, G2Affine base, G1Affine p)
            : current{std::move(current)}, base{std::move(base)}, p{std::move(p)} {}

    Fp12 doubling_step(Fp12 &f) override {
        auto coeffs = pairing::doubling_step(this->current);
        return ell(f, coeffs, this->p);
    }

    Fp12 addition_step(Fp12 &f) override {
        auto coeffs = pairing::addition_step(this->current, this->base);
        return ell(f, coeffs, this->p);
    }

    Fp12 square_output(Fp12 &f) override {
        return f.square();
    }

    Fp12 conjugate(Fp12 &f) override {
        return f.conjugate();
    }

    Fp12 one() override {
        return Fp12::one();
    }
};

Gt pairings(const G1Affine &p, const G2Affine &q) {
    bool const either_identity = p.is_identity() | q.is_identity();
    G1Affine const temp_p = either_identity ? G1Affine::generator() : p;
    G2Affine const temp_q = either_identity ? G2Affine::generator() : q;

    HelperPairing helper{G2Projective{temp_q}, temp_q, temp_p};
    Fp12 const temp = miller_loop(helper);
    MillerLoopResult temp2{either_identity ? Fp12::one() : temp};
    return temp2.final_exponentiation();
}

} // namespace bls12_381::pairing