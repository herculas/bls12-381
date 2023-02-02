#include "pairing/miller_loop_result.h"

#include <tuple>
#include "field/fp2.h"
#include "group/constant.h"

namespace bls12_381::pairing {

MillerLoopResult::MillerLoopResult() : data{field::Fp12::one()} {}

MillerLoopResult::MillerLoopResult(const field::Fp12 &data) : data{data} {}

MillerLoopResult::MillerLoopResult(field::Fp12 &&data) : data{data} {}

std::tuple<field::Fp2, field::Fp2> fp4_square(const field::Fp2 &a, const field::Fp2 &b) {
    field::Fp2 t0 = a.square();
    field::Fp2 t1 = b.square();
    field::Fp2 t2 = t1.mul_by_non_residue();
    field::Fp2 c0 = t2 + t0;
    t2 = a + b;
    t2 = t2.square();
    t2 -= t0;
    field::Fp2 c1 = t2 - t1;
    return {c0, c1};
}

field::Fp12 cyclotomic_square(const field::Fp12 &point) {
    field::Fp2 z0 = point.get_c0().get_c0();
    field::Fp2 z4 = point.get_c0().get_c1();
    field::Fp2 z3 = point.get_c0().get_c2();
    field::Fp2 z2 = point.get_c1().get_c0();
    field::Fp2 z1 = point.get_c1().get_c1();
    field::Fp2 z5 = point.get_c1().get_c2();

    auto [t0, t1] = fp4_square(z0, z1);

    // A
    z0 = t0 - z0;
    z0 = z0 + z0 + t0;
    z1 = t1 + z1;
    z1 = z1 + z1 + t1;

    std::tie(t0, t1) = fp4_square(z2, z3);
    auto [t2, t3] = fp4_square(z4, z5);

    // C
    z4 = t0 - z4;
    z4 = z4 + z4 + t0;
    z5 = t1 + z5;
    z5 = z5 + z5 + t1;

    // B
    t0 = t3.mul_by_non_residue();
    z2 = t0 + z2;
    z2 = z2 + z2 + t0;
    z3 = t2 - z3;
    z3 = z3 + z3 + t2;

    return field::Fp12{
            field::Fp6{z0, z4, z3},
            field::Fp6{z2, z1, z5},
    };
}

field::Fp12 cyclotomic_exp(const field::Fp12 &point) {
    uint64_t x = group::constant::BLS_X;
    field::Fp12 temp = field::Fp12::one();
    bool found_one = false;
    for (int i = 63; i >= 0; --i) {
        bool bit = ((x >> i) & 1) == 1;
        if (found_one) temp = cyclotomic_square(temp);
        else found_one = bit;
        if (bit) temp *= point;
    }
    return temp.conjugate();
}

group::Gt MillerLoopResult::final_exponentiation() {
    field::Fp12 f = this->data;
    field::Fp12 t0 = f.frobenius_map().frobenius_map().frobenius_map().frobenius_map().frobenius_map().frobenius_map();

    assert(f.invert().has_value());

    field::Fp12 t1 = f.invert().value();
    field::Fp12 t2 = t0 * t1;
    t1 = t2;
    t2 = t2.frobenius_map().frobenius_map();
    t2 *= t1;
    t1 = cyclotomic_square(t2).conjugate();
    field::Fp12 t3 = cyclotomic_exp(t2);
    field::Fp12 t4 = cyclotomic_square(t3);
    field::Fp12 t5 = t1 * t3;
    t1 = cyclotomic_exp(t5);
    t0 = cyclotomic_exp(t1);
    field::Fp12 t6 = cyclotomic_exp(t0);
    t6 *= t4;
    t4 = cyclotomic_exp(t6);
    t5 = t5.conjugate();
    t4 *= t5 * t2;
    t5 = t2.conjugate();
    t1 *= t2;
    t1 = t1.frobenius_map().frobenius_map().frobenius_map();
    t6 *= t5;
    t6 = t6.frobenius_map();
    t3 *= t0;
    t3 = t3.frobenius_map().frobenius_map();
    t3 *= t1;
    t3 *= t6;
    f = t3 * t4;

    return group::Gt{f};
}

MillerLoopResult &MillerLoopResult::operator+=(const MillerLoopResult &rhs) {
    *this = MillerLoopResult{this->data * rhs.data};
    return *this;
}

field::Fp12 MillerLoopResult::get_data() const {
    return this->data;
}

} // namespace bls12_381::pairing