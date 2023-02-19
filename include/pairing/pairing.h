#ifndef BLS12_381_PAIRING_H
#define BLS12_381_PAIRING_H

#include <tuple>
#include <vector>

#include "group/g1_affine.h"
#include "group/g2_prepared.h"
#include "pairing/miller_loop_result.h"
#include "pairing/miller_loop_driver.h"

namespace bls12_381::pairing {

template<typename T>
T miller_loop(MillerLoopDriver<T> &driver) {
    T f = driver.one();

    bool found_one = false;
    for (int i = 63; i >= 0; --i) {
        bool bit = (((group::constant::BLS_X >> 1) >> i) & 1) == 1;
        if (!found_one) {
            found_one = bit;
            continue;
        }
        f = driver.doubling_step(f);
        if (bit) f = driver.addition_step(f);
        f = driver.square_output(f);
    }
    f = driver.doubling_step(f);
    f = driver.conjugate(f);
    return f;
}

template<>
void miller_loop(MillerLoopDriver<void> &driver);

MillerLoopResult multi_miller_loop(const std::vector<std::tuple<group::G1Affine, group::G2Prepared>> &terms);
group::Gt pairings(const group::G1Affine &p, const group::G2Affine &q);
field::Fp12 ell(const field::Fp12 &f,
                const std::tuple<field::Fp2, field::Fp2, field::Fp2> &coefficients,
                const group::G1Affine &p);

} // namespace bls12_381::pairing

#endif //BLS12_381_PAIRING_H