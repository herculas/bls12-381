#ifndef BLS12_381_MILLER_LOOP_DRIVER_H
#define BLS12_381_MILLER_LOOP_DRIVER_H

#include <tuple>

#include "field/fp2.h"
#include "field/fp12.h"
#include "group/constant.h"
#include "group/g2_affine.h"
#include "group/g2_projective.h"

namespace bls12_381::pairing {

template<typename T>
struct MillerLoopDriver {
    virtual T doubling_step(T &f) = 0;
    virtual T addition_step(T &f) = 0;
    virtual T square_output(T &f) = 0;
    virtual T conjugate(T &f) = 0;
    virtual T one() = 0;
    virtual ~MillerLoopDriver() = default;
};

template<>
struct MillerLoopDriver<void> {
    virtual void doubling_step() = 0;
    virtual void addition_step() = 0;
    virtual void square_output() = 0;
    virtual void conjugate() = 0;
    virtual void one() = 0;
    virtual ~MillerLoopDriver() = default;
};

std::tuple<field::Fp2, field::Fp2, field::Fp2> doubling_step(group::G2Projective &point);
std::tuple<field::Fp2, field::Fp2, field::Fp2> addition_step(group::G2Projective &r, const group::G2Affine &q);

} // namespace bls12_381::pairing

#endif //BLS12_381_MILLER_LOOP_DRIVER_H