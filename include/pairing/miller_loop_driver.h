#ifndef BLS12_381_MILLER_LOOP_DRIVER_H
#define BLS12_381_MILLER_LOOP_DRIVER_H

#include <array>

#include "field/fp2.h"
#include "group/constant.h"
#include "group/g2_affine.h"
#include "group/g2_projective.h"

namespace bls12_381::pairing {

template<typename T>
struct MillerLoopDriver {
    virtual T doubling_step(T f) = 0;
    virtual T addition_step(T f) = 0;
    virtual T square_output(T f) = 0;
    virtual T conjugate(T f) = 0;
    virtual T one() = 0;
};

template<>
struct MillerLoopDriver<void> {
    virtual void doubling_step() = 0;
    virtual void addition_step() = 0;
    virtual void square_output() = 0;
    virtual void conjugate() = 0;
    virtual void one() = 0;
};

template<typename T>
T miller_loop(MillerLoopDriver<T> &driver) {
    T f = MillerLoopDriver<T>::one();

    bool found_one = false;
    for (int i = 63; i >= 0; --i) {
        bool bit = (((group::constant::BLS_X >> 1) >> i) & 1) == 1;
        if (!found_one) {
            found_one = bit;
            continue;
        }
        f = driver.doubling_step(f);
        if (bit) f = driver.addition_step(f);
        f = MillerLoopDriver<T>::square_output(f);
    }
    f = driver.doubling_step(f);
    f = MillerLoopDriver<T>::conjugate(f);
    return f;
}

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
    }
    driver.doubling_step();
}

std::array<field::Fp2, 3> doubling_step(group::G2Projective &point);
std::array<field::Fp2, 3> addition_step(group::G2Projective &r, const group::G2Affine &q);

} // namespace bls12_381::pairing

#endif //BLS12_381_MILLER_LOOP_DRIVER_H