#ifndef BLS12_381_PAIRING_H
#define BLS12_381_PAIRING_H

#include <tuple>
#include <vector>

#include "group/g1_affine.h"
#include "group/g2_prepared.h"
#include "pairing/miller_loop_result.h"
#include "pairing/miller_loop_driver.h"

namespace bls12_381::pairing {

/**
 * @brief A generic implementation of the Miller loop.
 * @tparam T The type of the result of the Miller loop.
 * @param driver The driver of the Miller loop.
 * @return The result of the Miller loop of type <tt>T</tt>.
 */
template<typename T>
T miller_loop(MillerLoopDriver<T> &driver) {
    T f = driver.one();

    bool found_one = false;
    for (int i = 63; i >= 0; --i) { // NOLINT(cppcoreguidelines-avoid-magic-numbers)
        bool const bit = (((group::constant::BLS_X >> 1) >> i) & 1) == 1;
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

/**
 * @brief Computes the sum of miller_loop(a_i, b_i) given a series of terms (a_1, b_1), ..., (a_n, b_n).
 * @param terms A series of terms (a_i, b_i) where a's are <tt>G1Affine</tt> elements and b's are <tt>G2Prepared</tt>
 *          elements.
 * @return The result of the multiple Miller loop.
 */
MillerLoopResult multi_miller_loop(const std::vector<std::tuple<group::G1Affine, group::G2Prepared>> &terms);

/**
 * @brief Invokes the pairing function without the use of pre-computation and other optimizations.
 * @param p A <tt>G1Affine</tt> element.
 * @param q A <tt>G2Affine</tt> element.
 * @return The result of the pairing function, a <tt>Gt</tt> element.
 */
group::Gt pairings(const group::G1Affine &p, const group::G2Affine &q);

} // namespace bls12_381::pairing

#endif //BLS12_381_PAIRING_H