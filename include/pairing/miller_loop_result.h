#ifndef BLS12_381_MILLER_LOOP_RESULT_H
#define BLS12_381_MILLER_LOOP_RESULT_H

#include "field/fp12.h"
#include "group/gt.h"

namespace bls12_381::pairing {

/**
 * @brief Represents results of a Miller loop, one of the most expensive operations in the pairing.
 * @note <tt>MillerLoopResult</tt>s cannot be compared with each other, until <tt>final_exponentiation</tt> is called,
 *          which is also very expensive.
 */
class MillerLoopResult {
private:
    field::Fp12 data;

public:
    MillerLoopResult();
    explicit MillerLoopResult(const field::Fp12 &data);
    explicit MillerLoopResult(field::Fp12 &&data);

    /**
     * @brief Performs a "final exponentiation" routine to convert the result of a Miller loop into a <tt>Gt</tt>
     *          element with help of efficient squaring operation in the so-called <tt>cyclotomic subgroup</tt> of
     *          <tt>Fp6</tt> so that it can be compared with other <tt>Gt</tt> elements.
     * @return The result of the final exponentiation, a <tt>Gt</tt> element.
     */
    group::Gt final_exponentiation();
    [[nodiscard]] const field::Fp12 &get_data() const;

public:
    MillerLoopResult &operator+=(const MillerLoopResult &rhs);

public:
    friend inline MillerLoopResult operator+(const MillerLoopResult &a, const MillerLoopResult &b) {
        return MillerLoopResult(a) += b;
    }
};

} // namespace bls12_381::pairing

#endif //BLS12_381_MILLER_LOOP_RESULT_H