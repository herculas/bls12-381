#ifndef BLS12_381_MILLER_LOOP_H
#define BLS12_381_MILLER_LOOP_H

#include "field/fp12.h"
#include "group/gt.h"

namespace bls12_381::pairing {

class MillerLoopResult {
private:
    field::Fp12 data;

public:
    MillerLoopResult();
    explicit MillerLoopResult(const field::Fp12 &data);
    explicit MillerLoopResult(field::Fp12 &&data);

    group::Gt final_exponentiation();

public:
    MillerLoopResult &operator+=(const MillerLoopResult &rhs);

public:
    friend inline MillerLoopResult operator+(const MillerLoopResult &a, const MillerLoopResult &b) {
        return MillerLoopResult(a) += b;
    }
};

} // namespace bls12_381::pairing

#endif //BLS12_381_MILLER_LOOP_H