#ifndef BLS12_381_G2_PREPARED_H
#define BLS12_381_G2_PREPARED_H

#include <array>
#include <vector>

#include "field/fp2.h"
#include "g2_affine.h"

namespace bls12_381::group {

class G2Prepared {
private:
    bool infinity;
    std::vector<std::array<field::Fp2, 3>> coefficients;

public:
    G2Prepared(bool infinity, const std::vector<std::array<field::Fp2, 3>> &coefficients);
    explicit G2Prepared(const G2Affine &point);
    explicit G2Prepared(G2Affine &&point);

};

} // namespace bls12_381::group

#endif //BLS12_381_G2_PREPARED_H