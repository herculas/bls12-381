#ifndef BLS12_381_G2_PREPARED_H
#define BLS12_381_G2_PREPARED_H

#include <tuple>
#include <vector>

#include "field/fp2.h"
#include "group/g2_affine.h"
#include "pairing/miller_loop_driver.h"

namespace bls12_381::group {

class G2Prepared {
private:
    bool infinity;
    std::vector<std::tuple<field::Fp2, field::Fp2, field::Fp2>> coefficients;

public:
    G2Prepared(bool infinity, std::vector<std::tuple<field::Fp2, field::Fp2, field::Fp2>> coefficients);
    explicit G2Prepared(const G2Affine &point);
    explicit G2Prepared(G2Affine &&point);

    [[nodiscard]] bool is_identity() const;
    [[nodiscard]] std::vector<std::tuple<field::Fp2, field::Fp2, field::Fp2>> get_coeffs() const;
};

} // namespace bls12_381::group

#endif //BLS12_381_G2_PREPARED_H