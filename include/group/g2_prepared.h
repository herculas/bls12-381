#ifndef BLS12_381_G2_PREPARED_H
#define BLS12_381_G2_PREPARED_H

#include <tuple>
#include <vector>

#include "field/fp2.h"
#include "group/g2_affine.h"
#include "pairing/miller_loop_driver.h"

namespace bls12_381::group {

/**
 * @brief This represents caches computations pertaining to a G2 element as part of the pairing function (specifically,
 *          the Miller loop)
 * @note This should be computed whenever a G2 element is being used in multiple pairings or is otherwise known in
 *          advance. This should be used in conjunction with the <tt>multi_miller_loop</tt> function provided in the
 *          <tt>pairing</tt> namespace.
 */
class G2Prepared {
private:
    bool infinity;
    std::vector<std::tuple<field::Fp2, field::Fp2, field::Fp2>> coefficients;

public:
    G2Prepared(bool infinity, std::vector<std::tuple<field::Fp2, field::Fp2, field::Fp2>> coefficients);
    explicit G2Prepared(const G2Affine &point);
    explicit G2Prepared(G2Affine &&point);

    [[nodiscard]] bool is_identity() const;
    [[nodiscard]] const std::vector<std::tuple<field::Fp2, field::Fp2, field::Fp2>> &get_coeffs() const;

    /**
     * @brief Creates a <tt>G2Prepared</tt> element from a set of bytes created by <tt>G2Prepared::to_raw_bytes</tt>.
     * @param bytes A vector of bytes representing a <tt>G2Prepared</tt> element.
     * @return A <tt>G2Prepared</tt> element, if exists.
     * @note No check is performed and no constant time is guaranteed. The <tt>infinity</tt> attribute is also lost.
     *          The expected usage of this function is for trusted sets of data where performance is critical.
     */
    static auto from_slice_unchecked(const std::vector<uint8_t> &bytes) -> G2Prepared;

    /**
     * @brief Converts a <tt>G2Prepared</tt> element into raw bytes representation.
     * @return A vector of bytes representing the <tt>G2Prepared</tt> element.
     * @note The intended usage of this function is for trusted sets of data where performance is critical. This way,
     *          the <tt>infinity</tt> attribute will not be serialized, and the <tt>coefficients</tt> attribute will be
     *          stored without any check.
     */
    [[nodiscard]] std::vector<uint8_t> to_raw_bytes() const;
};

} // namespace bls12_381::group

#endif //BLS12_381_G2_PREPARED_H