#ifndef BLS12_381_G2_PROJECTIVE_H
#define BLS12_381_G2_PROJECTIVE_H

#include <array>
#include <cstdint>
#include <vector>

#include "core/rng.h"

#include "field/fp2.h"

namespace bls12_381::scalar { class Scalar; }
namespace bls12_381::group { class G2Affine; }

namespace bls12_381::group {

/**
 * @brief An element on the G2 curve in projective coordinate space.
 */
class G2Projective {
private:
    field::Fp2 x;
    field::Fp2 y;
    field::Fp2 z;

public:
    G2Projective();

    G2Projective(const G2Projective &point);
    explicit G2Projective(const G2Affine &point);
    explicit G2Projective(const field::Fp2 &x, const field::Fp2 &y, const field::Fp2 &z);

    G2Projective(G2Projective &&point) noexcept;
    explicit G2Projective(G2Affine &&point);
    explicit G2Projective(field::Fp2 &&x, field::Fp2 &&y, field::Fp2 &&z);

    ~G2Projective();

    /**
     * @brief Returns the identity element of G2 in projective coordinate form.
     * @return The identity element of G2.
     * @note The identity element is the point at infinity.
     */
    static G2Projective identity() noexcept;

    /**
     * @brief Returns a fixed generator of G2 in projective coordinate form.
     * @return A chosen generator of G2.
     */
    static G2Projective generator() noexcept;
    static G2Projective random(rng::core::RngCore &rng);

    /**
     * @brief Converts a batch of <tt>G2Projective</tt> elements into <tt>G2Affine</tt> elements.
     * @param points The vector of <tt>G2Projective</tt> elements to be converted.
     * @return The vector of <tt>G2Affine</tt> elements.
     */
    static std::vector<G2Affine> batch_normalize(const std::vector<G2Projective> &points);

    [[nodiscard]] const field::Fp2 &get_x() const noexcept;
    [[nodiscard]] const field::Fp2 &get_y() const noexcept;
    [[nodiscard]] const field::Fp2 &get_z() const noexcept;

    /**
     * @brief Checks if the point is the identity element, i.e. the point at infinity.
     * @return <tt>true</tt> if the point is the identity element, <tt>false</tt> otherwise.
     */
    [[nodiscard]] bool is_identity() const;

    /**
     * @brief Checks if the point is on the curve.
     * @return <tt>true</tt> if the point is on the curve, <tt>false</tt> otherwise.
     * @note This should always return <tt>true</tt> unless an unchecked API was misused.
     */
    [[nodiscard]] bool is_on_curve() const;

    /**
     * @brief Doubles this point.
     * @return The doubled point.
     */
    [[nodiscard]] G2Projective doubles() const;

    /**
     * @brief Multiplies this point by <tt>BLS_X</tt>.
     * @return The multiplication result.
     */
    [[nodiscard]] G2Projective mul_by_x() const;

    /**
     * @brief Clear the cofactor using Budroni-Pintore algorithm.
     * @details This is equivalent to multiplying by h_eff = 3(z^2 - 1) * h_2, where h_2 is the cofactor of G2 and z
     *          is the parameter of BLS12-381.
     * @return The point multiplied by the cofactor.
     */
    [[nodiscard]] G2Projective clear_cofactor() const;
    [[nodiscard]] G2Projective psi() const;
    [[nodiscard]] G2Projective psi2() const;

private:
    /**
     * @brief Adds this point to another point.
     * @param rhs The point to add to this point.
     * @return The sum of the two points.
     */
    [[nodiscard]] G2Projective add(const G2Projective &rhs) const;

    /**
     * @brief Adds this point to another point in affine coordinate space.
     * @param rhs The point to be added to this point in affine coordinate space.
     * @return The sum of two G1 points in projective coordinate space.
     */
    [[nodiscard]] G2Projective add_mixed(const G2Affine &rhs) const;

    /**
     * @brief Multiplies this point by a byte array.
     * @param bytes The byte array to be multiplied to this point.
     * @return The multiplication result.
     */
    [[nodiscard]] G2Projective multiply(const std::array<uint8_t, 32> &bytes) const;

public:
    G2Projective operator-() const;
    G2Projective &operator=(const G2Projective &rhs);
    G2Projective &operator=(G2Projective &&rhs) noexcept;

    G2Projective &operator+=(const G2Projective &rhs);
    G2Projective &operator-=(const G2Projective &rhs);

    G2Projective &operator+=(const G2Affine &rhs);
    G2Projective &operator-=(const G2Affine &rhs);

    G2Projective &operator*=(const scalar::Scalar &rhs);

public:
    friend inline G2Projective operator+(const G2Projective &a, const G2Affine &b) { return G2Projective(a) += b; }
    friend inline G2Projective operator-(const G2Projective &a, const G2Affine &b) { return G2Projective(a) -= b; }

    friend inline G2Projective operator+(const G2Projective &a, const G2Projective &b) { return G2Projective(a) += b; }
    friend inline G2Projective operator-(const G2Projective &a, const G2Projective &b) { return G2Projective(a) -= b; }

    friend inline G2Projective operator*(const G2Projective &a, const scalar::Scalar &b) { return G2Projective(a) *= b; }

    friend inline bool operator==(const G2Projective &a, const G2Projective &b) {
        field::Fp2 const x1 = a.x * b.z;
        field::Fp2 const x2 = b.x * a.z;
        field::Fp2 const y1 = a.y * b.z;
        field::Fp2 const y2 = b.y * a.z;

        bool const a_is_zero = a.z.is_zero();
        bool const b_is_zero = b.z.is_zero();

        return (a_is_zero & b_is_zero) | ((!a_is_zero) & (!b_is_zero) & (x1 == x2) & (y1 == y2));
    }

    friend inline bool operator!=(const G2Projective &a, const G2Projective &b) {
        return !(a == b);
    }
};

} // namespace bls12_381::group

#endif //BLS12_381_G2_PROJECTIVE_H