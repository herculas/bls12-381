#ifndef BLS12_381_G1_PROJECTIVE_H
#define BLS12_381_G1_PROJECTIVE_H

#include <array>
#include <cstdint>
#include <vector>

#include "core/rng.h"

#include "field/fp.h"

namespace bls12_381::scalar { class Scalar; }
namespace bls12_381::group { class G1Affine; }

namespace bls12_381::group {

/**
 * @brief An element on the G1 curve in projective coordinate space.
 */
class G1Projective {
private:
    field::Fp x;
    field::Fp y;
    field::Fp z;

public:
    G1Projective();

    G1Projective(const G1Projective &point);
    explicit G1Projective(const G1Affine &point);
    explicit G1Projective(const field::Fp &x, const field::Fp &y, const field::Fp &z);

    G1Projective(G1Projective &&point) noexcept;
    explicit G1Projective(G1Affine &&point);
    explicit G1Projective(field::Fp &&x, field::Fp &&y, field::Fp &&z);

    ~G1Projective();

    /**
     * @brief Returns the identity element of G1 in projective coordinate form.
     * @return The identity element of G1.
     * @note The identity element is the point at infinity.
     */
    static G1Projective identity() noexcept;

    /**
     * @brief Returns a fixed generator of G1 in projective coordinate form.
     * @return A chosen generator of G1.
     */
    static G1Projective generator() noexcept;
    static G1Projective random(rng::core::RngCore &rng);

    /**
     * @brief Converts a batch of <tt>G1Projective</tt> elements into <tt>G1Affine</tt> elements.
     * @param points The vector of <tt>G1Projective</tt> elements to be converted.
     * @return The vector of <tt>G1Affine</tt> elements.
     */
    static std::vector<G1Affine> batch_normalize(const std::vector<G1Projective> &points);

    [[nodiscard]] const field::Fp &get_x() const noexcept;
    [[nodiscard]] const field::Fp &get_y() const noexcept;
    [[nodiscard]] const field::Fp &get_z() const noexcept;

    /**
     * @brief Checks if the point is the identity element, i.e., the point at infinity.
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
    [[nodiscard]] G1Projective doubles() const;

    /**
     * @brief Multiplies this point by <tt>BLS_X</tt>.
     * @return The multiplication result.
     */
    [[nodiscard]] G1Projective mul_by_x() const;

    /**
     * @brief Clears the cofactor of this point.
     * @details Multiplies this point by (1 - z) where z is the parameter of BLS12-381, which suffices to clear the
     *          cofactor and map elliptic curve points to elements of G1.
     * @return The multiplication result.
     */
    [[nodiscard]] G1Projective clear_cofactor() const;

private:
    /**
     * @brief Adds this point to another point.
     * @param rhs The point to be added to this point.
     * @return The sum of two G1 point.
     */
    [[nodiscard]] G1Projective add(const G1Projective &rhs) const;

    /**
     * @brief Adds this point to another point in affine coordinate space.
     * @param rhs The point to be added to this point in affine coordinate space.
     * @return The sum of two G1 points in projective coordinate space.
     */
    [[nodiscard]] G1Projective add_mixed(const G1Affine &rhs) const;

    /**
     * @brief Multiplies this point by a byte array.
     * @param bytes The byte array to be multiplied to this point.
     * @return The multiplication result.
     */
    [[nodiscard]] G1Projective multiply(const std::array<uint8_t, 32> &bytes) const;

public:
    G1Projective operator-() const;
    G1Projective &operator=(const G1Projective &rhs);
    G1Projective &operator=(G1Projective &&rhs) noexcept;

    G1Projective &operator+=(const G1Projective &rhs);
    G1Projective &operator-=(const G1Projective &rhs);
    G1Projective &operator+=(const G1Affine &rhs);
    G1Projective &operator-=(const G1Affine &rhs);
    G1Projective &operator*=(const scalar::Scalar &rhs);

    G1Projective operator+(const G1Projective &rhs) const;
    G1Projective operator-(const G1Projective &rhs) const;
    G1Projective operator+(const G1Affine &rhs) const;
    G1Projective operator-(const G1Affine &rhs) const;
    G1Projective operator*(const scalar::Scalar &rhs) const;

public:
    friend inline bool operator==(const G1Projective &a, const G1Projective &b) {
        field::Fp const x1 = a.x * b.z;
        field::Fp const x2 = b.x * a.z;
        field::Fp const y1 = a.y * b.z;
        field::Fp const y2 = b.y * a.z;

        bool const a_is_zero = a.z.is_zero();
        bool const b_is_zero = b.z.is_zero();

        return (a_is_zero & b_is_zero) | ((!a_is_zero) & (!b_is_zero) & (x1 == x2) & (y1 == y2));
    }

    friend inline bool operator!=(const G1Projective &a, const G1Projective &b) {
        return !(a == b);
    }
};

} // namespace bls12_381::group

#endif //BLS12_381_G1_PROJECTIVE_H