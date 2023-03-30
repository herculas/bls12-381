#ifndef BLS12_381_SCALAR_H
#define BLS12_381_SCALAR_H

#include <array>
#include <cstdint>
#include <optional>

#include "core/rng.h"

#include "group/g1_affine.h"
#include "group/g1_projective.h"
#include "group/g2_affine.h"
#include "group/g2_projective.h"

namespace bls12_381::scalar {

/**
 * @brief Represents an element of the scalar field Fq of the BLS12-381 elliptic curve.
 * @details The internal representation of the scalar is 4 64-bit unsigned integer in little-endian order.
 * @note <tt>Scalar</tt> values are always in Montgomery form, i.e., Scalar(a) = a * R mod q, where R = 2 ^ 256.
 */
class Scalar {
public:
    static constexpr size_t WIDTH = 4;
    static constexpr size_t BYTE_SIZE = WIDTH * sizeof(uint64_t);

private:
    std::array<uint64_t, Scalar::WIDTH> data;

public:
    Scalar() noexcept;

    Scalar(const Scalar &scalar) noexcept;
    explicit Scalar(uint64_t val) noexcept;
    explicit Scalar(const std::array<uint64_t, Scalar::WIDTH> &data) noexcept;

    Scalar(Scalar &&scalar) noexcept;
    explicit Scalar(std::array<uint64_t, Scalar::WIDTH> &&data) noexcept;

    ~Scalar() noexcept;

    /**
     * @brief Returns zero, the additive identity.
     * @return The additive identity.
     */
    static Scalar zero() noexcept;

    /**
     * @brief Returns one, the multiplicative identity.
     * @return The multiplicative identity.
     */
    static Scalar one() noexcept;
    static Scalar random(rng::core::RngCore &rng);

    static Scalar montgomery_reduce(const std::array<uint64_t, Scalar::WIDTH * 2> &rs);

    /**
     * @brief Converts an little-endian integer into its (congruent) <tt>Scalar</tt> representation.
     * @param values the little-endian integer to be converted.
     * @return the <tt>Scalar</tt> element.
     */
    static Scalar from_raw(const std::array<uint64_t, Scalar::WIDTH> &values);

    /**
     * @brief Converts a 512-bit little-endian integer into a <tt>Scalar</tt> by reducing by the modulus.
     * @param bytes the 512-bit little-endian integer to be converted.
     * @return the <tt>Scalar</tt> element.
     */
    static Scalar from_bytes_wide(const std::array<uint8_t, Scalar::BYTE_SIZE * 2> &bytes);

    /**
     * @brief Attempts to converts a byte array into a <tt>Scalar</tt> element.
     * @param bytes the byte array to be converted.
     * @return the <tt>Scalar</tt> element if the input is canonical.
     */
    static std::optional<Scalar> from_bytes(const std::array<uint8_t, Scalar::BYTE_SIZE> &bytes);

    [[nodiscard]] bool is_zero() const;
    [[nodiscard]] std::string to_hex_str() const;

    /**
     * @brief Converts an <tt>Scalar</tt> element into a byte array in little-endian order.
     * @return the byte representation in little-endian order.
     */
    [[nodiscard]] std::array<uint8_t, Scalar::BYTE_SIZE> to_bytes() const;

    /**
     * @brief Doubles the <tt>Scalar</tt> element.
     * @return the doubled element.
     */
    [[nodiscard]] Scalar doubles() const;

    /**
     * @brief Squares the <tt>Scalar</tt> element.
     * @return the squared element.
     */
    [[nodiscard]] Scalar square() const;
    [[nodiscard]] Scalar subtract_modulus() const;

    /**
     * @brief Computes the power of the <tt>Scalar</tt> element.
     * @param exp the exponent, a little-endian integer.
     * @return the power of this element.
     */
    [[nodiscard]] Scalar pow(const std::array<uint64_t, Scalar::WIDTH> &exp) const;

    /**
     * @brief Attempts to compute the square root of the <tt>Scalar</tt> element.
     * @return the square root of this element, if it exists.
     */
    [[nodiscard]] std::optional<Scalar> sqrt() const;

    /**
     * @brief Attempts to compute the inverse of the <tt>Scalar</tt> element.
     * @return the inverse of this element, if it exists.
     */
    [[nodiscard]] std::optional<Scalar> invert() const;

private:
    /**
     * @brief Reduces a 512-bit little-endian integer by decomposing it into two 256-bit digits with the higher bits
     *          multiplied by 2 ^ 256.
     * @param limbs the 512-bit little-endian integer to be reduced.
     * @return the reduced <tt>Scalar</tt> element.
     */
    static Scalar reduce(const std::array<uint64_t, Scalar::WIDTH * 2> &limbs);

public:
    Scalar operator-() const;
    Scalar &operator=(const Scalar &rhs);
    Scalar &operator=(Scalar &&rhs) noexcept;

    Scalar &operator+=(const Scalar &rhs);
    Scalar &operator-=(const Scalar &rhs);
    Scalar &operator*=(const Scalar &rhs);

public:
    friend inline Scalar operator+(const Scalar &a, const Scalar &b) { return Scalar(a) += b; }
    friend inline Scalar operator-(const Scalar &a, const Scalar &b) { return Scalar(a) -= b; }
    friend inline Scalar operator*(const Scalar &a, const Scalar &b) { return Scalar(a) *= b; }

    friend group::G1Projective operator*(const Scalar &a, const group::G1Affine &b);
    friend group::G1Projective operator*(const Scalar &a, const group::G1Projective &b);
    friend group::G2Projective operator*(const Scalar &a, const group::G2Affine &b);
    friend group::G2Projective operator*(const Scalar &a, const group::G2Projective &b);

    friend inline bool operator==(const Scalar &a, const Scalar &b) { return a.data == b.data; }
    friend inline bool operator!=(const Scalar &a, const Scalar &b) { return a.data != b.data; }
};

} // namespace bls12_381::scalar

#endif //BLS12_381_SCALAR_H
