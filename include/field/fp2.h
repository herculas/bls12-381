#ifndef BLS12_381_FP2_H
#define BLS12_381_FP2_H

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "core/rng.h"

#include "field/fp.h"

namespace bls12_381::field {

/**
 * @brief Represents an element of the field Fp2.
 */
class Fp2 {
private:
    Fp c0;
    Fp c1;

public:
    Fp2() noexcept;

    Fp2(const Fp2 &fp) noexcept;
    explicit Fp2(const Fp &fp) noexcept;
    explicit Fp2(const Fp &fp0, const Fp &fp1) noexcept;

    Fp2(Fp2 &&fp) noexcept;
    explicit Fp2(Fp &&fp) noexcept;
    explicit Fp2(Fp &&fp0, Fp &&fp1) noexcept;

    ~Fp2() noexcept;

    /**
     * @brief Returns zero, the additive identity.
     * @return The additive identity.
     */
    static Fp2 zero() noexcept;

    /**
     * @brief Returns one, the multiplicative identity.
     * @return The multiplicative identity.
     */
    static Fp2 one() noexcept;
    static Fp2 random(rng::core::RngCore &rng);

    [[nodiscard]] const Fp &get_c0() const noexcept;
    [[nodiscard]] const Fp &get_c1() const noexcept;

    [[nodiscard]] bool is_zero() const;

    /**
     * @brief Checks if this <tt>Fp2</tt> element is strictly lexicographically larger than its negation.
     * @return <tt>true</tt> if this element is lexicographically larger than its negation, <tt>false</tt> otherwise.
     */
    [[nodiscard]] bool lexicographically_largest() const;
    [[nodiscard]] std::string to_hex_str() const;

    [[nodiscard]] Fp2 square() const;
    [[nodiscard]] Fp2 conjugate() const;

    /**
     * @brief Raise this element to p.
     * @return The Frobenius map of this element.
     */
    [[nodiscard]] Fp2 frobenius_map() const;
    [[nodiscard]] Fp2 mul_by_non_residue() const;

    /**
     * @brief Computes the power of this <tt>Fp2</tt> element.
     * @param exp The exponent.
     * @return The power of this element to the given exponent.
     */
    [[nodiscard]] Fp2 pow(const std::array<uint64_t, Fp::WIDTH> &exp) const;
    [[nodiscard]] Fp2 pow_extended(const std::vector<uint64_t> &exp) const;

    /**
     * @brief Computes the square root of this <tt>Fp2</tt> element.
     * @return The square root of this element, if exists.
     */
    [[nodiscard]] auto sqrt() const -> std::optional<Fp2>;

    /**
     * @brief Computes the multiplicative inverse of this <tt>Fp2</tt> element.
     * @return The multiplicative inverse of this element, if exists.
     */
    [[nodiscard]] auto invert() const -> std::optional<Fp2>;

public:
    Fp2 operator-() const;
    Fp2 &operator=(const Fp2 &rhs);
    Fp2 &operator=(Fp2 &&rhs) noexcept;

    Fp2 &operator+=(const Fp2 &rhs);
    Fp2 &operator-=(const Fp2 &rhs);
    Fp2 &operator*=(const Fp2 &rhs);

    Fp2 operator+(const Fp2 &rhs) const;
    Fp2 operator-(const Fp2 &rhs) const;
    Fp2 operator*(const Fp2 &rhs) const;

public:
    friend inline bool operator==(const Fp2 &a, const Fp2 &b) { return a.c0 == b.c0 && a.c1 == b.c1; }
    friend inline bool operator!=(const Fp2 &a, const Fp2 &b) { return a.c0 != b.c0 || a.c1 != b.c1; }
};

} // namespace bls12_381::field

#endif //BLS12_381_FP2_H