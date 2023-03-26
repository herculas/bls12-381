#ifndef BLS12_381_FP6_H
#define BLS12_381_FP6_H

#include <optional>

#include "core/rng.h"

#include "field/fp.h"
#include "field/fp2.h"

namespace bls12_381::field {

/**
 * @brief Represents an element of the field Fp6.
 * @details This represents an element c0 + c1 * v + c2 * v^2 of Fp6 = Fp2 / (v^3 - u - 1), where v is the quadratic
 *          non-residue.
 */
class Fp6 {
private:
    Fp2 c0;
    Fp2 c1;
    Fp2 c2;

public:
    Fp6();

    Fp6(const Fp6 &fp);
    explicit Fp6(const Fp &fp);
    explicit Fp6(const Fp2 &fp);
    explicit Fp6(const Fp2 &fp0, const Fp2 &fp1, const Fp2 &fp2);

    Fp6(Fp6 &&fp) noexcept;
    explicit Fp6(Fp &&fp);
    explicit Fp6(Fp2 &&fp);
    explicit Fp6(Fp2 &&fp0, Fp2 &&fp1, Fp2 &&fp2);

    /**
     * @brief Returns zero, the additive identity.
     * @return The additive identity.
     */
    static Fp6 zero() noexcept;

    /**
     * @brief Returns one, the multiplicative identity.
     * @return The multiplicative identity.
     */
    static Fp6 one() noexcept;
    static Fp6 random(rng::core::RngCore &rng);

    [[nodiscard]] const Fp2 &get_c0() const noexcept;
    [[nodiscard]] const Fp2 &get_c1() const noexcept;
    [[nodiscard]] const Fp2 &get_c2() const noexcept;

    [[nodiscard]] bool is_zero() const;

    [[nodiscard]] Fp6 square() const;

    /**
     * @brief Raise this element to p.
     * @return The Frobenius map of this element.
     */
    [[nodiscard]] Fp6 frobenius_map() const;

    /**
     * @brief multiplies this element by the quadratic non-residue v.
     * @return The result of the multiplication.
     */
    [[nodiscard]] Fp6 mul_by_non_residue() const;
    [[nodiscard]] Fp6 mul_by_fp2(const Fp2 &fp) const;
    [[nodiscard]] Fp6 mul_by_fp2(const Fp2 &fp0, const Fp2 &fp1) const;

    /**
     * @brief Computes the multiplicative inverse of this <tt>Fp6</tt> element.
     * @return The multiplicative inverse of this element, if exists.
     */
    [[nodiscard]] auto invert() const -> std::optional<Fp6>;

private:
    /**
     * @brief Multiplies this element by another element using the Longa's algorithm.
     * @param b The other element.
     * @return The result of the multiplication.
     */
    [[nodiscard]] Fp6 mul_interleaved(const Fp6 &b) const;

public:
    Fp6 operator-() const;
    Fp6 &operator=(const Fp6 &rhs);
    Fp6 &operator=(Fp6 &&rhs) noexcept;

    Fp6 &operator+=(const Fp6 &rhs);
    Fp6 &operator-=(const Fp6 &rhs);
    Fp6 &operator*=(const Fp6 &rhs);

    Fp6 operator+(const Fp6 &rhs) const;
    Fp6 operator-(const Fp6 &rhs) const;
    Fp6 operator*(const Fp6 &rhs) const;

public:
    friend inline bool operator==(const Fp6 &a, const Fp6 &b) { return a.c0 == b.c0 && a.c1 == b.c1 && a.c2 == b.c2; }
    friend inline bool operator!=(const Fp6 &a, const Fp6 &b) { return a.c0 != b.c0 || a.c1 != b.c1 || a.c2 != b.c2; }
};

} // namespace bls12_381::field

#endif //BLS12_381_FP6_H