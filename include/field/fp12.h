#ifndef BLS12_381_FP12_H
#define BLS12_381_FP12_H

#include <optional>

#include "core/rng.h"

#include "field/fp.h"
#include "field/fp2.h"
#include "field/fp6.h"

namespace bls12_381::field {

/**
 * @brief Represents an element of the field Fp12.
 * @details This represents an element c0 + c1 * w of Fp12 = Fp6 / (w^2 - v), where w is the cubic non-residue.
 */
class Fp12 {
private:
    Fp6 c0;
    Fp6 c1;

public:
    Fp12() noexcept;

    Fp12(const Fp12 &fp) noexcept;
    explicit Fp12(const Fp &fp) noexcept;
    explicit Fp12(const Fp2 &fp) noexcept;
    explicit Fp12(const Fp6 &fp) noexcept;
    explicit Fp12(const Fp6 &fp0, const Fp6 &fp1) noexcept;

    Fp12(Fp12 &&fp) noexcept;
    explicit Fp12(Fp &&fp) noexcept;
    explicit Fp12(Fp2 &&fp) noexcept;
    explicit Fp12(Fp6 &&fp) noexcept;
    explicit Fp12(Fp6 &&fp0, Fp6 &&fp1) noexcept;

    ~Fp12() noexcept;

    /**
     * @brief Returns zero, the additive identity.
     * @return The additive identity.
     */
    static Fp12 zero() noexcept;

    /**
     * @brief Returns one, the multiplicative identity.
     * @return The multiplicative identity.
     */
    static Fp12 one() noexcept;
    static Fp12 random(rng::core::RngCore &rng);

    [[nodiscard]] const Fp6 &get_c0() const noexcept;
    [[nodiscard]] const Fp6 &get_c1() const noexcept;

    [[nodiscard]] bool is_zero() const;

    [[nodiscard]] Fp12 square() const;
    [[nodiscard]] Fp12 conjugate() const;
    [[nodiscard]] Fp12 frobenius_map() const;
    [[nodiscard]] Fp12 mul_by_fp2(const Fp2 &fp0, const Fp2 &fp1, const Fp2 &fp4) const;

    /**
     * @brief Computes the multiplicative inverse of this <tt>Fp12</tt> element.
     * @return The multiplicative inverse of this element, if exists.
     */
    [[nodiscard]] auto invert() const -> std::optional<Fp12>;

public:
    Fp12 operator-() const;
    Fp12 &operator=(const Fp12 &rhs);
    Fp12 &operator=(Fp12 &&rhs) noexcept;

    Fp12 &operator+=(const Fp12 &rhs);
    Fp12 &operator-=(const Fp12 &rhs);
    Fp12 &operator*=(const Fp12 &rhs);

    Fp12 operator+(const Fp12 &rhs) const;
    Fp12 operator-(const Fp12 &rhs) const;
    Fp12 operator*(const Fp12 &rhs) const;

private:
    friend inline bool operator==(const Fp12 &a, const Fp12 &b) { return a.c0 == b.c0 && a.c1 == b.c1; }
    friend inline bool operator!=(const Fp12 &a, const Fp12 &b) { return a.c0 != b.c0 || a.c1 != b.c1; }
};
} // namespace bls12_381::field

#endif //BLS12_381_FP12_H