#ifndef BLS12_381_GT_H
#define BLS12_381_GT_H

#include "core/rng.h"

#include "field/fp12.h"

namespace bls12_381::scalar { class Scalar; }

/**
 * @brief This is an element of Gt, the target group of the pairing function. As with G1 and G2, this group has order q.
 * @note Typically, Gt is multiplicatively written, but we write it additively to keep code and abstractions consistent.
 */
namespace bls12_381::group {

class Gt {
private:
    field::Fp12 data;

public:
    Gt();

    Gt(const Gt &point);
    explicit Gt(const field::Fp12 &point);

    Gt(Gt &&point) noexcept;
    explicit Gt(field::Fp12 &&point);

    ~Gt();

    /**
     * @brief Returns the identity element of Gt.
     * @return The identity element of Gt.
     * @note The identity element is 1.
     */
    static Gt identity() noexcept;

    /**
     * @brief Returns a generator of Gt.
     * @return The generator of Gt.
     */
    static Gt generator() noexcept;
    static Gt random(rng::core::RngCore &rng);

    [[nodiscard]] bool is_identity() const;

    /**
     * @brief Doubles this point.
     * @return The doubled point.
     */
    [[nodiscard]] Gt doubles() const;

public:
    Gt operator-() const;
    Gt &operator=(const Gt &rhs);
    Gt &operator=(Gt &&rhs) noexcept;

    Gt &operator+=(const Gt &rhs);
    Gt &operator-=(const Gt &rhs);

    Gt &operator*=(const scalar::Scalar &rhs);

public:
    friend inline Gt operator+(const Gt &a, const Gt &b) { return Gt(a) += b; }
    friend inline Gt operator-(const Gt &a, const Gt &b) { return Gt(a) -= b; }

    friend inline Gt operator*(const Gt &a, const scalar::Scalar &b) { return Gt(a) *= b; }

    friend inline bool operator==(const Gt &a, const Gt &b) { return a.data == b.data; }
    friend inline bool operator!=(const Gt &a, const Gt &b) { return a.data != b.data; }
};

} // namespace bls12_381::group

#endif //BLS12_381_GT_H