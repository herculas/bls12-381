#ifndef BLS12_381_FP12_H
#define BLS12_381_FP12_H

#include "fp6.h"

namespace bls12_381::field {

class Fp12 {
private:
    Fp6 c0;
    Fp6 c1;

public:
    Fp12();

    explicit Fp12(Fp &&fp);
    explicit Fp12(Fp2 &&fp);
    explicit Fp12(Fp6 &&fp);
    explicit Fp12(Fp6 &&fp0, Fp6 &&fp1);

    explicit Fp12(const Fp &fp);
    explicit Fp12(const Fp2 &fp);
    explicit Fp12(const Fp6 &fp);
    explicit Fp12(const Fp6 &fp0, const Fp6 &fp1);

    static Fp12 zero();
    static Fp12 one();
    static Fp12 random();

    [[nodiscard]] Fp6 get_c0() const;
    [[nodiscard]] Fp6 get_c1() const;

    [[nodiscard]] bool is_zero() const;

    [[nodiscard]] Fp12 square() const;
    [[nodiscard]] Fp12 conjugate() const;
    [[nodiscard]] Fp12 frobenius_map() const;
    [[nodiscard]] Fp12 mul_by_fp2(const Fp2 &fp0, const Fp2 &fp1, const Fp2 &fp4) const;

    [[nodiscard]] std::optional<Fp12> invert() const;

public:
    Fp12 operator-() const;
    Fp12 &operator=(const Fp12 &rhs);

    Fp12 &operator+=(const Fp12 &rhs);
    Fp12 &operator-=(const Fp12 &rhs);
    Fp12 &operator*=(const Fp12 &rhs);

private:
    friend inline Fp12 operator+(const Fp12 &a, const Fp12 &b) { return Fp12(a) += b; }
    friend inline Fp12 operator-(const Fp12 &a, const Fp12 &b) { return Fp12(a) -= b; }
    friend inline Fp12 operator*(const Fp12 &a, const Fp12 &b) { return Fp12(a) *= b; }

    friend inline bool operator==(const Fp12 &a, const Fp12 &b) { return a.c0 == b.c0 && a.c1 == b.c1; }
    friend inline bool operator!=(const Fp12 &a, const Fp12 &b) { return a.c0 != b.c0 || a.c1 != b.c1; }
};
} // namespace bls12_381::field

#endif //BLS12_381_FP12_H