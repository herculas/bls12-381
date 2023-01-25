#ifndef BLS12_381_FP6_H
#define BLS12_381_FP6_H

#include "fp2.h"

class Fp6 {
private:
    Fp2 c0;
    Fp2 c1;
    Fp2 c2;

public:
    Fp6();
    explicit Fp6(Fp fp);
    explicit Fp6(Fp2 fp);
    explicit Fp6(Fp2 fp0, Fp2 fp1, Fp2 fp2);

    static Fp6 zero();
    static Fp6 one();
    static Fp6 random();

    [[nodiscard]] bool is_zero() const;

    [[nodiscard]] Fp6 square() const;
    [[nodiscard]] Fp6 frobenius_map() const;
    [[nodiscard]] Fp6 mul_by_non_residue() const;
    [[nodiscard]] Fp6 mul_by_fp2(Fp2 fp) const;
    [[nodiscard]] Fp6 mul_by_fp2(Fp2 fp0, Fp2 fp1) const;

    [[nodiscard]] std::optional<Fp6> invert() const;

private:
    [[nodiscard]] Fp6 mul_interleaved(Fp6 b) const;

public:
    Fp6 &operator=(const Fp6 &rhs);
    Fp6 &operator+=(const Fp6 &rhs);
    Fp6 &operator-=(const Fp6 &rhs);
    Fp6 &operator*=(const Fp6 &rhs);

    Fp6 operator-() const;

public:
    friend inline Fp6 operator+(const Fp6 &a, const Fp6 &b) { return Fp6(a) += b; }
    friend inline Fp6 operator-(const Fp6 &a, const Fp6 &b) { return Fp6(a) -= b; }
    friend inline Fp6 operator*(const Fp6 &a, const Fp6 &b) { return Fp6(a) *= b; }

    friend inline bool operator==(const Fp6 &a, const Fp6 &b) { return a.c0 == b.c0 && a.c1 == b.c1 && a.c2 == b.c2; }
    friend inline bool operator!=(const Fp6 &a, const Fp6 &b) { return a.c0 != b.c0 || a.c1 != b.c1 || a.c2 != b.c2; }
};

#endif //BLS12_381_FP6_H
