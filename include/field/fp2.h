#ifndef BLS12_381_FP2_H
#define BLS12_381_FP2_H

#include "fp.h"

class Fp2 {
private:
    Fp c0;
    Fp c1;

public:
    Fp2();
    explicit Fp2(Fp fp);
    explicit Fp2(Fp fp0, Fp fp1);

    static Fp2 zero();
    static Fp2 one();
    static Fp2 random();

    [[nodiscard]] Fp getC0() const;
    [[nodiscard]] Fp getC1() const;

    [[nodiscard]] bool is_zero() const;
    [[nodiscard]] bool lexicographically_largest() const;
    [[nodiscard]] std::string getHex() const;

    [[nodiscard]] Fp2 square() const;
    [[nodiscard]] Fp2 conjugate() const;
    [[nodiscard]] Fp2 frobenius_map() const;
    [[nodiscard]] Fp2 mul_by_non_residue() const;
    [[nodiscard]] Fp2 pow_vartime(std::array<uint64_t, Fp::WIDTH> exp) const;
    [[nodiscard]] Fp2 pow_vartime_extended(std::vector<uint64_t> exp) const;

    [[nodiscard]] std::optional<Fp2> sqrt() const;
    [[nodiscard]] std::optional<Fp2> invert() const;

public:
    Fp2 &operator=(const Fp2 &rhs);
    Fp2 &operator+=(const Fp2 &rhs);
    Fp2 &operator-=(const Fp2 &rhs);
    Fp2 &operator*=(const Fp2 &rhs);

    Fp2 operator-() const;

public:
    friend inline Fp2 operator+(const Fp2 &a, const Fp2 &b) { return Fp2(a) += b; }
    friend inline Fp2 operator-(const Fp2 &a, const Fp2 &b) { return Fp2(a) -= b; }
    friend inline Fp2 operator*(const Fp2 &a, const Fp2 &b) { return Fp2(a) *= b; }

    friend inline bool operator==(const Fp2 &a, const Fp2 &b) { return a.c0 == b.c0 && a.c1 == b.c1; }
    friend inline bool operator!=(const Fp2 &a, const Fp2 &b) { return a.c0 != b.c0 || a.c1 != b.c1; }
};

#endif //BLS12_381_FP2_H
