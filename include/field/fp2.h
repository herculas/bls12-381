#ifndef BLS12_381_FP2_H
#define BLS12_381_FP2_H

#include "fp.h"

namespace bls12_381::field {

class Fp2 {
private:
    Fp c0;
    Fp c1;

public:
    Fp2();

    explicit Fp2(Fp &&fp);
    explicit Fp2(Fp &&fp0, Fp &&fp1);

    explicit Fp2(const Fp &fp);
    explicit Fp2(const Fp &fp0, const Fp &fp1);

    static Fp2 zero();
    static Fp2 one();
    static Fp2 random();

    [[nodiscard]] Fp get_c0() const;
    [[nodiscard]] Fp get_c1() const;

    [[nodiscard]] bool is_zero() const;
    [[nodiscard]] bool lexicographically_largest() const;
    [[nodiscard]] std::string getHex() const;

    [[nodiscard]] Fp2 square() const;
    [[nodiscard]] Fp2 conjugate() const;
    [[nodiscard]] Fp2 frobenius_map() const;
    [[nodiscard]] Fp2 mul_by_non_residue() const;
    [[nodiscard]] Fp2 pow_vartime(const std::array<uint64_t, Fp::WIDTH> &exp) const;
    [[nodiscard]] Fp2 pow_vartime_extended(const std::vector<uint64_t> &exp) const;

    [[nodiscard]] std::optional<Fp2> sqrt() const;
    [[nodiscard]] std::optional<Fp2> invert() const;

public:
    Fp2 operator-() const;
    Fp2 &operator=(const Fp2 &rhs);

    Fp2 &operator+=(const Fp2 &rhs);
    Fp2 &operator-=(const Fp2 &rhs);
    Fp2 &operator*=(const Fp2 &rhs);

public:
    friend inline Fp2 operator+(const Fp2 &a, const Fp2 &b) { return Fp2(a) += b; }
    friend inline Fp2 operator-(const Fp2 &a, const Fp2 &b) { return Fp2(a) -= b; }
    friend inline Fp2 operator*(const Fp2 &a, const Fp2 &b) { return Fp2(a) *= b; }

    friend inline bool operator==(const Fp2 &a, const Fp2 &b) { return a.c0 == b.c0 && a.c1 == b.c1; }
    friend inline bool operator!=(const Fp2 &a, const Fp2 &b) { return a.c0 != b.c0 || a.c1 != b.c1; }
};
} // namespace bls12_381::field

#endif //BLS12_381_FP2_H