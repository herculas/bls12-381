#ifndef BLS12_381_FP6_H
#define BLS12_381_FP6_H

#include <optional>

#include "core/rng.h"

#include "fp.h"
#include "fp2.h"

namespace bls12_381::field {

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

    static Fp6 zero() noexcept;
    static Fp6 one() noexcept;
    static Fp6 random(rng::core::RngCore &rng);

    [[nodiscard]] Fp2 get_c0() const noexcept;
    [[nodiscard]] Fp2 get_c1() const noexcept;
    [[nodiscard]] Fp2 get_c2() const noexcept;

    [[nodiscard]] bool is_zero() const;

    [[nodiscard]] Fp6 square() const;
    [[nodiscard]] Fp6 frobenius_map() const;
    [[nodiscard]] Fp6 mul_by_non_residue() const;
    [[nodiscard]] Fp6 mul_by_fp2(const Fp2 &fp) const;
    [[nodiscard]] Fp6 mul_by_fp2(const Fp2 &fp0, const Fp2 &fp1) const;

    [[nodiscard]] std::optional<Fp6> invert() const;

private:
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