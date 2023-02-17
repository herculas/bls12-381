#include "field/fp2.h"

namespace bls12_381::field {

Fp2::Fp2() : c0{Fp::zero()}, c1{Fp::zero()} {}

Fp2::Fp2(const Fp2 &fp) = default;

Fp2::Fp2(const Fp &fp) : c0{fp}, c1{Fp::zero()} {}

Fp2::Fp2(const Fp &fp0, const Fp &fp1) : c0{fp0}, c1{fp1} {}

Fp2::Fp2(Fp2 &&fp) noexcept = default;

Fp2::Fp2(Fp &&fp) : c0{std::move(fp)}, c1{Fp::zero()} {}

Fp2::Fp2(Fp &&fp0, Fp &&fp1) : c0{std::move(fp0)}, c1{std::move(fp1)} {}

Fp2 Fp2::zero() noexcept {
    return Fp2{
            Fp::zero(),
            Fp::zero()
    };
}

Fp2 Fp2::one() noexcept {
    return Fp2{
            Fp::one(),
            Fp::zero()
    };
}

Fp2 Fp2::random() {
    return Fp2{
            Fp::random(),
            Fp::random()
    };
}

Fp Fp2::get_c0() const noexcept {
    return this->c0;
}

Fp Fp2::get_c1() const noexcept {
    return this->c1;
}

bool Fp2::is_zero() const {
    return this->c0.is_zero() && this->c1.is_zero();
}

bool Fp2::lexicographically_largest() const {
    return this->c1.lexicographically_largest() | (this->c1.is_zero() & this->c0.lexicographically_largest());
}

std::string Fp2::to_hex_str() const {
    return this->c0.to_hex_str() + " + " + this->c1.to_hex_str() + " * u";
}

Fp2 Fp2::square() const {
    const Fp a = this->c0 + this->c1;
    const Fp b = this->c0 - this->c1;
    const Fp c = this->c0 + this->c0;

    return Fp2{
            a * b,
            c * this->c1
    };
}

Fp2 Fp2::conjugate() const {
    return Fp2{
            this->c0,
            -this->c1
    };
}

Fp2 Fp2::frobenius_map() const {
    return this->conjugate();
}

Fp2 Fp2::mul_by_non_residue() const {
    return Fp2{
            this->c0 - this->c1,
            this->c0 + this->c1
    };
}

Fp2 Fp2::pow(const std::array<uint64_t, Fp::WIDTH> &exp) const {
    Fp2 res = Fp2::one();
    for (int i = Fp::WIDTH - 1; i >= 0; --i) {
        for (int j = 63; j >= 0; --j) {
            res = res.square();
            if (((exp[i] >> j) & 0x01) == 0x01) {
                res *= *this;
            }
        }
    }
    return res;
}

Fp2 Fp2::pow_extended(const std::vector<uint64_t> &exp) const {
    const auto len = static_cast<int32_t>(exp.size());
    Fp2 res = Fp2::one();
    for (int i = len - 1; i >= 0; --i) {
        for (int j = 63; j >= 0; --j) {
            res = res.square();
            if (((exp[i] >> j) & 0x01) == 0x01) {
                res *= *this;
            }
        }
    }
    return res;
}

std::optional<Fp2> Fp2::sqrt() const {
    const std::array<uint64_t, Fp::WIDTH> exp = {
            0xee7fbfffffffeaaa, 0x07aaffffac54ffff, 0xd9cc34a83dac3d89,
            0xd91dd2e13ce144af, 0x92c6e9ed90d2eb35, 0x0680447a8e5ff9a6,
    };
    const Fp2 a1 = this->pow(exp);
    const Fp2 alpha = a1.square() * *this;
    const Fp2 x0 = a1 * *this;

    Fp2 sqrt;
    if (alpha == -Fp2::one()) {
        sqrt = Fp2{-x0.c1, x0.c0};
    } else {
        const std::array<uint64_t, Fp::WIDTH> exp2 = {
                0xdcff7fffffffd555, 0x0f55ffff58a9ffff, 0xb39869507b587b12,
                0xb23ba5c279c2895f, 0x258dd3db21a5d66b, 0x0d0088f51cbff34d,
        };
        sqrt = (alpha + Fp2::one()).pow(exp2) * x0;
    }

    if (sqrt.square() == *this) {
        return sqrt;
    } else {
        return std::nullopt;
    }
}

std::optional<Fp2> Fp2::invert() const {
    const auto temp = (this->c0.square() + this->c1.square()).invert();
    if (temp.has_value()) {
        return Fp2{
                this->c0 * temp.value(),
                this->c1 * (-temp.value())
        };
    } else {
        return std::nullopt;
    }
}

Fp2 &Fp2::operator=(const Fp2 &rhs) = default;

Fp2 &Fp2::operator=(Fp2 &&rhs) noexcept = default;

Fp2 &Fp2::operator+=(const Fp2 &rhs) {
    *this = Fp2{
            this->c0 + rhs.c0,
            this->c1 + rhs.c1
    };
    return *this;
}

Fp2 &Fp2::operator-=(const Fp2 &rhs) {
    *this = Fp2{
            this->c0 - rhs.c0,
            this->c1 - rhs.c1
    };
    return *this;
}

Fp2 &Fp2::operator*=(const Fp2 &rhs) {
    *this = Fp2{
            Fp::sum_of_products({this->c0, -this->c1}, {rhs.c0, rhs.c1}),
            Fp::sum_of_products({this->c0, this->c1}, {rhs.c1, rhs.c0})
    };
    return *this;
}

Fp2 Fp2::operator-() const {
    return Fp2{
            -this->c0,
            -this->c1
    };
}

Fp2 Fp2::operator+(const Fp2 &rhs) const {
    return Fp2{this->c0 + rhs.c0, this->c1 + rhs.c1};
}

Fp2 Fp2::operator-(const Fp2 &rhs) const {
    return Fp2{this->c0 - rhs.c0, this->c1 - rhs.c1};
}

Fp2 Fp2::operator*(const Fp2 &rhs) const {
    return Fp2{
            Fp::sum_of_products({this->c0, -this->c1}, {rhs.c0, rhs.c1}),
            Fp::sum_of_products({this->c0, this->c1}, {rhs.c1, rhs.c0})
    };
}

} // namespace bls12_381::field