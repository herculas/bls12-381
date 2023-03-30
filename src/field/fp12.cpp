#include "field/fp12.h"

namespace bls12_381::field {

using rng::core::RngCore;

Fp12::Fp12() noexcept: c0{Fp6::zero()}, c1{Fp6::zero()} {}

Fp12::Fp12(const Fp12 &fp) noexcept = default;

Fp12::Fp12(const Fp &fp) noexcept: c0{Fp6(fp)}, c1{Fp6::zero()} {}

Fp12::Fp12(const Fp2 &fp) noexcept: c0{Fp6(fp)}, c1{Fp6::zero()} {}

Fp12::Fp12(const Fp6 &fp) noexcept: c0{fp}, c1{Fp6::zero()} {}

Fp12::Fp12(const Fp6 &fp0, const Fp6 &fp1) noexcept: c0{fp0}, c1{fp1} {}

Fp12::Fp12(Fp12 &&fp) noexcept = default;

Fp12::Fp12(Fp &&fp) noexcept: c0{Fp6{fp}}, c1{Fp6::zero()} {}

Fp12::Fp12(Fp2 &&fp) noexcept: c0{Fp6(fp)}, c1{Fp6::zero()} {}

Fp12::Fp12(Fp6 &&fp) noexcept: c0{std::move(fp)}, c1{Fp6::zero()} {}

Fp12::Fp12(Fp6 &&fp0, Fp6 &&fp1) noexcept: c0{std::move(fp0)}, c1{std::move(fp1)} {}

Fp12::~Fp12() noexcept = default;

Fp12 Fp12::zero() noexcept {
    return Fp12{
            Fp6::zero(),
            Fp6::zero(),
    };
}

Fp12 Fp12::one() noexcept {
    return Fp12{
            Fp6::one(),
            Fp6::zero(),
    };
}

Fp12 Fp12::random(RngCore &rng) {
    return Fp12{
            Fp6::random(rng),
            Fp6::random(rng),
    };
}

const Fp6 &Fp12::get_c0() const noexcept {
    return this->c0;
}

const Fp6 &Fp12::get_c1() const noexcept {
    return this->c1;
}

bool Fp12::is_zero() const {
    return this->c0.is_zero() && this->c1.is_zero();
}

Fp12 Fp12::square() const {
    const Fp6 mul = this->c0 * this->c1;
    const Fp6 add = this->c0 + this->c1;

    Fp6 const s0 = (this->c1.mul_by_non_residue() + this->c0) * add - mul - mul.mul_by_non_residue();
    Fp6 const s1 = mul + mul;

    return Fp12{s0, s1};
}

Fp12 Fp12::conjugate() const {
    return Fp12{
            this->c0,
            -this->c1,
    };
}

Fp12 Fp12::frobenius_map() const {
    Fp6 const s0 = this->c0.frobenius_map();
    Fp6 const s1 = this->c1.frobenius_map();
    Fp6 const temp{
            Fp2{
                    Fp({
                               0x07089552b319d465, 0xc6695f92b50a8313, 0x97e83cccd117228f,
                               0xa35baecab2dc29ee, 0x1ce393ea5daace4d, 0x08f2220fb0fb66eb,
                       }),
                    Fp({
                               0xb2f66aad4ce5d646, 0x5842a06bfc497cec, 0xcf4895d42599d394,
                               0xc11b9cba40a8e8d0, 0x2e3813cbe5a0de89, 0x110eefda88847faf,
                       }),
            }
    };
    return Fp12{s0, s1 * temp};
}

Fp12 Fp12::mul_by_fp2(const Fp2 &fp0, const Fp2 &fp1, const Fp2 &fp4) const {
    const Fp6 aa = this->c0.mul_by_fp2(fp0, fp1);
    const Fp6 bb = this->c1.mul_by_fp2(fp4);

    Fp6 const s0 = bb.mul_by_non_residue() + aa;
    Fp6 const s1 = (this->c1 + this->c0).mul_by_fp2(fp0, fp1 + fp4) - aa - bb;

    return Fp12{s0, s1};
}

std::optional<Fp12> Fp12::invert() const {
    Fp6 const temp = this->c0.square() - this->c1.square().mul_by_non_residue();
    std::optional<Fp6> res = temp.invert();
    if (!res.has_value()) return std::nullopt;
    return Fp12{
            this->c0 * res.value(),
            this->c1 * (-res.value())
    };
}

Fp12 &Fp12::operator=(const Fp12 &rhs) = default;

Fp12 &Fp12::operator=(Fp12 &&rhs) noexcept = default;

Fp12 &Fp12::operator+=(const Fp12 &rhs) {
    *this = Fp12{
            this->c0 + rhs.c0,
            this->c1 + rhs.c1,
    };
    return *this;
}

Fp12 &Fp12::operator-=(const Fp12 &rhs) {
    *this = Fp12{
            this->c0 - rhs.c0,
            this->c1 - rhs.c1,
    };
    return *this;
}

Fp12 &Fp12::operator*=(const Fp12 &rhs) {
    const Fp6 aa = this->c0 * rhs.c0;
    const Fp6 bb = this->c1 * rhs.c1;

    Fp6 const s0 = bb.mul_by_non_residue() + aa;
    Fp6 const s1 = (this->c1 + this->c0) * (rhs.c0 + rhs.c1) - aa - bb;

    *this = Fp12{s0, s1};
    return *this;
}

Fp12 Fp12::operator-() const {
    return Fp12{
            -this->c0,
            -this->c1,
    };
}

Fp12 Fp12::operator+(const Fp12 &rhs) const {
    return Fp12{this->c0 + rhs.c0, this->c1 + rhs.c1};
}

Fp12 Fp12::operator-(const Fp12 &rhs) const {
    return Fp12{this->c0 - rhs.c0, this->c1 - rhs.c1};
}

Fp12 Fp12::operator*(const Fp12 &rhs) const {
    const Fp6 aa = this->c0 * rhs.c0;
    const Fp6 bb = this->c1 * rhs.c1;

    Fp6 const s0 = bb.mul_by_non_residue() + aa;
    Fp6 const s1 = (this->c1 + this->c0) * (rhs.c0 + rhs.c1) - aa - bb;

    return Fp12{s0, s1};
}

} // namespace bls12_381::field