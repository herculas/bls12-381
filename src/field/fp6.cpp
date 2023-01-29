#include "field/fp6.h"

namespace bls12_381::field {

Fp6::Fp6() : c0{Fp2::zero()}, c1{Fp2::zero()}, c2{Fp2::zero()} {}

Fp6::Fp6(Fp &&fp) : c0{Fp2(fp)}, c1{Fp2::zero()}, c2{Fp2::zero()} {}

Fp6::Fp6(Fp2 &&fp) : c0{fp}, c1{Fp2::zero()}, c2{Fp2::zero()} {}

Fp6::Fp6(Fp2 &&fp0, Fp2 &&fp1, Fp2 &&fp2) : c0{fp0}, c1{fp1}, c2{fp2} {}

Fp6::Fp6(const Fp &fp) : c0{Fp2(fp)}, c1{Fp2::zero()}, c2{Fp2::zero()} {}

Fp6::Fp6(const Fp2 &fp) : c0{fp}, c1{Fp2::zero()}, c2{Fp2::zero()} {}

Fp6::Fp6(const Fp2 &fp0, const Fp2 &fp1, const Fp2 &fp2) : c0{fp0}, c1{fp1}, c2{fp2} {}

Fp6 Fp6::zero() {
    return Fp6{
            Fp2::zero(),
            Fp2::zero(),
            Fp2::zero()
    };
}

Fp6 Fp6::one() {
    return Fp6{
            Fp2::one(),
            Fp2::zero(),
            Fp2::zero()
    };
}

Fp6 Fp6::random() {
    return Fp6{
            Fp2::random(),
            Fp2::random(),
            Fp2::random()
    };
}

bool Fp6::is_zero() const {
    return this->c0.is_zero() && this->c1.is_zero() && this->c2.is_zero();
}

Fp6 Fp6::square() const {
    Fp2 s0 = this->c0.square();
    Fp2 ab = this->c0 * this->c1;
    Fp2 s1 = ab + ab;
    Fp2 s2 = (this->c0 - this->c1 + this->c2).square();
    Fp2 bc = this->c1 * this->c2;
    Fp2 s3 = bc + bc;
    Fp2 s4 = this->c2.square();

    return Fp6{
            s3.mul_by_non_residue() + s0,
            s4.mul_by_non_residue() + s1,
            s1 + s2 + s3 - s0 - s4,
    };
}

Fp6 Fp6::mul_by_fp2(const Fp2 &fp) const {
    return Fp6{
            (this->c2 * fp).mul_by_non_residue(),
            this->c0 * fp,
            this->c1 * fp
    };
}

Fp6 Fp6::mul_by_fp2(const Fp2 &fp0, const Fp2 &fp1) const {
    Fp2 a_a = this->c0 * fp0;
    Fp2 b_b = this->c1 * fp1;

    Fp2 t1 = (this->c2 * fp1).mul_by_non_residue() + a_a;
    Fp2 t2 = (fp0 + fp1) * (this->c0 + this->c1) - a_a - b_b;
    Fp2 t3 = this->c2 * fp0 + b_b;

    return Fp6{t1, t2, t3};
}

Fp6 Fp6::mul_interleaved(const Fp6 &b) const {
    Fp6 a = *this;
    Fp b10_p_b11 = b.c1.get_c0() + b.c1.get_c1();
    Fp b10_m_b11 = b.c1.get_c0() - b.c1.get_c1();
    Fp b20_p_b21 = b.c2.get_c0() + b.c2.get_c1();
    Fp b20_m_b21 = b.c2.get_c0() - b.c2.get_c1();

    auto c000 = {a.c0.get_c0(), -a.c0.get_c1(), a.c1.get_c0(), -a.c1.get_c1(), a.c2.get_c0(), -a.c2.get_c1()};
    auto c001 = {b.c0.get_c0(), b.c0.get_c1(), b20_m_b21, b20_p_b21, b10_m_b11, b10_p_b11};
    auto c010 = {a.c0.get_c0(), a.c0.get_c1(), a.c1.get_c0(), a.c1.get_c1(), a.c2.get_c0(), a.c2.get_c1()};
    auto c011 = {b.c0.get_c1(), b.c0.get_c0(), b20_p_b21, b20_m_b21, b10_p_b11, b10_m_b11};

    auto c100 = {a.c0.get_c0(), -a.c0.get_c1(), a.c1.get_c0(), -a.c1.get_c1(), a.c2.get_c0(), -a.c2.get_c1()};
    auto c101 = {b.c1.get_c0(), b.c1.get_c1(), b.c0.get_c0(), b.c0.get_c1(), b20_m_b21, b20_p_b21};
    auto c110 = {a.c0.get_c0(), a.c0.get_c1(), a.c1.get_c0(), a.c1.get_c1(), a.c2.get_c0(), a.c2.get_c1()};
    auto c111 = {b.c1.get_c1(), b.c1.get_c0(), b.c0.get_c1(), b.c0.get_c0(), b20_p_b21, b20_m_b21};

    auto c200 = {a.c0.get_c0(), -a.c0.get_c1(), a.c1.get_c0(), -a.c1.get_c1(), a.c2.get_c0(), -a.c2.get_c1()};
    auto c201 = {b.c2.get_c0(), b.c2.get_c1(), b.c1.get_c0(), b.c1.get_c1(), b.c0.get_c0(), b.c0.get_c1()};
    auto c210 = {a.c0.get_c0(), a.c0.get_c1(), a.c1.get_c0(), a.c1.get_c1(), a.c2.get_c0(), a.c2.get_c1()};
    auto c211 = {b.c2.get_c1(), b.c2.get_c0(), b.c1.get_c1(), b.c1.get_c0(), b.c0.get_c1(), b.c0.get_c0()};

    return Fp6{
            Fp2{
                    Fp::sum_of_products(c000, c001),
                    Fp::sum_of_products(c010, c011),
            },
            Fp2{
                    Fp::sum_of_products(c100, c101),
                    Fp::sum_of_products(c110, c111),
            },
            Fp2{
                    Fp::sum_of_products(c200, c201),
                    Fp::sum_of_products(c210, c211),
            }
    };
}

Fp6 Fp6::frobenius_map() const {
    Fp2 fp0 = this->c0.frobenius_map();
    Fp2 fp1 = this->c1.frobenius_map();
    Fp2 fp2 = this->c2.frobenius_map();

    fp1 = fp1 * Fp2{
            Fp::zero(),
            Fp({
                       0xcd03c9e48671f071, 0x5dab22461fcda5d2, 0x587042afd3851b95,
                       0x8eb60ebe01bacb9e, 0x03f97d6e83d050d2, 0x18f0206554638741,
               })
    };

    fp2 = fp2 * Fp2{
            Fp({
                       0x890dc9e4867545c3, 0x2af322533285a5d5, 0x50880866309b7e2c,
                       0xa20d1b8c7e881024, 0x14e4f04fe2db9068, 0x14e56d3f1564853a,
               }),
            Fp::zero()
    };
    return Fp6{fp0, fp1, fp2};
}

Fp6 Fp6::mul_by_non_residue() const {
    return Fp6{
            this->c2.mul_by_non_residue(),
            this->c0,
            this->c1
    };
}

std::optional<Fp6> Fp6::invert() const {
    Fp2 s0 = this->c0.square() - (this->c1 * this->c2).mul_by_non_residue();
    Fp2 s1 = this->c2.square().mul_by_non_residue() - (this->c0 * this->c1);
    Fp2 s2 = this->c1.square() - (this->c0 * this->c2);

    Fp2 temp = ((this->c1 * s2) + (this->c2 * s1)).mul_by_non_residue() + (this->c0 * s0);
    std::optional<Fp2> res = temp.invert();
    if (!res.has_value()) return std::nullopt;
    return Fp6{
            res.value() * s0,
            res.value() * s1,
            res.value() * s2,
    };
}

Fp6 &Fp6::operator=(const Fp6 &rhs) {
    if (this == &rhs) return *this;
    this->c0 = rhs.c0;
    this->c1 = rhs.c1;
    this->c2 = rhs.c2;
    return *this;
}

Fp6 &Fp6::operator+=(const Fp6 &rhs) {
    *this = Fp6{
            this->c0 + rhs.c0,
            this->c1 + rhs.c1,
            this->c2 + rhs.c2
    };
    return *this;
}

Fp6 &Fp6::operator-=(const Fp6 &rhs) {
    *this = Fp6{
            this->c0 - rhs.c0,
            this->c1 - rhs.c1,
            this->c2 - rhs.c2
    };
    return *this;
}

Fp6 &Fp6::operator*=(const Fp6 &rhs) {
    *this = this->mul_interleaved(rhs);
    return *this;
}

Fp6 Fp6::operator-() const {
    return Fp6{
            -this->c0,
            -this->c1,
            -this->c2
    };
}

} // namespace bls12_381::field