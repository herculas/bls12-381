#include "group/g2_projective.h"

#include <cassert>

#include "field/constant.h"
#include "group/constant.h"
#include "group/g2_affine.h"
#include "scalar/scalar.h"

namespace bls12_381::group {

using rng::core::RngCore;

using field::Fp;
using field::Fp2;
using scalar::Scalar;

G2Projective::G2Projective() : x{Fp2::zero()}, y{Fp2::one()}, z{Fp2::zero()} {}

G2Projective::G2Projective(const G2Projective &point) = default;

G2Projective::G2Projective(const G2Affine &point)
        : x{point.get_x()}, y{point.get_y()}, z{point.is_identity() ? Fp2::zero() : Fp2::one()} {}

G2Projective::G2Projective(const Fp2 &x, const Fp2 &y, const Fp2 &z)  // NOLINT(bugprone-easily-swappable-parameters)
        : x{x}, y{y}, z{z} {}

G2Projective::G2Projective(G2Projective &&point) noexcept = default;

G2Projective::G2Projective(G2Affine &&point)
        : x{point.get_x()}, y{point.get_y()}, z{point.is_identity() ? Fp2::zero() : Fp2::one()} {}

G2Projective::G2Projective(Fp2 &&x, Fp2 &&y, Fp2 &&z)
        : x{std::move(x)}, y{std::move(y)}, z{std::move(z)} {}

G2Projective::~G2Projective() = default;

G2Projective G2Projective::identity() noexcept {
    return G2Projective{};
}

G2Projective G2Projective::generator() noexcept {
    return G2Projective{
            Fp2{
                    Fp({
                               0xf5f28fa202940a10, 0xb3f5fb2687b4961a, 0xa1a893b53e2ae580,
                               0x9894999d1a3caee9, 0x6f67b7631863366b, 0x058191924350bcd7,
                       }),
                    Fp({
                               0xa5a9c0759e23f606, 0xaaa0c59dbccd60c3, 0x3bb17e18e2867806,
                               0x1b1ab6cc8541b367, 0xc2b6ed0ef2158547, 0x11922a097360edf3,
                       }),
            },
            Fp2{
                    Fp({
                               0x4c730af860494c4a, 0x597cfa1f5e369c5a, 0xe7e6856caa0a635a,
                               0xbbefb5e96e0d495f, 0x07d3a975f0ef25a2, 0x0083fd8e7e80dae5,
                       }),
                    Fp({
                               0xadc0fc92df64b05d, 0x18aa270a2b1461dc, 0x86adac6a3be4eba0,
                               0x79495c4ec93da33a, 0xe7175850a43ccaed, 0x0b2bc2a163de1bf2,
                       }),
            },
            Fp2::one(),
    };
}

G2Projective G2Projective::random(RngCore &rng) {
    while (true) {
        const Fp2 rx = Fp2::random(rng);
        const bool flip_sign = rng.get_uint32() % 2 != 0;
        const auto temp = (rx.square() * rx + field::constant::B2).sqrt();
        if (!temp.has_value()) continue;
        const Fp2 &ry = temp.value();
        const G2Affine point{rx, flip_sign ? -ry : ry, false};
        const G2Projective curve(point);
        G2Projective res = curve.clear_cofactor();
        if (!res.is_identity()) return res;
    }
}

std::vector<G2Affine> G2Projective::batch_normalize(const std::vector<G2Projective> &points) {
    std::vector<G2Affine> results(points.size());
    std::vector<Fp2> temp_xs(points.size());

    Fp2 acc = Fp2::one();
    for (int i = 0; i < points.size(); ++i) {
        temp_xs[i] = acc;
        acc = points[i].is_identity() ? acc : (acc * points[i].z);
    }

    assert(acc.invert().has_value());
    acc = acc.invert().value();

    for (int i = static_cast<int32_t>(points.size()) - 1; i >= 0; --i) {
        Fp2 const temp = temp_xs[i] * acc;
        acc = points[i].is_identity() ? acc : (acc * points[i].z);
        G2Affine const r{points[i].x * temp, points[i].y * temp, false};
        results[i] = points[i].is_identity() ? G2Affine::identity() : r;
    }

    return results;
}

const Fp2 &G2Projective::get_x() const noexcept {
    return this->x;
}

const Fp2 &G2Projective::get_y() const noexcept {
    return this->y;
}

const Fp2 &G2Projective::get_z() const noexcept {
    return this->z;
}

bool G2Projective::is_identity() const {
    return this->z.is_zero();
}

bool G2Projective::is_on_curve() const {
    return ((this->y.square() * this->z) ==
            (this->x.square() * this->x + this->z.square() * this->z * field::constant::B2)) ||
           this->z.is_zero();
}

Fp2 mul_by_3b(Fp2 &a) {
    return a * field::constant::B3;
}

G2Projective G2Projective::doubles() const {
    Fp2 t0 = this->y.square();
    Fp2 z3 = t0 + t0;
    z3 = z3 + z3;
    z3 = z3 + z3;
    Fp2 t1 = this->y * this->z;
    Fp2 t2 = this->z.square();
    t2 = mul_by_3b(t2);
    Fp2 x3 = t2 * z3;
    Fp2 y3 = t0 + t2;
    z3 = t1 * z3;
    t1 = t2 + t2;
    t2 = t1 + t2;
    t0 = t0 - t2;
    y3 = t0 * y3;
    y3 = x3 + y3;
    t1 = this->x * this->y;
    x3 = t0 * t1;
    x3 = x3 + x3;

    G2Projective const temp{x3, y3, z3};
    return this->is_identity() ? G2Projective::identity() : temp;
}

G2Projective G2Projective::add(const G2Projective &rhs) const {
    Fp2 t0 = this->x * rhs.x;
    Fp2 t1 = this->y * rhs.y;
    Fp2 t2 = this->z * rhs.z;
    Fp2 t3 = this->x + this->y;
    Fp2 t4 = rhs.x + rhs.y;
    t3 = t3 * t4;
    t4 = t0 + t1;
    t3 = t3 - t4;
    t4 = this->y + this->z;
    Fp2 x3 = rhs.y + rhs.z;
    t4 = t4 * x3;
    x3 = t1 + t2;
    t4 = t4 - x3;
    x3 = this->x + this->z;
    Fp2 y3 = rhs.x + rhs.z;
    x3 = x3 * y3;
    y3 = t0 + t2;
    y3 = x3 - y3;
    x3 = t0 + t0;
    t0 = x3 + t0;
    t2 = mul_by_3b(t2);
    Fp2 z3 = t1 + t2;
    t1 = t1 - t2;
    y3 = mul_by_3b(y3);
    x3 = t4 * y3;
    t2 = t3 * t1;
    x3 = t2 - x3;
    y3 = y3 * t0;
    t1 = t1 * z3;
    y3 = t1 + y3;
    t0 = t0 * t3;
    z3 = z3 * t4;
    z3 = z3 + t0;
    return G2Projective{x3, y3, z3};
}

G2Projective G2Projective::add_mixed(const G2Affine &rhs) const {
    Fp2 t0 = this->x * rhs.get_x();
    Fp2 t1 = this->y * rhs.get_y();
    Fp2 t3 = rhs.get_x() + rhs.get_y();
    Fp2 t4 = this->x + this->y;
    t3 = t3 * t4;
    t4 = t0 + t1;
    t3 = t3 - t4;
    t4 = rhs.get_y() * this->z;
    t4 = t4 + this->y;
    Fp2 y3 = rhs.get_x() * this->z;
    y3 = y3 + this->x;
    Fp2 x3 = t0 + t0;
    t0 = x3 + t0;
    Fp2 te = this->z;
    Fp2 t2 = mul_by_3b(te);
    Fp2 z3 = t1 + t2;
    t1 = t1 - t2;
    y3 = mul_by_3b(y3);
    x3 = t4 * y3;
    t2 = t3 * t1;
    x3 = t2 - x3;
    y3 = y3 * t0;
    t1 = t1 * z3;
    y3 = t1 + y3;
    t0 = t0 * t3;
    z3 = z3 * t4;
    z3 = z3 + t0;

    G2Projective const temp{x3, y3, z3};
    return rhs.is_identity() ? *this : temp;
}

G2Projective G2Projective::mul_by_x() const {
    G2Projective this_x = G2Projective::identity();
    G2Projective temp = *this;
    uint64_t sx = constant::BLS_X >> 1;
    while (sx != 0) {
        temp = temp.doubles();
        if (sx % 2 == 1) {
            this_x += temp;
        }
        sx >>= 1;
    }
    return -this_x;
}

G2Projective G2Projective::clear_cofactor() const {
    G2Projective const t1 = this->mul_by_x();
    G2Projective const t2 = this->psi();
    return this->doubles().psi2() + (t1 + t2).mul_by_x() - t1 - t2 - *this;
}

G2Projective G2Projective::psi() const {
    Fp2 const psi_coeff_x{
            Fp::zero(),
            Fp(
                    {
                            0x890dc9e4867545c3, 0x2af322533285a5d5, 0x50880866309b7e2c,
                            0xa20d1b8c7e881024, 0x14e4f04fe2db9068, 0x14e56d3f1564853a,
                    }
            )
    };
    Fp2 const psi_coeff_y{
            Fp(
                    {
                            0x3e2f585da55c9ad1, 0x4294213d86c18183, 0x382844c88b623732,
                            0x92ad2afd19103e18, 0x1d794e4fac7cf0b9, 0x0bd592fc7d825ec8,
                    }
            ),
            Fp(
                    {
                            0x7bcfa7a25aa30fda, 0xdc17dec12a927e7c, 0x2f088dd86b4ebef1,
                            0xd1ca2087da74d4a7, 0x2da2596696cebc1d, 0x0e2b7eedbbfd87d2,
                    }
            )
    };
    return G2Projective{
            this->x.frobenius_map() * psi_coeff_x,
            this->y.frobenius_map() * psi_coeff_y,
            this->z.frobenius_map()
    };
}

G2Projective G2Projective::psi2() const {
    Fp2 const psi2_coeff_x{
            Fp(
                    {
                            0xcd03c9e48671f071, 0x5dab22461fcda5d2, 0x587042afd3851b95,
                            0x8eb60ebe01bacb9e, 0x03f97d6e83d050d2, 0x18f0206554638741,
                    }
            ),
            Fp::zero(),
    };
    return G2Projective{
            this->x * psi2_coeff_x,
            -this->y,
            this->z
    };
}

G2Projective G2Projective::multiply(const std::array<uint8_t, 32> &bytes) const {
    G2Projective acc = G2Projective::identity();
    for (auto iter = bytes.rbegin(); iter != bytes.rend(); ++iter) {
        for (int i = 7; i >= 0; --i) {
            if (iter == bytes.rbegin() && i == 7) continue;
            uint8_t const bit = (*iter >> i) & static_cast<uint8_t>(1);
            acc = acc.doubles();
            if (bit != 0) acc = acc + *this;
        }
    }
    return acc;
}

G2Projective G2Projective::operator-() const {
    return G2Projective{
            this->x,
            -this->y,
            this->z,
    };
}

G2Projective &G2Projective::operator=(const G2Projective &rhs) = default;

G2Projective &G2Projective::operator=(G2Projective &&rhs) noexcept = default;

G2Projective &G2Projective::operator+=(const G2Projective &rhs) {
    *this = this->add(rhs);
    return *this;
}

G2Projective &G2Projective::operator-=(const G2Projective &rhs) {
    *this = this->add(-rhs);
    return *this;
}

G2Projective &G2Projective::operator+=(const G2Affine &rhs) {
    *this = this->add_mixed(rhs);
    return *this;
}

G2Projective &G2Projective::operator-=(const G2Affine &rhs) {
    *this = this->add_mixed(-rhs);
    return *this;
}

G2Projective &G2Projective::operator*=(const Scalar &rhs) {
    *this = this->multiply(rhs.to_bytes());
    return *this;
}

} // namespace bls12_381::group