#include "group/g1_projective.h"

#include <cassert>

#ifdef __cpp_lib_ranges
#include <ranges>
#endif

#include "field/constant.h"
#include "group/constant.h"
#include "group/g1_affine.h"
#include "scalar/scalar.h"

namespace bls12_381::group {

using rng::core::RngCore;

using field::Fp;
using scalar::Scalar;

G1Projective::G1Projective() : x{Fp::zero()}, y{Fp::one()}, z{Fp::zero()} {}

G1Projective::G1Projective(const G1Projective &point) = default;

G1Projective::G1Projective(const G1Affine &point)
        : x{point.get_x()}, y{point.get_y()}, z{point.is_identity() ? Fp::zero() : Fp::one()} {}

G1Projective::G1Projective(const Fp &x, const Fp &y, const Fp &z)  // NOLINT(bugprone-easily-swappable-parameters)
        : x{x}, y{y}, z{z} {}

G1Projective::G1Projective(G1Projective &&point) noexcept = default;

G1Projective::G1Projective(G1Affine &&point)
        : x{point.get_x()}, y{point.get_y()}, z{point.is_identity() ? Fp::zero() : Fp::one()} {}

G1Projective::G1Projective(Fp &&x, Fp &&y, Fp &&z) : x{std::move(x)}, y{std::move(y)}, z{std::move(z)} {}

G1Projective::~G1Projective() = default;

G1Projective G1Projective::identity() noexcept {
    return G1Projective{};
}

G1Projective G1Projective::generator() noexcept {
    return G1Projective{
            Fp({
                       0x5cb38790fd530c16, 0x7817fc679976fff5, 0x154f95c7143ba1c1,
                       0xf0ae6acdf3d0e747, 0xedce6ecc21dbf440, 0x120177419e0bfb75,
               }),
            Fp({
                       0xbaac93d50ce72271, 0x8c22631a7918fd8e, 0xdd595f13570725ce,
                       0x51ac582950405194, 0x0e1c8c3fad0059c0, 0x0bbc3efc5008a26a,
               }),
            Fp::one(),
    };
}

G1Projective G1Projective::random(RngCore &rng) {
    while (true) {
        const Fp rx = Fp::random(rng);
        const bool flip_sign = rng.get_uint32() % 2 != 0;
        const auto temp = (rx.square() * rx + field::constant::B).sqrt();
        if (!temp.has_value()) continue;
        const Fp &ry = temp.value();
        const G1Affine point{rx, flip_sign ? -ry : ry, false};
        const G1Projective curve(point);
        G1Projective res = curve.clear_cofactor();
        if (!res.is_identity()) return res;
    }
}

std::vector<G1Affine> G1Projective::batch_normalize(const std::vector<G1Projective> &points) {
    std::vector<G1Affine> results(points.size());
    std::vector<Fp> temp_xs(points.size());

    Fp acc = Fp::one();
    for (int i = 0; i < points.size(); ++i) {
        temp_xs[i] = acc;
        acc = points[i].is_identity() ? acc : (acc * points[i].z);
    }

    assert(acc.invert().has_value());
    acc = acc.invert().value();

    for (int i = static_cast<int32_t>(points.size()) - 1; i >= 0; --i) {
        Fp const temp = temp_xs[i] * acc;
        acc = points[i].is_identity() ? acc : (acc * points[i].z);
        G1Affine const r{
                points[i].x * temp,
                points[i].y * temp,
                false,
        };
        results[i] = points[i].is_identity() ? G1Affine::identity() : r;
    }

    return results;
}

const Fp &G1Projective::get_x() const noexcept {
    return this->x;
}

const Fp &G1Projective::get_y() const noexcept {
    return this->y;
}

const Fp &G1Projective::get_z() const noexcept {
    return this->z;
}

bool G1Projective::is_identity() const {
    return this->z.is_zero();
}

bool G1Projective::is_on_curve() const {
    // y ^ 2 * z = x ^ 3 + b * z ^ 3.
    return ((this->y.square() * this->z) ==
            (this->x.square() * this->x + this->z.square() * this->z * field::constant::B)) ||
           this->z.is_zero();
}

Fp mul_by_3b(Fp &a) {
    a = a + a;
    a = a + a;
    return a + a + a;
}

G1Projective G1Projective::doubles() const {
    Fp t0 = this->y.square();
    Fp z3 = t0 + t0;
    z3 = z3 + z3;
    z3 = z3 + z3;
    Fp t1 = this->y * this->z;
    Fp t2 = this->z.square();
    t2 = mul_by_3b(t2);
    Fp x3 = t2 * z3;
    Fp y3 = t0 + t2;
    z3 = t1 * z3;
    t1 = t2 + t2;
    t2 = t1 + t2;
    t0 = t0 - t2;
    y3 = t0 * y3;
    y3 = x3 + y3;
    t1 = this->x * this->y;
    x3 = t0 * t1;
    x3 = x3 + x3;

    G1Projective const temp{x3, y3, z3};
    return this->is_identity() ? G1Projective::identity() : temp;
}

G1Projective G1Projective::add(const G1Projective &rhs) const {
    Fp t0 = this->x * rhs.x;
    Fp t1 = this->y * rhs.y;
    Fp t2 = this->z * rhs.z;
    Fp t3 = this->x + this->y;
    Fp t4 = rhs.x + rhs.y;
    t3 = t3 * t4;
    t4 = t0 + t1;
    t3 = t3 - t4;
    t4 = this->y + this->z;
    Fp x3 = rhs.y + rhs.z;
    t4 = t4 * x3;
    x3 = t1 + t2;
    t4 = t4 - x3;
    x3 = this->x + this->z;
    Fp y3 = rhs.x + rhs.z;
    x3 = x3 * y3;
    y3 = t0 + t2;
    y3 = x3 - y3;
    x3 = t0 + t0;
    t0 = x3 + t0;
    t2 = mul_by_3b(t2);
    Fp z3 = t1 + t2;
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
    return G1Projective{x3, y3, z3};
}

G1Projective G1Projective::add_mixed(const G1Affine &rhs) const {
    Fp t0 = this->x * rhs.get_x();
    Fp t1 = this->y * rhs.get_y();
    Fp t3 = rhs.get_x() + rhs.get_y();
    Fp t4 = this->x + this->y;
    t3 = t3 * t4;
    t4 = t0 + t1;
    t3 = t3 - t4;
    t4 = rhs.get_y() * this->z;
    t4 = t4 + this->y;
    Fp y3 = rhs.get_x() * this->z;
    y3 = y3 + this->x;
    Fp x3 = t0 + t0;
    t0 = x3 + t0;
    Fp te = this->z;
    Fp t2 = mul_by_3b(te);
    Fp z3 = t1 + t2;
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

    G1Projective const temp{x3, y3, z3};
    return rhs.is_identity() ? *this : temp;
}

G1Projective G1Projective::multiply(const std::array<uint8_t, 32> &bytes) const {
    G1Projective acc = G1Projective::identity();
#ifdef __cpp_lib_ranges
    for (const uint8_t &bit: bytes
                          | std::views::reverse
                          | std::views::transform([](uint8_t byte) {
        std::vector<uint8_t> res{};
        res.reserve(8);
        for (const auto &i: std::views::iota(0, 8) | std::views::reverse) {
            res.push_back((byte >> i) & static_cast<uint8_t>(1));
        }
        return res;
    })
                          | std::views::join
                          | std::views::drop(1)
            ) {
        acc = acc.doubles();
        if (bit != 0) acc = acc + *this;
    }
#else
    for (auto iter = bytes.rbegin(); iter != bytes.rend(); ++iter) {
        for (int i = 7; i >= 0; --i) {
            if (iter == bytes.rbegin() && i == 7) continue;
            uint8_t const bit = (*iter >> i) & static_cast<uint8_t>(1);
            acc = acc.doubles();
            if (bit != 0) acc = acc + *this;
        }
    }
#endif
    return acc;
}

G1Projective G1Projective::mul_by_x() const {
    G1Projective this_x = G1Projective::identity();
    G1Projective temp = *this;
    uint64_t sx = group::constant::BLS_X >> 1;
    while (sx != 0) {
        temp = temp.doubles();
        if (sx % 2 == 1)
            this_x += temp;
        sx >>= 1;
    }
    return -this_x;
}

G1Projective G1Projective::clear_cofactor() const {
    return *this - this->mul_by_x();
}

G1Projective G1Projective::operator-() const {
    return G1Projective{
            this->x,
            -this->y,
            this->z,
    };
}

G1Projective &G1Projective::operator=(const G1Projective &rhs) = default;

G1Projective &G1Projective::operator=(G1Projective &&rhs) noexcept = default;

G1Projective &G1Projective::operator+=(const G1Projective &rhs) {
    *this = this->add(rhs);
    return *this;
}

G1Projective &G1Projective::operator-=(const G1Projective &rhs) {
    *this = this->add(-rhs);
    return *this;
}

G1Projective &G1Projective::operator+=(const G1Affine &rhs) {
    *this = this->add_mixed(rhs);
    return *this;
}

G1Projective &G1Projective::operator-=(const G1Affine &rhs) {
    *this = this->add_mixed(-rhs);
    return *this;
}

G1Projective &G1Projective::operator*=(const Scalar &rhs) {
    *this = this->multiply(rhs.to_bytes());
    return *this;
}

G1Projective G1Projective::operator+(const G1Projective &rhs) const {
    return this->add(rhs);
}

G1Projective G1Projective::operator-(const G1Projective &rhs) const {
    return this->add(-rhs);
}

G1Projective G1Projective::operator+(const G1Affine &rhs) const {
    return this->add_mixed(rhs);
}

G1Projective G1Projective::operator-(const G1Affine &rhs) const {
    return this->add_mixed(-rhs);
}

G1Projective G1Projective::operator*(const Scalar &rhs) const {
    return this->multiply(rhs.to_bytes());
}

} // namespace bls12_381::group