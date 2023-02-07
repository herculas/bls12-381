#ifndef BLS12_381_G2_PROJECTIVE_H
#define BLS12_381_G2_PROJECTIVE_H

#include "field/fp2.h"

namespace bls12_381::scalar { class Scalar; }
namespace bls12_381::group { class G2Affine; }

namespace bls12_381::group {

class G2Projective {
private:
    field::Fp2 x;
    field::Fp2 y;
    field::Fp2 z;

public:
    G2Projective();

    G2Projective(const G2Projective &point);
    explicit G2Projective(const G2Affine &point);
    explicit G2Projective(const field::Fp2 &x, const field::Fp2 &y, const field::Fp2 &z);

    G2Projective(G2Projective &&point) noexcept;
    explicit G2Projective(G2Affine &&point);
    explicit G2Projective(field::Fp2 &&x, field::Fp2 &&y, field::Fp2 &&z);

    static G2Projective identity() noexcept;
    static G2Projective generator() noexcept;
    static G2Projective random();

    static std::vector<G2Affine> batch_normalize(const std::vector<G2Projective> &points);

    [[nodiscard]] field::Fp2 get_x() const noexcept;
    [[nodiscard]] field::Fp2 get_y() const noexcept;
    [[nodiscard]] field::Fp2 get_z() const noexcept;

    [[nodiscard]] bool is_identity() const;
    [[nodiscard]] bool is_on_curve() const;

    [[nodiscard]] G2Projective doubles() const;
    [[nodiscard]] G2Projective add(const G2Projective &rhs) const;
    [[nodiscard]] G2Projective add_mixed(const G2Affine &rhs) const;
    [[nodiscard]] G2Projective mul_by_x() const;
    [[nodiscard]] G2Projective clear_cofactor() const;
    [[nodiscard]] G2Projective psi() const;
    [[nodiscard]] G2Projective psi2() const;

private:
    [[nodiscard]] G2Projective multiply(const std::array<uint8_t, 32> &bytes) const;

public:
    G2Projective operator-() const;
    G2Projective &operator=(const G2Projective &rhs);
    G2Projective &operator=(G2Projective &&rhs) noexcept;

    G2Projective &operator+=(const G2Projective &rhs);
    G2Projective &operator-=(const G2Projective &rhs);

    G2Projective &operator+=(const G2Affine &rhs);
    G2Projective &operator-=(const G2Affine &rhs);

    G2Projective &operator*=(const scalar::Scalar &rhs);

public:
    friend inline G2Projective operator+(const G2Projective &a, const G2Affine &b) { return G2Projective(a) += b; }
    friend inline G2Projective operator-(const G2Projective &a, const G2Affine &b) { return G2Projective(a) -= b; }

    friend inline G2Projective operator+(const G2Projective &a, const G2Projective &b) { return G2Projective(a) += b; }
    friend inline G2Projective operator-(const G2Projective &a, const G2Projective &b) { return G2Projective(a) -= b; }

    friend inline G2Projective operator*(const G2Projective &a, const scalar::Scalar &b) { return G2Projective(a) *= b; }

    friend inline bool operator==(const G2Projective &a, const G2Projective &b) {
        field::Fp2 x1 = a.x * b.z;
        field::Fp2 x2 = b.x * a.z;
        field::Fp2 y1 = a.y * b.z;
        field::Fp2 y2 = b.y * a.z;

        bool a_is_zero = a.z.is_zero();
        bool b_is_zero = b.z.is_zero();

        return (a_is_zero & b_is_zero) | ((!a_is_zero) & (!b_is_zero) & (x1 == x2) & (y1 == y2));
    }

    friend inline bool operator!=(const G2Projective &a, const G2Projective &b) {
        return !(a == b);
    }
};

} // namespace bls12_381::group

#endif //BLS12_381_G2_PROJECTIVE_H