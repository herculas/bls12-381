#ifndef BLS12_381_G1_PROJECTIVE_H
#define BLS12_381_G1_PROJECTIVE_H

#include "field/fp.h"

class Scalar;
class G1Affine;

class G1Projective {
private:
    Fp x;
    Fp y;
    Fp z;

public:
    G1Projective();

    explicit G1Projective(G1Affine &&point);
    explicit G1Projective(Fp &&x, Fp &&y, Fp &&z);

    explicit G1Projective(const G1Affine &point);
    explicit G1Projective(const Fp &x, const Fp &y, const Fp &z);

    static G1Projective identity();
    static G1Projective generator();
    static G1Projective random();

    static std::vector<G1Affine> batch_normalize(const std::vector<G1Projective> &points);

    [[nodiscard]] Fp getX() const;
    [[nodiscard]] Fp getY() const;
    [[nodiscard]] Fp getZ() const;

    [[nodiscard]] bool is_identity() const;
    [[nodiscard]] bool is_on_curve() const;

    [[nodiscard]] G1Projective doubles() const;
    [[nodiscard]] G1Projective add(const G1Projective &rhs) const;
    [[nodiscard]] G1Projective add_mixed(const G1Affine &rhs) const;
    [[nodiscard]] G1Projective multiply(const std::array<uint8_t, 32> &bytes) const;
    [[nodiscard]] G1Projective mul_by_x() const;
    [[nodiscard]] G1Projective clear_cofactor() const;

public:
    G1Projective operator-() const;
    G1Projective &operator=(const G1Projective &rhs);

    G1Projective &operator+=(const G1Projective &rhs);
    G1Projective &operator-=(const G1Projective &rhs);

    G1Projective &operator+=(const G1Affine &rhs);
    G1Projective &operator-=(const G1Affine &rhs);

    G1Projective &operator*=(const Scalar &rhs);

public:
    friend inline G1Projective operator+(const G1Projective &a, const G1Affine &b) { return G1Projective(a) += b; }
    friend inline G1Projective operator-(const G1Projective &a, const G1Affine &b) { return G1Projective(a) -= b; }

    friend inline G1Projective operator+(const G1Projective &a, const G1Projective &b) { return G1Projective(a) += b; }
    friend inline G1Projective operator-(const G1Projective &a, const G1Projective &b) { return G1Projective(a) -= b; }

    friend inline G1Projective operator*(const G1Projective &a, const Scalar &b) { return G1Projective(a) *= b; }

    friend inline bool operator==(const G1Projective &a, const G1Projective &b) {
        Fp x1 = a.x * b.z;
        Fp x2 = b.x * a.z;
        Fp y1 = a.y * b.z;
        Fp y2 = b.y * a.z;

        bool a_is_zero = a.z.is_zero();
        bool b_is_zero = b.z.is_zero();

        return (a_is_zero & b_is_zero) | ((!a_is_zero) & (!b_is_zero) & (x1 == x2) & (y1 == y2));
    }

    friend inline bool operator!=(const G1Projective &a, const G1Projective &b) {
        return !(a == b);
    }
};

#endif //BLS12_381_G1_PROJECTIVE_H