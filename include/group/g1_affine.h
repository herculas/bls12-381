#ifndef BLS12_381_G1_AFFINE_H
#define BLS12_381_G1_AFFINE_H

#include "field/fp.h"

namespace bls12_381 {
namespace scalar { class Scalar; }
namespace group { class G1Projective; }
}

namespace bls12_381::group {

class G1Affine {
private:
    field::Fp x;
    field::Fp y;
    bool infinity;

public:
    G1Affine();

    explicit G1Affine(G1Projective &&point);
    explicit G1Affine(field::Fp &&x, field::Fp &&y, bool infinity);

    explicit G1Affine(const G1Projective &point);
    explicit G1Affine(const field::Fp &x, const field::Fp &y, bool infinity);

    static G1Affine identity();
    static G1Affine generator();

    static std::optional<G1Affine> from_compressed(const std::array<uint8_t, field::Fp::WIDTH * sizeof(uint64_t)> &bytes);
    static std::optional<G1Affine> from_compressed_unchecked(const std::array<uint8_t, field::Fp::WIDTH * sizeof(uint64_t)> &bytes);
    static std::optional<G1Affine> from_uncompressed(const std::array<uint8_t, field::Fp::WIDTH * sizeof(uint64_t) * 2> &bytes);
    static std::optional<G1Affine> from_uncompressed_unchecked(const std::array<uint8_t, field::Fp::WIDTH * sizeof(uint64_t) * 2> &bytes);

    [[nodiscard]] field::Fp getX() const;
    [[nodiscard]] field::Fp getY() const;

    [[nodiscard]] bool is_identity() const;
    [[nodiscard]] bool is_on_curve() const;
    [[nodiscard]] bool is_torsion_free() const;

    [[nodiscard]] std::array<uint8_t, field::Fp::WIDTH * sizeof(uint64_t)> to_compressed() const;
    [[nodiscard]] std::array<uint8_t, field::Fp::WIDTH * sizeof(uint64_t) * 2> to_uncompressed() const;

    [[nodiscard]] G1Affine endomorphism() const;

public:
    G1Affine operator-() const;
    G1Affine &operator=(const G1Affine &rhs);

public:
    friend G1Projective operator+(const G1Affine &a, const G1Projective &b);
    friend G1Projective operator-(const G1Affine &a, const G1Projective &b);
    friend G1Projective operator*(const G1Affine &a, const scalar::Scalar &b);

    friend inline bool operator==(const G1Affine &a, const G1Affine &b) {
        return (a.infinity & b.infinity) | ((!a.infinity) & (!b.infinity) & (a.x == b.x) & (a.y == b.y));
    }

    friend inline bool operator!=(const G1Affine &a, const G1Affine &b) {
        return !(a == b);
    }
};
} // namespace bls12_381::group

#endif //BLS12_381_G1_AFFINE_H