#ifndef BLS12_381_G2_AFFINE_H
#define BLS12_381_G2_AFFINE_H

#include "field/fp.h"
#include "scalar/scalar.h"

//class G2Projective;
//
//class G2Affine {
//private:
//    Fp x;
//    Fp y;
//    bool infinity;
//
//public:
//    G2Affine();
//    explicit G2Affine(G2Projective &&point);
//    explicit G2Affine(Fp &&x, Fp &&y, bool infinity);
//
//    explicit G2Affine(const G2Projective &point);
//    explicit G2Affine(const Fp &x, const Fp &y, bool infinity);
//
//    static G2Affine identity();
//    static G2Affine generator();
//
//    static std::optional<G2Affine> from_compressed(const std::array<uint8_t, Fp::WIDTH * sizeof(uint64_t)> &bytes);
//    static std::optional<G2Affine> from_compressed_unchecked(const std::array<uint8_t, Fp::WIDTH * sizeof(uint64_t)> &bytes);
//    static std::optional<G2Affine> from_uncompressed(const std::array<uint8_t, Fp::WIDTH * sizeof(uint64_t) * 2> &bytes);
//    static std::optional<G2Affine> from_uncompressed_unchecked(const std::array<uint8_t, Fp::WIDTH * sizeof(uint64_t) * 2> &bytes);
//
//    [[nodiscard]] Fp getX() const;
//    [[nodiscard]] Fp getY() const;
//
//    [[nodiscard]] bool is_identity() const;
//    [[nodiscard]] bool is_on_curve() const;
//    [[nodiscard]] bool is_torsion_free() const;
//
//    [[nodiscard]] std::array<uint8_t, Fp::WIDTH * sizeof(uint64_t)> to_compressed() const;
//    [[nodiscard]] std::array<uint8_t, Fp::WIDTH * sizeof(uint64_t) * 2> to_uncompressed() const;
//
//    [[nodiscard]] G2Affine endomorphism() const;
//
//public:
//    G2Affine operator-() const;
//    G2Affine &operator=(const G2Affine &rhs);
//public:
//    friend G2Projective operator+(const G2Affine &a, const G2Projective &b);
//    friend G2Projective operator-(const G2Affine &a, const G2Projective &b);
//    friend G2Projective operator*(const G2Affine &a, const Scalar &b);
//
//    friend inline bool operator==(const G2Affine &a, const G2Affine &b) {
//        return (a.infinity & b.infinity) | ((!a.infinity) & (!b.infinity) & (a.x == b.x) & (a.y == b.y));
//    }
//    friend inline bool operator!=(const G2Affine &a, const G2Affine &b) {
//        return !(a == b);
//    }
//};

#endif //BLS12_381_G2_AFFINE_H
