#ifndef BLS12_381_G2_AFFINE_H
#define BLS12_381_G2_AFFINE_H

#include <array>
#include <cstdint>
#include <optional>

#include "field/fp2.h"

namespace bls12_381::scalar { class Scalar; }
namespace bls12_381::group { class G2Projective; }

namespace bls12_381::group {

class G2Affine {
public:
    static constexpr int32_t WIDTH = field::Fp::WIDTH * 2;
    static constexpr int32_t BYTE_SIZE = WIDTH * sizeof(uint64_t);

private:
    field::Fp2 x;
    field::Fp2 y;
    bool infinity;

public:
    G2Affine();

    G2Affine(const G2Affine &point);
    explicit G2Affine(const G2Projective &point);
    explicit G2Affine(const field::Fp2 &x, const field::Fp2 &y, bool infinity);

    G2Affine(G2Affine &&point) noexcept;
    explicit G2Affine(G2Projective &&point);
    explicit G2Affine(field::Fp2 &&x, field::Fp2 &&y, bool infinity);

    static G2Affine identity() noexcept;
    static G2Affine generator() noexcept;

    static std::optional<G2Affine> from_compressed(const std::array<uint8_t, G2Affine::BYTE_SIZE> &bytes);
    static std::optional<G2Affine> from_compressed_unchecked(const std::array<uint8_t, G2Affine::BYTE_SIZE> &bytes);
    static std::optional<G2Affine> from_uncompressed(const std::array<uint8_t, G2Affine::BYTE_SIZE * 2> &bytes);
    static std::optional<G2Affine> from_uncompressed_unchecked(const std::array<uint8_t, G2Affine::BYTE_SIZE * 2> &bytes);

    [[nodiscard]] field::Fp2 get_x() const noexcept;
    [[nodiscard]] field::Fp2 get_y() const noexcept;

    [[nodiscard]] bool is_identity() const;
    [[nodiscard]] bool is_on_curve() const;
    [[nodiscard]] bool is_torsion_free() const;

    [[nodiscard]] std::array<uint8_t, G2Affine::BYTE_SIZE> to_compressed() const;
    [[nodiscard]] std::array<uint8_t, G2Affine::BYTE_SIZE * 2> to_uncompressed() const;

public:
    G2Affine operator-() const;
    G2Affine &operator=(const G2Affine &rhs);
    G2Affine &operator=(G2Affine &&rhs) noexcept;

public:
    friend G2Projective operator+(const G2Affine &a, const G2Projective &b);
    friend G2Projective operator-(const G2Affine &a, const G2Projective &b);
    friend G2Projective operator*(const G2Affine &a, const scalar::Scalar &b);

    friend inline bool operator==(const G2Affine &a, const G2Affine &b) {
        return (a.infinity & b.infinity) | ((!a.infinity) & (!b.infinity) & (a.x == b.x) & (a.y == b.y));
    }
    friend inline bool operator!=(const G2Affine &a, const G2Affine &b) {
        return !(a == b);
    }
};

} // namespace bls12_381::group

#endif //BLS12_381_G2_AFFINE_H