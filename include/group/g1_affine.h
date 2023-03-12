#ifndef BLS12_381_G1_AFFINE_H
#define BLS12_381_G1_AFFINE_H

#include <array>
#include <cstdint>
#include <optional>

#include "field/fp.h"

namespace bls12_381::scalar { class Scalar; }
namespace bls12_381::group { class G1Projective; }

namespace bls12_381::group {

class G1Affine {
public:
    static constexpr int32_t WIDTH = field::Fp::WIDTH;
    static constexpr int32_t BYTE_SIZE = WIDTH * sizeof(uint64_t);
    static constexpr int32_t RAW_SIZE = BYTE_SIZE * 2 + 1;

private:
    field::Fp x;
    field::Fp y;
    bool infinity;

public:
    G1Affine();

    G1Affine(const G1Affine &point);
    explicit G1Affine(const G1Projective &point);
    explicit G1Affine(const field::Fp &x, const field::Fp &y, bool infinity);

    G1Affine(G1Affine &&point) noexcept;
    explicit G1Affine(G1Projective &&point);
    explicit G1Affine(field::Fp &&x, field::Fp &&y, bool infinity);

    static G1Affine identity() noexcept;
    static G1Affine generator() noexcept;

    static std::optional<G1Affine> from_compressed(const std::array<uint8_t, G1Affine::BYTE_SIZE> &bytes);
    static std::optional<G1Affine> from_compressed_unchecked(const std::array<uint8_t, G1Affine::BYTE_SIZE> &bytes);
    static std::optional<G1Affine> from_uncompressed(const std::array<uint8_t, G1Affine::BYTE_SIZE * 2> &bytes);
    static std::optional<G1Affine> from_uncompressed_unchecked(const std::array<uint8_t, G1Affine::BYTE_SIZE * 2> &bytes);
    static G1Affine from_slice_unchecked(const std::vector<uint8_t> &bytes);

    [[nodiscard]] const field::Fp &get_x() const noexcept;
    [[nodiscard]] const field::Fp &get_y() const noexcept;

    [[nodiscard]] bool is_identity() const;
    [[nodiscard]] bool is_on_curve() const;
    [[nodiscard]] bool is_torsion_free() const;

    [[nodiscard]] std::array<uint8_t, G1Affine::RAW_SIZE> to_raw_bytes() const;
    [[nodiscard]] std::array<uint8_t, G1Affine::BYTE_SIZE> to_compressed() const;
    [[nodiscard]] std::array<uint8_t, G1Affine::BYTE_SIZE * 2> to_uncompressed() const;

private:
    [[nodiscard]] G1Affine endomorphism() const;

public:
    G1Affine operator-() const;

    G1Affine &operator=(const G1Affine &rhs);
    G1Affine &operator=(G1Affine &&rhs) noexcept;

    G1Projective operator+(const G1Projective &rhs) const;
    G1Projective operator-(const G1Projective &rhs) const;
    G1Projective operator*(const scalar::Scalar &rhs) const;

public:
    friend inline bool operator==(const G1Affine &a, const G1Affine &b) {
        return (a.infinity & b.infinity) | ((!a.infinity) & (!b.infinity) & (a.x == b.x) & (a.y == b.y));
    }

    friend inline bool operator!=(const G1Affine &a, const G1Affine &b) {
        return !(a == b);
    }
};
} // namespace bls12_381::group

#endif //BLS12_381_G1_AFFINE_H