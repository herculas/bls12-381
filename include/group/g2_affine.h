#ifndef BLS12_381_G2_AFFINE_H
#define BLS12_381_G2_AFFINE_H

#include <array>
#include <cstdint>
#include <optional>

#include "field/fp2.h"

namespace bls12_381::scalar { class Scalar; }
namespace bls12_381::group { class G2Projective; }

namespace bls12_381::group {

/**
 * @brief An element on the G2 curve in affine coordinate space.
 * @details The group G2 uses Fp2 elements for coordinates.
 * @note   Values of <tt>G2Affine</tt> are guaranteed to be in a order q subgroup unless the unchecked API was misused.
 *          It is ideal to keep elements in this representation to reduce memory usage and improve performance through
 *          the use of mixed curve model arithmetic.
 */
class G2Affine {
public:
    static constexpr int32_t WIDTH = field::Fp::WIDTH * 2;
    static constexpr int32_t BYTE_SIZE = WIDTH * sizeof(uint64_t);
    static constexpr int32_t RAW_SIZE = BYTE_SIZE * 2 + 1;

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

    /**
     * @brief Returns the identity element of G2 in affine coordinate form.
     * @return The identity element of G2.
     * @note The identity element is the point at infinity.
     */
    static G2Affine identity() noexcept;

    /**
     * @brief Returns a fixed generator of G2 in affine coordinate form.
     * @return A chosen generator of G2.
     * @note This generator is chosen to be a simple point on the curve. To derive the generator, the lexicographically
     *          smallest valid x-coordinate and the smallest corresponding y-coordinate is chosen. The result is then
     *          scaled by the cofactor to ensure the point is not the identity.
     */
    static G2Affine generator() noexcept;

    /**
     * @brief Attempts to convert a set of bytes created by <tt>G2Affine::to_raw_bytes</tt> into an element of <tt>G1Affine</tt>.
     * @param bytes The byte vector in big-endian order of size 193 bytes.
     * @return The <tt>G2Affine</tt> value, if exists.
     * @note No check is performed and no constant time is guaranteed. The expected usage of this function is for
     *          trusted bytes where performance is critical. For secure serialization, please refer to
     *          <tt>G2Affine::from_compressed</tt> or <tt>G2Affine::from_uncompressed</tt>. After deserialization,
     *          you can check the point using <tt>G2Affine::is_on_curve</tt> and <tt>G2Affine::is_torsion_free</tt>.
     */
    static auto from_slice_unchecked(const std::vector<uint8_t> &bytes) -> G2Affine;

    /**
     * @brief Attempts to convert a byte array into an element of <tt>G2Affine</tt>.
     * @param bytes The compressed byte array in big-endian order of size 96 bytes.
     * @return The <tt>G2Affine</tt> value, if exists.
     */
    static std::optional<G2Affine> from_compressed(const std::array<uint8_t, G2Affine::BYTE_SIZE> &bytes);
    static std::optional<G2Affine> from_compressed_unchecked(const std::array<uint8_t, G2Affine::BYTE_SIZE> &bytes);

    /**
     * @brief Attempts to convert a byte array into an element of <tt>G2Affine</tt>.
     * @param bytes The uncompressed byte array in big-endian order of size 192 bytes.
     * @return The <tt>G2Affine</tt> value, if exists.
     */
    static std::optional<G2Affine> from_uncompressed(const std::array<uint8_t, G2Affine::BYTE_SIZE * 2> &bytes);
    static std::optional<G2Affine> from_uncompressed_unchecked(const std::array<uint8_t, G2Affine::BYTE_SIZE * 2> &bytes);

    [[nodiscard]] const field::Fp2 &get_x() const noexcept;
    [[nodiscard]] const field::Fp2 &get_y() const noexcept;

    /**
     * @brief Checks if the point is the identity element, i.e., the point at infinity.
     * @return <tt>true</tt> if the point is the identity element, <tt>false</tt> otherwise.
     */
    [[nodiscard]] bool is_identity() const;

    /**
     * @brief Checks if the point is on the curve.
     * @return <tt>true</tt> if the point is on the curve, <tt>false</tt> otherwise.
     * @note This should always return <tt>true</tt> unless an unchecked API was misused.
     */
    [[nodiscard]] bool is_on_curve() const;

    /**
     * @brief Checks if the point is free of an h-torsion component, and so is in the order q subgroup of G2.
     * @return <tt>true</tt> if the point is in the order q subgroup, <tt>false</tt> otherwise.
     * @note This should always return <tt>true</tt> unless an unchecked API was misused.
     */
    [[nodiscard]] bool is_torsion_free() const;

    /**
     * @brief Convert the an element of <tt>G2Affine</tt> into a byte array in raw representation.
     * @details The raw representation is a compressed form with an additional byte to indicate the infinity flag.
     * @return A byte array in big-endian order of size 193 bytes.
     * @note The intended usage of this function is for trusted sets of data where performance is critical.
     */
    [[nodiscard]] std::array<uint8_t, G2Affine::RAW_SIZE> to_raw_bytes() const;

    /**
     * @brief Convert the an element of <tt>G2Affine</tt> into a compressed byte array.
     * @return A byte array in big-endian order of size 96 bytes.
     * @note In compressed form, only the x-coordinate is serialized.
     */
    [[nodiscard]] std::array<uint8_t, G2Affine::BYTE_SIZE> to_compressed() const;

    /**
     * @brief Convert the an element of <tt>G2Affine</tt> into an uncompressed byte array.
     * @return A byte array in big-endian order of size 192 bytes.
     * @note In uncompressed form, both the x- and y-coordinates are serialized.
     */
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