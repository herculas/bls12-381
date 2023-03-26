#ifndef BLS12_381_G1_AFFINE_H
#define BLS12_381_G1_AFFINE_H

#include <array>
#include <cstdint>
#include <optional>

#include "field/fp.h"

namespace bls12_381::scalar { class Scalar; }
namespace bls12_381::group { class G1Projective; }

namespace bls12_381::group {

/**
 * @brief An element on the G1 curve in affine coordinate space.
 * @details The group G1 uses Fp elements for coordinates.
 * @note Values of <tt>G1Affine</tt> are guaranteed to be in a order q subgroup unless the unchecked API was misused.
 *          It is ideal to keep elements in this representation to reduce memory usage and improve performance through
 *          the use of mixed curve model arithmetic.
 */
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

    /**
     * @brief Returns the identity element of G1 in affine coordinate form.
     * @return The identity element of G1.
     * @note The identity element is the point at infinity.
     */
    static G1Affine identity() noexcept;

    /**
     * @brief Returns a fixed generator of G1 in affine coordinate form.
     * @return A chosen generator of G1.
     * @note This generator is chosen to be a simple point on the curve. To derive the generator, the lexicographically
     *          smallest valid x-coordinate and the smallest corresponding y-coordinate is chosen. The result is then
     *          scaled by the cofactor to ensure the point is not the identity.
     */
    static G1Affine generator() noexcept;

    /**
     * @brief Attempts to convert a set of bytes created by <tt>G1Affine::to_raw_bytes</tt> into an element of
     *          <tt>G1Affine</tt>.
     * @param bytes The byte vector in big-endian order of size 97 bytes.
     * @return The <tt>G1Affine</tt> value, if exists.
     * @note No check is performed and no constant time is guaranteed. The expected usage of this function is for
     *          trusted bytes where performance is critical. For secure serialization, please refer to
     *          <tt>G1Affine::from_compressed</tt> or <tt>G1Affine::from_uncompressed</tt>. After deserialization,
     *          you can check the point using <tt>G1Affine::is_on_curve</tt> and <tt>G1Affine::is_torsion_free</tt>.
     */
    static auto from_slice_unchecked(const std::vector<uint8_t> &bytes) -> G1Affine;

    /**
     * @brief Attempts to convert a byte array into an element of <tt>G1Affine</tt>.
     * @param bytes The compressed byte array in big-endian order of size 48 bytes.
     * @return The <tt>G1Affine</tt> value, if exists.
     */
    static auto from_compressed(const std::array<uint8_t, G1Affine::BYTE_SIZE> &bytes) -> std::optional<G1Affine>;
    static auto from_compressed_unchecked(const std::array<uint8_t, G1Affine::BYTE_SIZE> &bytes) -> std::optional<G1Affine>;

    /**
     * @brief Attempts to convert a byte array into an element of <tt>G1Affine</tt>.
     * @param bytes The uncompressed byte array in big-endian order of size 96 bytes.
     * @return The <tt>G1Affine</tt> value, if exists.
     */
    static auto from_uncompressed(const std::array<uint8_t, G1Affine::BYTE_SIZE * 2> &bytes) -> std::optional<G1Affine>;
    static auto from_uncompressed_unchecked(const std::array<uint8_t, G1Affine::BYTE_SIZE * 2> &bytes) -> std::optional<G1Affine>;

    [[nodiscard]] const field::Fp &get_x() const noexcept;
    [[nodiscard]] const field::Fp &get_y() const noexcept;

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
     * @brief Checks if the point is free of an h-torsion component, and so is in the order q subgroup of G1.
     * @return <tt>true</tt> if the point is in the order q subgroup, <tt>false</tt> otherwise.
     * @note This should always return <tt>true</tt> unless an unchecked API was misused.
     */
    [[nodiscard]] bool is_torsion_free() const;

    /**
     * @brief Convert the an element of <tt>G1Affine</tt> into a byte array in raw representation.
     * @details The raw representation is a compressed form with an additional byte to indicate the infinity flag.
     * @return A byte array in big-endian order of size 97 bytes.
     * @note The intended usage of this function is for trusted sets of data where performance is critical.
     */
    [[nodiscard]] std::array<uint8_t, G1Affine::RAW_SIZE> to_raw_bytes() const;

    /**
     * @brief Convert the an element of <tt>G1Affine</tt> into a compressed byte array.
     * @return A byte array in big-endian order of size 48 bytes.
     * @note In compressed form, only the x-coordinate is serialized.
     */
    [[nodiscard]] std::array<uint8_t, G1Affine::BYTE_SIZE> to_compressed() const;

    /**
     * @brief Convert the an element of <tt>G1Affine</tt> into an uncompressed byte array.
     * @return A byte array in big-endian order of size 96 bytes.
     * @note In uncompressed form, both the x- and y-coordinates are serialized.
     */
    [[nodiscard]] std::array<uint8_t, G1Affine::BYTE_SIZE * 2> to_uncompressed() const;

private:
    /**
     * @brief Computes the endomorphism of this element.
     * @details The endomorphism is defined as endomorphism(x, y) = (BETA * x, y) where BETA is a non-trivial cubic
     *          root of unity in Fp.
     * @return The endomorphism of this element.
     */
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