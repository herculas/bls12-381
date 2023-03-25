#ifndef BLS12_381_FP_H
#define BLS12_381_FP_H

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "core/rng.h"

namespace bls12_381::field {

/**
 * @brief Represents an element of the field Fp.
 * @details The internal representation is six 64-bit unsigned integers in little-endian order. <tt>Fp</tt> values are
 *          always in Montgomery form, i.e., Fp(x) = xR mod p, where R = 2 ^ 384.
 */
class Fp {
public:
    static constexpr int32_t WIDTH = 6;
    static constexpr int32_t BYTE_SIZE = Fp::WIDTH * sizeof(uint64_t);

private:
    std::array<uint64_t, Fp::WIDTH> data;

public:
    Fp();

    Fp(const Fp &fp);
    explicit Fp(uint64_t val);
    explicit Fp(const std::array<uint64_t, Fp::WIDTH> &data);

    Fp(Fp &&fp) noexcept;
    explicit Fp(std::array<uint64_t, Fp::WIDTH> &&data);

    /**
     * @brief Returns zero, the additive identity.
     * @return The additive identity.
     */
    static Fp zero() noexcept;

    /**
     * @brief Returns one, the multiplicative identity.
     * @return The multiplicative identity.
     */
    static Fp one() noexcept;
    static Fp random(rng::core::RngCore &rng);

    static Fp montgomery_reduce(const std::array<uint64_t, Fp::WIDTH * 2> &ts);

    /**
     * @brief Computes the inner-product of two series of <tt>Fp</tt> elements.
     * @param a A series of <tt>Fp</tt> elements.
     * @param b A series of <tt>Fp</tt> elements.
     * @return The inner-product of the A and B.
     */
    static Fp sum_of_products(const std::vector<Fp> &a, const std::vector<Fp> &b);

    /**
     * @brief Attempts to convert byte array into an <tt>Fp</tt> value, failing if the input is not canonical.
     * @param bytes The byte array in big-endian order of size 48 bytes.
     * @return The <tt>Fp</tt> value, if the input is canonical.
     */
    static auto from_bytes(const std::array<uint8_t, Fp::BYTE_SIZE> &bytes) -> std::optional<Fp>;

    [[nodiscard]] const std::array<uint64_t, Fp::WIDTH> &get_data() const;

    [[nodiscard]] bool is_zero() const;

    /**
     * @brief Checks if this <tt>Fp</tt> element is strictly lexicographically larger than its negation.
     * @return <tt>true</tt> if this element is lexicographically larger than its negation, <tt>false</tt> otherwise.
     */
    [[nodiscard]] bool lexicographically_largest() const;

    [[nodiscard]] std::string to_hex_str() const;

    /**
     * @brief Converts an element of <tt>Fp</tt> into a byte array.
     * @return A byte array in big-endian order of size 48 bytes.
     */
    [[nodiscard]] auto to_bytes() const -> std::array<uint8_t, Fp::BYTE_SIZE>;

    [[nodiscard]] Fp square() const;
    [[nodiscard]] Fp subtract_modulus() const;

    /**
     * @brief Computes the power of this <tt>Fp</tt> element.
     * @param exp The exponent.
     * @return The power of this element to the given exponent.
     */
    [[nodiscard]] Fp pow(const std::array<uint64_t, Fp::WIDTH> &exp) const;

    /**
     * @brief Computes the square root of this <tt>Fp</tt> element.
     * @details This method uses the Tonelli-Shank's algorithm, as p = 3 (mod 4). This means that we only need to
     *          exponentiate the element by (p + 1) / 4. This only works if the element is a quadratic residue, so we
     *          should check if we get a correct result at the end.
     * @return The square root of this element, if exists.
     */
    [[nodiscard]] auto sqrt() const -> std::optional<Fp>;

    /**
     * @brief Computes the multiplicative inverse of this <tt>Fp</tt> element.
     * @return The multiplicative inverse of this element, if this element is non-zero.
     */
    [[nodiscard]] auto invert() const -> std::optional<Fp>;

private:
    static Fp reduce(const std::array<uint64_t, Fp::WIDTH * 2> &limbs);

public:
    Fp operator-() const;
    Fp &operator=(Fp &&rhs) noexcept;
    Fp &operator=(const Fp &rhs);

    Fp &operator+=(const Fp &rhs);
    Fp &operator-=(const Fp &rhs);
    Fp &operator*=(const Fp &rhs);

    Fp operator+(const Fp &rhs) const;
    Fp operator-(const Fp &rhs) const;
    Fp operator*(const Fp &rhs) const;

public:
    friend inline bool operator==(const Fp &a, const Fp &b) { return a.data == b.data; }
    friend inline bool operator!=(const Fp &a, const Fp &b) { return a.data != b.data; }
};

} // namespace bls12_381::field

#endif //BLS12_381_FP_H