#ifndef BLS12_381_SCALAR_H
#define BLS12_381_SCALAR_H

#include <array>
#include <cstdint>
#include <optional>

#include "core/rng.h"

#include "group/g1_affine.h"
#include "group/g1_projective.h"
#include "group/g2_affine.h"
#include "group/g2_projective.h"

namespace bls12_381::scalar {

class Scalar {
public:
    static constexpr int32_t WIDTH = 4;
    static constexpr int32_t BYTE_SIZE = WIDTH * sizeof(uint64_t);

private:
    std::array<uint64_t, Scalar::WIDTH> data;

public:
    Scalar();

    Scalar(const Scalar &scalar);
    explicit Scalar(uint64_t val);
    explicit Scalar(const std::array<uint64_t, Scalar::WIDTH> &data);

    Scalar(Scalar &&scalar) noexcept;
    explicit Scalar(std::array<uint64_t, Scalar::WIDTH> &&data);

    static Scalar zero() noexcept;
    static Scalar one() noexcept;
    static Scalar random(rng::core::RngCore &rng);

    static Scalar montgomery_reduce(const std::array<uint64_t, Scalar::WIDTH * 2> &rs);
    static Scalar from_raw(const std::array<uint64_t, Scalar::WIDTH> &values);
    static Scalar from_bytes_wide(const std::array<uint8_t, Scalar::BYTE_SIZE * 2> &bytes);

    static std::optional<Scalar> from_bytes(const std::array<uint8_t, Scalar::BYTE_SIZE> &bytes);

    [[nodiscard]] bool is_zero() const;
    [[nodiscard]] std::string get_hex() const;
    [[nodiscard]] std::array<uint8_t, Scalar::BYTE_SIZE> to_bytes() const;

    [[nodiscard]] Scalar doubles() const;
    [[nodiscard]] Scalar square() const;
    [[nodiscard]] Scalar subtract_modulus() const;
    [[nodiscard]] Scalar pow(const std::array<uint64_t, Scalar::WIDTH> &exp) const;

    [[nodiscard]] std::optional<Scalar> sqrt() const;
    [[nodiscard]] std::optional<Scalar> invert() const;

private:
    static Scalar reduce(const std::array<uint64_t, Scalar::WIDTH * 2> &limbs);

public:
    Scalar operator-() const;
    Scalar &operator=(const Scalar &rhs);
    Scalar &operator=(Scalar &&rhs) noexcept;

    Scalar &operator+=(const Scalar &rhs);
    Scalar &operator-=(const Scalar &rhs);
    Scalar &operator*=(const Scalar &rhs);

public:
    friend inline Scalar operator+(const Scalar &a, const Scalar &b) { return Scalar(a) += b; }
    friend inline Scalar operator-(const Scalar &a, const Scalar &b) { return Scalar(a) -= b; }
    friend inline Scalar operator*(const Scalar &a, const Scalar &b) { return Scalar(a) *= b; }

    friend group::G1Projective operator*(const Scalar &a, const group::G1Affine &b);
    friend group::G1Projective operator*(const Scalar &a, const group::G1Projective &b);
    friend group::G2Projective operator*(const Scalar &a, const group::G2Affine &b);
    friend group::G2Projective operator*(const Scalar &a, const group::G2Projective &b);

    friend inline bool operator==(const Scalar &a, const Scalar &b) { return a.data == b.data; }
    friend inline bool operator!=(const Scalar &a, const Scalar &b) { return a.data != b.data; }
};
} // namespace bls12_381::scalar

#endif //BLS12_381_SCALAR_H
