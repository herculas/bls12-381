#ifndef BLS12_381_SCALAR_H
#define BLS12_381_SCALAR_H

#include <array>
#include <cstdint>
#include <optional>

#include "group/g1_affine.h"
#include "group/g1_projective.h"

namespace bls12_381::scalar {

class Scalar {
public:
    static constexpr int32_t WIDTH = 4;
private:
    std::array<uint64_t, Scalar::WIDTH> data;

public:
    Scalar();

    explicit Scalar(uint64_t val);
    explicit Scalar(std::array<uint64_t, Scalar::WIDTH> &&data);
    explicit Scalar(const std::array<uint64_t, Scalar::WIDTH> &data);

    static Scalar zero();
    static Scalar one();
    static Scalar random();

    static Scalar montgomery_reduce(const std::array<uint64_t, Scalar::WIDTH * 2> &rs);
    static Scalar from_raw(const std::array<uint64_t, Scalar::WIDTH> &values);
    static Scalar from_bytes_wide(const std::array<uint8_t, Scalar::WIDTH * sizeof(uint64_t) * 2> &bytes);

    static std::optional<Scalar> from_bytes(const std::array<uint8_t, Scalar::WIDTH * sizeof(uint64_t)> &bytes);

    [[nodiscard]] bool is_zero() const;
    [[nodiscard]] std::string getHex() const;
    [[nodiscard]] std::array<uint8_t, Scalar::WIDTH * sizeof(uint64_t)> to_bytes() const;

    [[nodiscard]] Scalar doubles() const;
    [[nodiscard]] Scalar square() const;
    [[nodiscard]] Scalar subtract_modulus() const;
    [[nodiscard]] Scalar pow_vartime(const std::array<uint64_t, Scalar::WIDTH> &exp) const;

    [[nodiscard]] std::optional<Scalar> sqrt() const;
    [[nodiscard]] std::optional<Scalar> invert() const;

private:
    static Scalar reduce(const std::array<uint64_t, Scalar::WIDTH * 2> &limbs);

public:
    Scalar operator-() const;
    Scalar &operator=(const Scalar &rhs);

    Scalar &operator+=(const Scalar &rhs);
    Scalar &operator-=(const Scalar &rhs);
    Scalar &operator*=(const Scalar &rhs);

public:
    friend inline Scalar operator+(const Scalar &a, const Scalar &b) { return Scalar(a) += b; }
    friend inline Scalar operator-(const Scalar &a, const Scalar &b) { return Scalar(a) -= b; }
    friend inline Scalar operator*(const Scalar &a, const Scalar &b) { return Scalar(a) *= b; }

    friend group::G1Projective operator*(const Scalar &a, const group::G1Affine &b);
    friend group::G1Projective operator*(const Scalar &a, const group::G1Projective &b);

    friend inline bool operator==(const Scalar &a, const Scalar &b) { return a.data == b.data; }
    friend inline bool operator!=(const Scalar &a, const Scalar &b) { return a.data != b.data; }
};
} // namespace bls12_381::scalar

#endif //BLS12_381_SCALAR_H
