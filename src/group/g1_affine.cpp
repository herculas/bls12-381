#include "group/g1_affine.h"

#include "utils/bit.h"

#include "field/constant.h"
#include "group/g1_projective.h"
#include "scalar/scalar.h"

namespace bls12_381::group {

using rng::util::bit::from_le_bytes;
using rng::util::bit::to_le_bytes;

using field::Fp;

G1Affine::G1Affine() : x{Fp::zero()}, y{Fp::one()}, infinity{true} {}

G1Affine::G1Affine(const G1Affine &point) = default;

G1Affine::G1Affine(const Fp &x, const Fp &y, bool infinity) : x{x}, y{y}, infinity{infinity} {}

G1Affine::G1Affine(const G1Projective &point) : x{Fp::zero()}, y{Fp::one()}, infinity{true} {
    const Fp z_inv = point.get_z().invert().value_or(Fp::zero());
    const Fp rx = point.get_x() * z_inv;
    const Fp ry = point.get_y() * z_inv;
    G1Affine temp{rx, ry, false};
    if (!z_inv.is_zero()) *this = temp;
}

G1Affine::G1Affine(G1Affine &&point) noexcept = default;

G1Affine::G1Affine(G1Projective &&point) : x{Fp::zero()}, y{Fp::one()}, infinity{true} {
    const Fp z_inv = point.get_z().invert().value_or(Fp::zero());
    const Fp rx = point.get_x() * z_inv;
    const Fp ry = point.get_y() * z_inv;
    G1Affine temp{rx, ry, false};
    if (!z_inv.is_zero()) *this = temp;
}

G1Affine::G1Affine(Fp &&x, Fp &&y, bool infinity) : x{std::move(x)}, y{std::move(y)}, infinity{infinity} {}

G1Affine G1Affine::identity() noexcept {
    return G1Affine{};
}

G1Affine G1Affine::generator() noexcept {
    return G1Affine{
            Fp({
                       0x5cb38790fd530c16, 0x7817fc679976fff5, 0x154f95c7143ba1c1,
                       0xf0ae6acdf3d0e747, 0xedce6ecc21dbf440, 0x120177419e0bfb75,
               }),
            Fp({
                       0xbaac93d50ce72271, 0x8c22631a7918fd8e, 0xdd595f13570725ce,
                       0x51ac582950405194, 0x0e1c8c3fad0059c0, 0x0bbc3efc5008a26a,
               }),
            false,
    };
}

std::optional<G1Affine>
G1Affine::from_compressed(const std::array<uint8_t, G1Affine::BYTE_SIZE> &bytes) {
    const std::optional<G1Affine> res = G1Affine::from_compressed_unchecked(bytes);
    if (!res.has_value()) return std::nullopt;
    if (!res.value().is_torsion_free()) return std::nullopt;
    return res.value();
}

std::optional<G1Affine> G1Affine::from_compressed_unchecked(const std::array<uint8_t, G1Affine::BYTE_SIZE> &bytes) {
    const bool compression_flag_set = (bytes[0] >> 7) & 1;
    const bool infinity_flat_set = (bytes[0] >> 6) & 1;
    const bool sort_flat_set = (bytes[0] >> 5) & 1;

    // try to decode the x-coordinate
    std::array<uint8_t, G1Affine::BYTE_SIZE> temp{bytes};
    temp[0] &= 0b00011111;
    std::optional<Fp> sx = Fp::from_bytes(temp);

    if (!sx.has_value()) return std::nullopt;
    if (infinity_flat_set & compression_flag_set & (!sort_flat_set) & sx.value().is_zero()) return G1Affine::identity();

    const Fp x_cord = sx.value();

    std::optional<Fp> sy = (x_cord.square() * x_cord + field::constant::B).sqrt();
    if (!sy.has_value()) return std::nullopt;

    Fp y_cord = sy.value();
    if (y_cord.lexicographically_largest() ^ sort_flat_set) y_cord = -y_cord;

    if (infinity_flat_set || (!compression_flag_set)) return std::nullopt;
    return G1Affine{x_cord, y_cord, infinity_flat_set};
}

std::optional<G1Affine> G1Affine::from_uncompressed(const std::array<uint8_t, G1Affine::BYTE_SIZE * 2> &bytes) {
    const std::optional<G1Affine> res = G1Affine::from_uncompressed_unchecked(bytes);
    if (!res.has_value()) return std::nullopt;
    if (!res.value().is_on_curve() || !res.value().is_torsion_free()) return std::nullopt;
    return res.value();
}

std::optional<G1Affine>
G1Affine::from_uncompressed_unchecked(const std::array<uint8_t, G1Affine::BYTE_SIZE * 2> &bytes) {
    const bool compression_flag_set = (bytes[0] >> 7) & 1;
    const bool infinity_flat_set = (bytes[0] >> 6) & 1;
    const bool sort_flat_set = (bytes[0] >> 5) & 1;

    // try to decode the x-coordinate
    std::array<uint8_t, G1Affine::BYTE_SIZE> temp{0};
    std::copy(bytes.begin(), bytes.begin() + G1Affine::BYTE_SIZE, temp.begin());
    temp[0] &= 0b00011111;
    const std::optional<Fp> sx = Fp::from_bytes(temp);

    // try to decode the y-coordinate
    std::copy(bytes.begin() + G1Affine::BYTE_SIZE, bytes.end(), temp.begin());
    const std::optional<Fp> sy = Fp::from_bytes(temp);

    if (!sx.has_value() || !sy.has_value()) return std::nullopt;
    if (compression_flag_set || sort_flat_set) return std::nullopt;
    if (infinity_flat_set && (!sx.value().is_zero() || !sy.value().is_zero())) return std::nullopt;

    G1Affine res{
            sx.value(),
            sy.value(),
            infinity_flat_set,
    };
    if (infinity_flat_set) {
        return G1Affine::identity();
    } else {
        return res;
    }
}

G1Affine G1Affine::from_slice_unchecked(const std::vector<uint8_t> &bytes) {
    std::array<uint64_t, Fp::WIDTH> x_bytes{};
    std::array<uint64_t, Fp::WIDTH> y_bytes{};

    std::array<uint8_t, sizeof(uint64_t)> temp_bytes{};

    for (int i = 0; i < G1Affine::RAW_SIZE - 1; i += 8) {
        const int count = i / 8;
        std::copy(bytes.begin() + i, bytes.begin() + i + 8, temp_bytes.begin());
        if (count < Fp::WIDTH)
            x_bytes[count] = from_le_bytes<uint64_t>(temp_bytes);
        else
            y_bytes[count - Fp::WIDTH] = from_le_bytes<uint64_t>(temp_bytes);
    }

    const Fp x = Fp{x_bytes};
    const Fp y = Fp{y_bytes};
    bool infinity = false;

    if (bytes.size() >= G1Affine::RAW_SIZE)
        infinity = static_cast<bool>(bytes[G1Affine::RAW_SIZE - 1]);

    return G1Affine{x, y, infinity};
}

const Fp &G1Affine::get_x() const noexcept {
    return this->x;
}

const Fp &G1Affine::get_y() const noexcept {
    return this->y;
}

bool G1Affine::is_identity() const {
    return this->infinity;
}

bool G1Affine::is_on_curve() const {
    return (this->y.square() - (this->x.square() * this->x)) == field::constant::B | this->infinity;
}

bool G1Affine::is_torsion_free() const {
    const G1Projective minus_x_squared_times_p = -G1Projective(*this).mul_by_x().mul_by_x();
    const G1Affine endomorphism_p = this->endomorphism();
    return minus_x_squared_times_p == G1Projective(endomorphism_p);
}

std::array<uint8_t, G1Affine::RAW_SIZE> G1Affine::to_raw_bytes() const {
    std::array<uint8_t, G1Affine::RAW_SIZE> bytes{};
    for (int i = 0; i < Fp::WIDTH; ++i) {
        const std::array<uint8_t, 8> x_bytes = to_le_bytes<uint64_t>(this->x.get_data()[i]);
        const std::array<uint8_t, 8> y_bytes = to_le_bytes<uint64_t>(this->y.get_data()[i]);

        std::copy(x_bytes.begin(), x_bytes.end(), bytes.begin() + i * 8);
        std::copy(y_bytes.begin(), y_bytes.end(), bytes.begin() + i * 8 + Fp::BYTE_SIZE);
    }
    bytes[G1Affine::RAW_SIZE - 1] = static_cast<uint8_t>(this->infinity);
    return bytes;
}

std::array<uint8_t, G1Affine::BYTE_SIZE> G1Affine::to_compressed() const {
    std::array<uint8_t, G1Affine::BYTE_SIZE> bytes = (this->infinity ? Fp::zero() : this->x).to_bytes();

    bytes[0] |= (static_cast<uint8_t>(1) << 7);// compression flag
    bytes[0] |= (this->infinity ? (static_cast<uint8_t>(1) << 6) : static_cast<uint8_t>(0)); // infinity flag
    bytes[0] |= (((!this->infinity) && this->y.lexicographically_largest()) ? (static_cast<uint8_t>(1) << 5)
                                                                            : static_cast<uint8_t>(0));// sort flag
    return bytes;
}

std::array<uint8_t, G1Affine::BYTE_SIZE * 2> G1Affine::to_uncompressed() const {
    std::array<uint8_t, G1Affine::BYTE_SIZE * 2> bytes{};
    const std::array<uint8_t, G1Affine::BYTE_SIZE> x_bytes = (this->infinity ? Fp::zero() : this->x).to_bytes();
    const std::array<uint8_t, G1Affine::BYTE_SIZE> y_bytes = (this->infinity ? Fp::zero() : this->y).to_bytes();

    std::copy(x_bytes.begin(), x_bytes.end(), bytes.begin());
    std::copy(y_bytes.begin(), y_bytes.end(), bytes.begin() + G1Affine::BYTE_SIZE);

    bytes[0] |= (this->infinity ? (static_cast<uint8_t>(1) << 6) : static_cast<uint8_t>(0)); // infinity flag

    return bytes;
}

G1Affine G1Affine::endomorphism() const {
    G1Affine res = *this;
    res.x *= field::constant::BETA;
    return res;
}

G1Affine &G1Affine::operator=(const G1Affine &rhs) = default;

G1Affine &G1Affine::operator=(G1Affine &&rhs) noexcept = default;

G1Affine G1Affine::operator-() const {
    return G1Affine{
            this->x,
            this->infinity ? Fp::one() : (-this->y),
            this->infinity,
    };
}

G1Projective G1Affine::operator+(const G1Projective &rhs) const {
    return G1Projective(rhs) += *this;
}

G1Projective G1Affine::operator-(const G1Projective &rhs) const {
    return -G1Projective(rhs) += *this;
}

G1Projective G1Affine::operator*(const scalar::Scalar &rhs) const {
    return G1Projective{*this} *= rhs;
}

} // namespace bls12_381::group