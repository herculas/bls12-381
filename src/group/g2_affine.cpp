#include "group/g2_affine.h"

#include "field/constant.h"
#include "group/g2_projective.h"
#include "scalar/scalar.h"

namespace bls12_381::group {

G2Affine::G2Affine() : x{field::Fp2::zero()}, y{field::Fp2::one()}, infinity{true} {}

G2Affine::G2Affine(G2Projective &&point) : x{field::Fp2::zero()}, y{field::Fp2::one()}, infinity{true} {
    field::Fp2 z_inv = point.get_z().invert().value_or(field::Fp2::zero());
    field::Fp2 rx = point.get_x() * z_inv;
    field::Fp2 ry = point.get_y() * z_inv;
    G2Affine temp{rx, ry, false};
    if (!z_inv.is_zero()) *this = temp;
}

G2Affine::G2Affine(field::Fp2 &&x, field::Fp2 &&y, bool infinity) : x{x}, y{y}, infinity{infinity} {}

G2Affine::G2Affine(const G2Projective &point) : x{field::Fp2::zero()}, y{field::Fp2::one()}, infinity{true} {
    field::Fp2 z_inv = point.get_z().invert().value_or(field::Fp2::zero());
    field::Fp2 rx = point.get_x() * z_inv;
    field::Fp2 ry = point.get_y() * z_inv;
    G2Affine temp{rx, ry, false};
    if (!z_inv.is_zero()) *this = temp;
}

G2Affine::G2Affine(const field::Fp2 &x, const field::Fp2 &y, bool infinity) : x{x}, y{y}, infinity{infinity} {}

G2Affine G2Affine::identity() {
    return G2Affine{};
}

G2Affine G2Affine::generator() {
    return G2Affine{
            field::Fp2{
                    field::Fp({
                                      0xf5f28fa202940a10, 0xb3f5fb2687b4961a, 0xa1a893b53e2ae580,
                                      0x9894999d1a3caee9, 0x6f67b7631863366b, 0x058191924350bcd7,
                              }),
                    field::Fp({
                                      0xa5a9c0759e23f606, 0xaaa0c59dbccd60c3, 0x3bb17e18e2867806,
                                      0x1b1ab6cc8541b367, 0xc2b6ed0ef2158547, 0x11922a097360edf3,
                              }),
            },
            field::Fp2{
                    field::Fp({
                                      0x4c730af860494c4a, 0x597cfa1f5e369c5a, 0xe7e6856caa0a635a,
                                      0xbbefb5e96e0d495f, 0x07d3a975f0ef25a2, 0x0083fd8e7e80dae5,
                              }),
                    field::Fp({
                                      0xadc0fc92df64b05d, 0x18aa270a2b1461dc, 0x86adac6a3be4eba0,
                                      0x79495c4ec93da33a, 0xe7175850a43ccaed, 0x0b2bc2a163de1bf2,
                              }),
            },
            false
    };
}

std::optional<G2Affine> G2Affine::from_compressed(const std::array<uint8_t, G2Affine::BYTE_SIZE> &bytes) {
    std::optional<G2Affine> res = G2Affine::from_compressed_unchecked(bytes);
    if (!res.has_value()) return std::nullopt;
    if (!res.value().is_torsion_free()) return std::nullopt;
    return res.value();
}

std::optional<G2Affine> G2Affine::from_compressed_unchecked(const std::array<uint8_t, G2Affine::BYTE_SIZE> &bytes) {
    bool compression_flag_set = (bytes[0] >> 7) & 1;
    bool infinity_flat_set = (bytes[0] >> 6) & 1;
    bool sort_flat_set = (bytes[0] >> 5) & 1;

    // try to decode the x-coordinate
    std::array<uint8_t, field::Fp::BYTE_SIZE> temp{0};
    std::copy(bytes.begin(), bytes.begin() + field::Fp::BYTE_SIZE, temp.begin());
    temp[0] &= 0b00011111;
    std::optional<field::Fp> sx_c1 = field::Fp::from_bytes(temp);

    std::copy(bytes.begin() + field::Fp::BYTE_SIZE, bytes.end(), temp.begin());
    std::optional<field::Fp> sx_c0 = field::Fp::from_bytes(temp);

    if (!sx_c0.has_value() || !sx_c1.has_value()) return std::nullopt;
    field::Fp2 x_cord{sx_c0.value(), sx_c1.value()};
    if (infinity_flat_set & compression_flag_set & (!sort_flat_set) & x_cord.is_zero()) return G2Affine::identity();

    std::optional<field::Fp2> sy = ((x_cord.square() * x_cord) + field::constant::B2).sqrt();
    if (!sy.has_value()) return std::nullopt;

    field::Fp2 y_cord = sy.value();
    if (y_cord.lexicographically_largest() ^ sort_flat_set) y_cord = -y_cord;

    if (infinity_flat_set || (!compression_flag_set)) return std::nullopt;
    return G2Affine{x_cord, y_cord, infinity_flat_set};
}

std::optional<G2Affine> G2Affine::from_uncompressed(const std::array<uint8_t, G2Affine::BYTE_SIZE * 2> &bytes) {
    std::optional<G2Affine> res = G2Affine::from_uncompressed_unchecked(bytes);
    if (!res.has_value()) return std::nullopt;
    if (!res.value().is_on_curve() || !res.value().is_torsion_free()) return std::nullopt;
    return res.value();
}

std::optional<G2Affine>
G2Affine::from_uncompressed_unchecked(const std::array<uint8_t, G2Affine::BYTE_SIZE * 2> &bytes) {
    bool compression_flag_set = (bytes[0] >> 7) & 1;
    bool infinity_flat_set = (bytes[0] >> 6) & 1;
    bool sort_flat_set = (bytes[0] >> 5) & 1;

    // try to decode the x-coordinate
    std::array<uint8_t, field::Fp::BYTE_SIZE> temp{0};
    std::copy(bytes.begin(), bytes.begin() + field::Fp::BYTE_SIZE, temp.begin());
    temp[0] &= 0b00011111;
    std::optional<field::Fp> sx_c1 = field::Fp::from_bytes(temp);
    std::copy(bytes.begin() + field::Fp::BYTE_SIZE, bytes.begin() + G2Affine::BYTE_SIZE, temp.begin());
    std::optional<field::Fp> sx_c0 = field::Fp::from_bytes(temp);

    // try to decode the y-coordinate
    std::copy(bytes.begin() + G2Affine::BYTE_SIZE, bytes.begin() + field::Fp::BYTE_SIZE + G2Affine::BYTE_SIZE,
              temp.begin());
    std::optional<field::Fp> sy_c1 = field::Fp::from_bytes(temp);
    std::copy(bytes.begin() + field::Fp::BYTE_SIZE + G2Affine::BYTE_SIZE, bytes.end(), temp.begin());
    std::optional<field::Fp> sy_c0 = field::Fp::from_bytes(temp);

    if (!sx_c0.has_value() || !sx_c1.has_value() || !sy_c0.has_value() || !sy_c1.has_value()) return std::nullopt;
    field::Fp2 x_cord{sx_c0.value(), sx_c1.value()};
    field::Fp2 y_cord{sy_c0.value(), sy_c1.value()};

    if (compression_flag_set || sort_flat_set) return std::nullopt;
    if (infinity_flat_set && (!x_cord.is_zero() || !y_cord.is_zero())) return std::nullopt;

    G2Affine res{x_cord, y_cord, infinity_flat_set};
    if (infinity_flat_set) {
        return G2Affine::identity();
    } else {
        return res;
    }
}

field::Fp2 G2Affine::get_x() const {
    return this->x;
}

field::Fp2 G2Affine::get_y() const {
    return this->y;
}

bool G2Affine::is_identity() const {
    return this->infinity;
}

bool G2Affine::is_on_curve() const {
    return (this->y.square() - (this->x.square() * this->x)) == field::constant::B2 | this->infinity;
}

bool G2Affine::is_torsion_free() const {
    G2Projective p(*this);
    return p.psi() == p.mul_by_x();
}

std::array<uint8_t, G2Affine::BYTE_SIZE> G2Affine::to_compressed() const {
    field::Fp2 x_cord = (this->infinity ? field::Fp2::zero() : this->x);
    std::array<uint8_t, field::Fp::BYTE_SIZE> temp_x0 = x_cord.get_c0().to_bytes();
    std::array<uint8_t, field::Fp::BYTE_SIZE> temp_x1 = x_cord.get_c1().to_bytes();

    std::array<uint8_t, G2Affine::BYTE_SIZE> bytes{};
    std::copy(temp_x1.begin(), temp_x1.end(), bytes.begin());
    std::copy(temp_x0.begin(), temp_x0.end(), bytes.begin() + field::Fp::BYTE_SIZE);

    bytes[0] |= (static_cast<uint8_t>(1) << 7); // compression flag
    bytes[0] |= (this->infinity ? (static_cast<uint8_t>(1) << 6) : static_cast<uint8_t>(0)); // infinity flag
    bytes[0] |= (((!this->infinity) && this->y.lexicographically_largest()) ? (static_cast<uint8_t>(1) << 5)
                                                                            : static_cast<uint8_t>(0));// sort flag

    return bytes;
}

std::array<uint8_t, G2Affine::BYTE_SIZE * 2> G2Affine::to_uncompressed() const {
    field::Fp2 x_cord = (this->infinity ? field::Fp2::zero() : this->x);
    field::Fp2 y_cord = (this->infinity ? field::Fp2::zero() : this->y);

    std::array<uint8_t, field::Fp::BYTE_SIZE> temp_x0 = x_cord.get_c0().to_bytes();
    std::array<uint8_t, field::Fp::BYTE_SIZE> temp_x1 = x_cord.get_c1().to_bytes();
    std::array<uint8_t, field::Fp::BYTE_SIZE> temp_y0 = y_cord.get_c0().to_bytes();
    std::array<uint8_t, field::Fp::BYTE_SIZE> temp_y1 = y_cord.get_c1().to_bytes();

    std::array<uint8_t, G2Affine::BYTE_SIZE * 2> bytes{};

    std::copy(temp_x1.begin(), temp_x1.end(), bytes.begin());
    std::copy(temp_x0.begin(), temp_x0.end(), bytes.begin() + field::Fp::BYTE_SIZE);
    std::copy(temp_y1.begin(), temp_y1.end(), bytes.begin() + G2Affine::BYTE_SIZE);
    std::copy(temp_y0.begin(), temp_y0.end(), bytes.begin() + G2Affine::BYTE_SIZE + field::Fp::BYTE_SIZE);

    bytes[0] |= (this->infinity ? (static_cast<uint8_t>(1) << 6) : static_cast<uint8_t>(0)); // infinity flag

    return bytes;
}

G2Affine G2Affine::operator-() const {
    return G2Affine{
            this->x,
            this->infinity ? field::Fp2::one() : (-this->y),
            this->infinity,
    };
}

G2Affine &G2Affine::operator=(const G2Affine &rhs) {
    if (this == &rhs) return *this;
    this->x = rhs.x;
    this->y = rhs.y;
    this->infinity = rhs.infinity;
    return *this;
}

G2Projective operator+(const G2Affine &a, const G2Projective &b) {
    return G2Projective(b) += a;
}

G2Projective operator-(const G2Affine &a, const G2Projective &b) {
    return -G2Projective(b) += a;
}

G2Projective operator*(const G2Affine &a, const scalar::Scalar &b) {
    return G2Projective(a) *= b;
}

} // namespace bls12_381::group