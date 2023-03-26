#include "group/g2_affine.h"

#include "utils/bit.h"

#include "field/constant.h"
#include "group/g2_projective.h"
#include "scalar/scalar.h"

namespace bls12_381::group {

using rng::util::bit::to_le_bytes;
using rng::util::bit::from_le_bytes;

using field::Fp;
using field::Fp2;
using scalar::Scalar;

G2Affine::G2Affine() : x{Fp2::zero()}, y{Fp2::one()}, infinity{true} {}

G2Affine::G2Affine(const G2Affine &point) = default;

G2Affine::G2Affine(const G2Projective &point) : x{Fp2::zero()}, y{Fp2::one()}, infinity{true} {
    Fp2 z_inv = point.get_z().invert().value_or(Fp2::zero());
    Fp2 rx = point.get_x() * z_inv;
    Fp2 ry = point.get_y() * z_inv;
    G2Affine temp{rx, ry, false};
    if (!z_inv.is_zero()) *this = temp;
}

G2Affine::G2Affine(const Fp2 &x, const Fp2 &y, bool infinity) : x{x}, y{y}, infinity{infinity} {}

G2Affine::G2Affine(G2Affine &&point) noexcept = default;

G2Affine::G2Affine(G2Projective &&point) : x{Fp2::zero()}, y{Fp2::one()}, infinity{true} {
    const Fp2 z_inv = point.get_z().invert().value_or(Fp2::zero());
    const Fp2 rx = point.get_x() * z_inv;
    const Fp2 ry = point.get_y() * z_inv;
    G2Affine temp{rx, ry, false};
    if (!z_inv.is_zero()) *this = temp;
}

G2Affine::G2Affine(Fp2 &&x, Fp2 &&y, bool infinity) : x{std::move(x)}, y{std::move(y)}, infinity{infinity} {}

G2Affine G2Affine::identity() noexcept {
    return G2Affine{};
}

G2Affine G2Affine::generator() noexcept {
    return G2Affine{
            Fp2{
                    Fp({
                               0xf5f28fa202940a10, 0xb3f5fb2687b4961a, 0xa1a893b53e2ae580,
                               0x9894999d1a3caee9, 0x6f67b7631863366b, 0x058191924350bcd7,
                       }),
                    Fp({
                               0xa5a9c0759e23f606, 0xaaa0c59dbccd60c3, 0x3bb17e18e2867806,
                               0x1b1ab6cc8541b367, 0xc2b6ed0ef2158547, 0x11922a097360edf3,
                       }),
            },
            Fp2{
                    Fp({
                               0x4c730af860494c4a, 0x597cfa1f5e369c5a, 0xe7e6856caa0a635a,
                               0xbbefb5e96e0d495f, 0x07d3a975f0ef25a2, 0x0083fd8e7e80dae5,
                       }),
                    Fp({
                               0xadc0fc92df64b05d, 0x18aa270a2b1461dc, 0x86adac6a3be4eba0,
                               0x79495c4ec93da33a, 0xe7175850a43ccaed, 0x0b2bc2a163de1bf2,
                       }),
            },
            false
    };
}

std::optional<G2Affine> G2Affine::from_compressed(const std::array<uint8_t, G2Affine::BYTE_SIZE> &bytes) {
    const std::optional<G2Affine> res = G2Affine::from_compressed_unchecked(bytes);
    if (!res.has_value()) return std::nullopt;
    if (!res.value().is_torsion_free()) return std::nullopt;
    return res.value();
}

std::optional<G2Affine> G2Affine::from_compressed_unchecked(const std::array<uint8_t, G2Affine::BYTE_SIZE> &bytes) {
    const bool compression_flag_set = (bytes[0] >> 7) & 1;
    const bool infinity_flat_set = (bytes[0] >> 6) & 1;
    const bool sort_flat_set = (bytes[0] >> 5) & 1;

    // try to decode the x-coordinate
    std::array<uint8_t, Fp::BYTE_SIZE> temp{0};
    std::copy(bytes.begin(), bytes.begin() + Fp::BYTE_SIZE, temp.begin());
    temp[0] &= 0b00011111;
    const std::optional<Fp> sx_c1 = Fp::from_bytes(temp);

    std::copy(bytes.begin() + Fp::BYTE_SIZE, bytes.end(), temp.begin());
    const std::optional<Fp> sx_c0 = Fp::from_bytes(temp);

    if (!sx_c0.has_value() || !sx_c1.has_value()) return std::nullopt;
    const Fp2 x_cord{sx_c0.value(), sx_c1.value()};
    if (infinity_flat_set & compression_flag_set & (!sort_flat_set) & x_cord.is_zero()) return G2Affine::identity();

    const std::optional<Fp2> sy = ((x_cord.square() * x_cord) + field::constant::B2).sqrt();
    if (!sy.has_value()) return std::nullopt;

    Fp2 y_cord = sy.value();
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
    const bool compression_flag_set = (bytes[0] >> 7) & 1;
    const bool infinity_flat_set = (bytes[0] >> 6) & 1;
    const bool sort_flat_set = (bytes[0] >> 5) & 1;

    // try to decode the x-coordinate
    std::array<uint8_t, Fp::BYTE_SIZE> temp{0};
    std::copy(bytes.begin(), bytes.begin() + Fp::BYTE_SIZE, temp.begin());
    temp[0] &= 0b00011111;
    const std::optional<Fp> sx_c1 = Fp::from_bytes(temp);
    std::copy(bytes.begin() + Fp::BYTE_SIZE, bytes.begin() + G2Affine::BYTE_SIZE, temp.begin());
    const std::optional<Fp> sx_c0 = Fp::from_bytes(temp);

    // try to decode the y-coordinate
    std::copy(bytes.begin() + G2Affine::BYTE_SIZE, bytes.begin() + Fp::BYTE_SIZE + G2Affine::BYTE_SIZE,
              temp.begin());
    const std::optional<Fp> sy_c1 = Fp::from_bytes(temp);
    std::copy(bytes.begin() + Fp::BYTE_SIZE + G2Affine::BYTE_SIZE, bytes.end(), temp.begin());
    const std::optional<Fp> sy_c0 = Fp::from_bytes(temp);

    if (!sx_c0.has_value() || !sx_c1.has_value() || !sy_c0.has_value() || !sy_c1.has_value()) return std::nullopt;
    const Fp2 x_cord{sx_c0.value(), sx_c1.value()};
    const Fp2 y_cord{sy_c0.value(), sy_c1.value()};

    if (compression_flag_set || sort_flat_set) return std::nullopt;
    if (infinity_flat_set && (!x_cord.is_zero() || !y_cord.is_zero())) return std::nullopt;

    G2Affine res{x_cord, y_cord, infinity_flat_set};
    if (infinity_flat_set) {
        return G2Affine::identity();
    } else {
        return res;
    }
}

const Fp2 &G2Affine::get_x() const noexcept {
    return this->x;
}

const Fp2 &G2Affine::get_y() const noexcept {
    return this->y;
}

bool G2Affine::is_identity() const {
    return this->infinity;
}

bool G2Affine::is_on_curve() const {
    return (this->y.square() - (this->x.square() * this->x)) == field::constant::B2 | this->infinity;
}

bool G2Affine::is_torsion_free() const {
    const G2Projective p(*this);
    return p.psi() == p.mul_by_x();
}

std::array<uint8_t, G2Affine::BYTE_SIZE> G2Affine::to_compressed() const {
    const Fp2 x_cord = (this->infinity ? Fp2::zero() : this->x);
    const std::array<uint8_t, Fp::BYTE_SIZE> temp_x0 = x_cord.get_c0().to_bytes();
    const std::array<uint8_t, Fp::BYTE_SIZE> temp_x1 = x_cord.get_c1().to_bytes();

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
    const Fp2 x_cord = (this->infinity ? Fp2::zero() : this->x);
    const Fp2 y_cord = (this->infinity ? Fp2::zero() : this->y);

    const std::array<uint8_t, Fp::BYTE_SIZE> temp_x0 = x_cord.get_c0().to_bytes();
    const std::array<uint8_t, Fp::BYTE_SIZE> temp_x1 = x_cord.get_c1().to_bytes();
    const std::array<uint8_t, Fp::BYTE_SIZE> temp_y0 = y_cord.get_c0().to_bytes();
    const std::array<uint8_t, Fp::BYTE_SIZE> temp_y1 = y_cord.get_c1().to_bytes();

    std::array<uint8_t, G2Affine::BYTE_SIZE * 2> bytes{};

    std::copy(temp_x1.begin(), temp_x1.end(), bytes.begin());
    std::copy(temp_x0.begin(), temp_x0.end(), bytes.begin() + Fp::BYTE_SIZE);
    std::copy(temp_y1.begin(), temp_y1.end(), bytes.begin() + G2Affine::BYTE_SIZE);
    std::copy(temp_y0.begin(), temp_y0.end(), bytes.begin() + G2Affine::BYTE_SIZE + Fp::BYTE_SIZE);

    bytes[0] |= (this->infinity ? (static_cast<uint8_t>(1) << 6) : static_cast<uint8_t>(0)); // infinity flag

    return bytes;
}

G2Affine G2Affine::operator-() const {
    return G2Affine{
            this->x,
            this->infinity ? Fp2::one() : (-this->y),
            this->infinity,
    };
}

G2Affine &G2Affine::operator=(const G2Affine &rhs) = default;

G2Affine &G2Affine::operator=(G2Affine &&rhs) noexcept = default;

G2Projective operator+(const G2Affine &a, const G2Projective &b) {
    return G2Projective(b) += a;
}

G2Projective operator-(const G2Affine &a, const G2Projective &b) {
    return -G2Projective(b) += a;
}

G2Projective operator*(const G2Affine &a, const Scalar &b) {
    return G2Projective(a) *= b;
}

auto G2Affine::from_slice_unchecked(const std::vector<uint8_t> &bytes) -> G2Affine {
    std::array<uint64_t, Fp::WIDTH> x0_bytes{};
    std::array<uint64_t, Fp::WIDTH> x1_bytes{};
    std::array<uint64_t, Fp::WIDTH> y0_bytes{};
    std::array<uint64_t, Fp::WIDTH> y1_bytes{};

    std::array<uint8_t, sizeof(uint64_t)> temp_bytes{};

for (int i = 0; i < G2Affine::RAW_SIZE - 1; i += 8) {
        const int count = i / 8;
        std::copy(bytes.begin() + i, bytes.begin() + i + 8, temp_bytes.begin());
        if (count < Fp::WIDTH)
            x0_bytes[count] = from_le_bytes<uint64_t>(temp_bytes);
        else if (count < Fp::WIDTH * 2)
            x1_bytes[count - Fp::WIDTH] = from_le_bytes<uint64_t>(temp_bytes);
        else if (count < Fp::WIDTH * 3)
            y0_bytes[count - Fp::WIDTH * 2] = from_le_bytes<uint64_t>(temp_bytes);
        else
            y1_bytes[count - Fp::WIDTH * 3] = from_le_bytes<uint64_t>(temp_bytes);
    }

    const Fp2 x = Fp2{Fp{x0_bytes}, Fp{x1_bytes}};
    const Fp2 y = Fp2{Fp{y0_bytes}, Fp{y1_bytes}};
    bool infinity = false;

    if (bytes.size() >= G2Affine::RAW_SIZE)
        infinity = static_cast<bool>(bytes[G2Affine::RAW_SIZE - 1]);

    return G2Affine{x, y, infinity};
}

std::array<uint8_t, G2Affine::RAW_SIZE> G2Affine::to_raw_bytes() const {
    std::array<uint8_t, G2Affine::RAW_SIZE> bytes{};
    for (int i = 0; i < Fp::WIDTH; ++i) {
        const std::array<uint8_t, 8> x0_bytes = to_le_bytes<uint64_t>(this->x.get_c0().get_data()[i]);
        const std::array<uint8_t, 8> x1_bytes = to_le_bytes<uint64_t>(this->x.get_c1().get_data()[i]);
        const std::array<uint8_t, 8> y0_bytes = to_le_bytes<uint64_t>(this->y.get_c0().get_data()[i]);
        const std::array<uint8_t, 8> y1_bytes = to_le_bytes<uint64_t>(this->y.get_c1().get_data()[i]);

        std::copy(x0_bytes.begin(), x0_bytes.end(), bytes.begin() + i * 8);
        std::copy(x1_bytes.begin(), x1_bytes.end(), bytes.begin() + i * 8 + Fp::BYTE_SIZE);
        std::copy(y0_bytes.begin(), y0_bytes.end(), bytes.begin() + i * 8 + Fp::BYTE_SIZE * 2);
        std::copy(y1_bytes.begin(), y1_bytes.end(), bytes.begin() + i * 8 + Fp::BYTE_SIZE * 3);
    }
    bytes[G2Affine::RAW_SIZE - 1] = static_cast<uint8_t>(this->infinity);
    return bytes;
}

} // namespace bls12_381::group