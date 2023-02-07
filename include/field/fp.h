#ifndef BLS12_381_FP_H
#define BLS12_381_FP_H

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace bls12_381::field {

class Fp {
public:
    static constexpr int32_t WIDTH = 6;
    static constexpr int32_t BYTE_SIZE = Fp::WIDTH * sizeof(uint64_t);

private:
    std::array<uint64_t, Fp::WIDTH> data;

public:
    Fp();
    Fp(const Fp &fp);
    Fp(Fp &&fp) noexcept;

    explicit Fp(uint64_t val);
    explicit Fp(std::array<uint64_t, Fp::WIDTH> &&data);
    explicit Fp(const std::array<uint64_t, Fp::WIDTH> &data);

    static Fp zero() noexcept;
    static Fp one() noexcept;
    static Fp random() noexcept;

    static Fp montgomery_reduce(const std::array<uint64_t, Fp::WIDTH * 2> &ts);
    static Fp sum_of_products(const std::vector<Fp> &a, const std::vector<Fp> &b);
    static std::optional<Fp> from_bytes(const std::array<uint8_t, Fp::BYTE_SIZE> &bytes);

    [[nodiscard]] bool is_zero() const;
    [[nodiscard]] bool lexicographically_largest() const;

    [[nodiscard]] std::string getHex() const;
    [[nodiscard]] std::array<uint8_t, Fp::BYTE_SIZE> to_bytes() const;

    [[nodiscard]] Fp square() const;
    [[nodiscard]] Fp subtract_modulus() const;
    [[nodiscard]] Fp pow_vartime(const std::array<uint64_t, Fp::WIDTH> &exp) const;

    [[nodiscard]] std::optional<Fp> sqrt() const;
    [[nodiscard]] std::optional<Fp> invert() const;

private:
    static Fp reduce(const std::array<uint64_t, Fp::WIDTH * 2> &limbs);

public:
    Fp operator-() const;
    Fp &operator=(Fp &&rhs) noexcept;
    Fp &operator=(const Fp &rhs);
    Fp &operator+=(const Fp &rhs);
    Fp &operator-=(const Fp &rhs);
    Fp &operator*=(const Fp &rhs);

public:
    friend inline Fp operator+(const Fp &a, const Fp &b) { return Fp(a) += b; }
    friend inline Fp operator-(const Fp &a, const Fp &b) { return Fp(a) -= b; }
    friend inline Fp operator*(const Fp &a, const Fp &b) { return Fp(a) *= b; }

    friend inline bool operator==(const Fp &a, const Fp &b) { return a.data == b.data; }
    friend inline bool operator!=(const Fp &a, const Fp &b) { return a.data != b.data; }
};

} // namespace bls12_381::field

#endif //BLS12_381_FP_H