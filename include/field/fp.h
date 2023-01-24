#ifndef BLS12_381_FP_H
#define BLS12_381_FP_H

#include <cstdint>
#include <vector>
#include <optional>
#include <span>

class Fp {
private:
    static constexpr int32_t WIDTH = 6;
    uint64_t data[WIDTH];

public:
    Fp();
    explicit Fp(uint64_t val);
    explicit Fp(const std::vector<uint64_t> &vch);

    static Fp zero();
    static Fp one();
    static Fp random();
    static Fp montgomery_reduce(std::span<uint64_t> ts);
    static std::optional<Fp> from_bytes(std::span<uint8_t> bytes);

    [[nodiscard]] bool is_zero() const;
    [[nodiscard]] bool lexicographically_largest() const;
    [[nodiscard]] uint8_t *to_bytes(std::span<uint8_t> bytes) const;
    [[nodiscard]] std::string getHex() const;

    [[nodiscard]] Fp square() const;
    [[nodiscard]] Fp subtract_modulus() const;
    [[nodiscard]] Fp pow_vartime(std::span<uint64_t> exp) const;

    [[nodiscard]] std::optional<Fp> sqrt() const;
    [[nodiscard]] std::optional<Fp> invert() const;

public:
    template<uint32_t N>
    Fp static sum_of_products(std::span<Fp> a, std::span<Fp> b);

public:
    Fp &operator=(const Fp &rhs);
    Fp &operator+=(const Fp &rhs);
    Fp &operator-=(const Fp &rhs);
    Fp &operator*=(const Fp &rhs);

    Fp operator-() const;

public:
    friend inline Fp operator+(const Fp &a, const Fp &b) { return Fp(a) += b; }
    friend inline Fp operator-(const Fp &a, const Fp &b) { return Fp(a) -= b; }
    friend inline Fp operator*(const Fp &a, const Fp &b) { return Fp(a) *= b; }

    friend inline bool operator==(const Fp &a, const Fp &b) { return std::memcmp(a.data, b.data, sizeof(a.data)) == 0; }
    friend inline bool operator!=(const Fp &a, const Fp &b) { return std::memcmp(a.data, b.data, sizeof(a.data)) != 0; }

private:
    static Fp reduce(std::span<uint64_t> limbs);
};

#endif //BLS12_381_FP_H
