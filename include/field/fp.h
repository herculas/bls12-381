#ifndef BLS12_381_FP_H
#define BLS12_381_FP_H

#include <cstdint>
#include <vector>

class Fp {
private:
    uint64_t data[6]{};
public:
    Fp();
    explicit Fp(uint64_t val);
    explicit Fp(const std::vector<uint64_t> &vch);

    static Fp zero();
    static Fp one();

    bool is_zero();
};

#endif //BLS12_381_FP_H
