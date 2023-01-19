#ifndef BLS12_381_G1_H
#define BLS12_381_G1_H

#include "field/fp.h"

class G1Affine {
    Fp x;
    Fp y;
    bool infinity;
};

#endif //BLS12_381_G1_H
