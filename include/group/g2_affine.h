#ifndef BLS12_381_G2_AFFINE_H
#define BLS12_381_G2_AFFINE_H

#include "field/fp.h"

class G2Affine {
    Fp x;
    Fp y;
    bool infinity;
};

#endif //BLS12_381_G2_AFFINE_H
