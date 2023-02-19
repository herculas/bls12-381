#include "pairing/miller_loop_driver.h"

namespace bls12_381::pairing {

using field::Fp2;
using group::G2Affine;
using group::G2Projective;

std::tuple<Fp2, Fp2, Fp2> doubling_step(G2Projective &point) {
    Fp2 tmp0 = point.get_x().square();
    Fp2 tmp1 = point.get_y().square();
    Fp2 tmp2 = tmp1.square();
    Fp2 tmp3 = (tmp1 + point.get_x()).square() - tmp0 - tmp2;
    tmp3 = tmp3 + tmp3;
    Fp2 tmp4 = tmp0 + tmp0 + tmp0;
    Fp2 tmp6 = point.get_x() + tmp4;
    Fp2 tmp5 = tmp4.square();
    Fp2 z_squared = point.get_z().square();

    Fp2 temp_x = tmp5 - tmp3 - tmp3;
    Fp2 temp_z = (point.get_z() + point.get_y()).square() - tmp1 - z_squared;
    Fp2 temp_y = (tmp3 - temp_x) * tmp4;

    tmp2 = tmp2 + tmp2;
    tmp2 = tmp2 + tmp2;
    tmp2 = tmp2 + tmp2;

    temp_y -= tmp2;

    tmp3 = tmp4 * z_squared;
    tmp3 = tmp3 + tmp3;
    tmp3 = -tmp3;
    tmp6 = tmp6.square() - tmp0 - tmp5;
    tmp1 = tmp1 + tmp1;
    tmp1 = tmp1 + tmp1;
    tmp6 = tmp6 - tmp1;

    tmp0 = temp_z * z_squared;

    tmp0 = tmp0 + tmp0;

    point = G2Projective{temp_x, temp_y, temp_z};
    return {tmp0, tmp3, tmp6};
}

std::tuple<Fp2, Fp2, Fp2> addition_step(G2Projective &r, const G2Affine &q) {
    Fp2 z_squared = r.get_z().square();
    Fp2 y_squared = q.get_y().square();
    Fp2 t0 = z_squared * q.get_x();
    Fp2 t1 = ((q.get_y() + r.get_z()).square() - y_squared - z_squared) * z_squared;
    Fp2 t2 = t0 - r.get_x();
    Fp2 t3 = t2.square();
    Fp2 t4 = t3 + t3;
    t4 = t4 + t4;
    Fp2 t5 = t4 * t2;
    Fp2 t6 = t1 - r.get_y() - r.get_y();
    Fp2 t9 = t6 * q.get_x();
    Fp2 t7 = t4 * r.get_x();

    Fp2 temp_r_x = t6.square() - t5 - t7 - t7;
    Fp2 temp_r_z = (r.get_z() + t2).square() - z_squared - t3;

    Fp2 t10 = q.get_y() + temp_r_z;
    Fp2 t8 = (t7 - temp_r_x) * t6;
    t0 = r.get_y() * t5;
    t0 = t0 + t0;

    Fp2 temp_r_y = t8 - t0;

    t10 = t10.square() - y_squared;
    Fp2 zt_squared = temp_r_z.square();
    t10 = t10 - zt_squared;
    t9 = t9 + t9 - t10;
    t10 = temp_r_z + temp_r_z;
    t6 = -t6;
    t1 = t6 + t6;

    r = G2Projective{temp_r_x, temp_r_y, temp_r_z};
    return {t10, t1, t9};
}

} // namespace bls12_381::pairing