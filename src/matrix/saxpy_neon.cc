
/*
 * saxpy.c
 *
 *  Created on: 2015-11-8
 *      Author: CharlesZhang
 */
/*  -- translated by f2c (version 19940927).
 You must link the resulting object file with the libraries:
 -lf2c -lm   (in that order)
 */
#include "saxpy_neon.h"

#ifdef _NEON_

#include <arm_neon.h>

#endif

/* Subroutine */
int accelerate_saxpy(const int n, const float sa, const float *sx, const int incx,
        float *sy, const int incy) {

    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i, m, ix, iy, mp1;
#ifdef _NEON_
    float32x4_t vresult, vsx, vsy, vsa;
    int i1, i2;
    vsa = vdupq_n_f32(sa); //load same literal value to vsa
#endif
    /*     constant times a vector plus a vector.
     uses unrolled loop for increments equal to one.
     jack dongarra, linpack, 3/11/78.
     modified 12/3/93, array(1) declarations changed to array(*)



     Parameter adjustments
     Function Body */
#define SY(I) sy[(I)-1]
#define SX(I) sx[(I)-1]

    if (n <= 0) {
        return 0;
    }
    if (sa == 0.f) {
        return 0;
    }
    if (incx == 1 && incy == 1) {
        goto L20;
    }

    /*        code for unequal increments or equal increments
     not equal to 1 */

    ix = 1;
    iy = 1;
    if (incx < 0) {
        ix = (-(n) + 1) * incx + 1;
    }
    if (incy < 0) {
        iy = (-(n) + 1) * incy + 1;
    }
    i__1 = n;
    for (i = 1; i <= n; ++i) {
        SY(iy) += sa * SX(ix);
        ix += incx;
        iy += incy;
        /* L10: */
    }
    return 0;

    /*        code for both increments equal to 1


     clean-up loop */

L20:
#ifdef _NEON_
    i1 = n / 4;
    i2 = n % 4;
    for (i = 0; i < i1; ++i) {
        vsx = vld1q_f32(sx); //load 4 values to sx and sy
        vsy = vld1q_f32(sy);
        vresult = vmlaq_f32(vsy, vsa, vsx); //result = a + b*c
        vst1q_f32(sy, vresult); //store 4 values back to SY
        sx += 4;
        sy += 4;
    }
    //calc remain
    for (i = 0; i < i2; ++i) {
        (*sy) += sa * (*sx);
        sx++;
        sy++;
    }
#else
    m = n % 4;
    if (m == 0) {
        goto L40;
    }
    i__1 = m;
    for (i = 1; i <= m; ++i) {
        SY(i) += sa * SX(i);

        /* L30: */
    }
    if (n < 4) {
        return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = n;
    for (i = mp1; i <= n; i += 4) {
        SY(i) += sa * SX(i);
        SY(i + 1) += sa * SX(i + 1);
        SY(i + 2) += sa * SX(i + 2);
        SY(i + 3) += sa * SX(i + 3);
        /* L50: */
    }
#endif
    return 0;
} /* saxpy_ */