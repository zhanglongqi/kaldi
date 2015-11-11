/*
 * sgemv.cc
 *
 *  Created on: 2015-11-9
 *      Author: CharlesZhang
 */
#include "sgemv_neon.h"
#ifdef _NEON_

#include <arm_neon.h>

#endif

/* Subroutine */
int accelerate_sgemv(char trans, int m, int n, float alpha,
        const float *a, int lda, const float *x, int incx, float beta, float *y,
        int incy) {
#ifdef _NEON_
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    int i__, j, ix, iy, jx, jy, kx, ky, info;
    float temp;
    int lenx, leny;
    //    extern logical lsame_(char *, char *);
    //    extern /* Subroutine */ int xerbla_(char *, int *);

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*  SGEMV  performs one of the matrix-vector operations */

    /*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y, */

    /*  where alpha and beta are scalars, x and y are vectors and A is an */
    /*  m by n matrix. */

    /*  Arguments */
    /*  ========== */

    /*  TRANS  - CHARACTER*1. */
    /*           On entry, TRANS specifies the operation to be performed as */
    /*           follows: */

    /*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y. */

    /*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y. */

    /*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y. */

    /*           Unchanged on exit. */

    /*  M      - INTEGER. */
    /*           On entry, M specifies the number of rows of the matrix A. */
    /*           M must be at least zero. */
    /*           Unchanged on exit. */

    /*  N      - INTEGER. */
    /*           On entry, N specifies the number of columns of the matrix A. */
    /*           N must be at least zero. */
    /*           Unchanged on exit. */

    /*  ALPHA  - REAL            . */
    /*           On entry, ALPHA specifies the scalar alpha. */
    /*           Unchanged on exit. */

    /*  A      - REAL             array of DIMENSION ( LDA, n ). */
    /*           Before entry, the leading m by n part of the array A must */
    /*           contain the matrix of coefficients. */
    /*           Unchanged on exit. */

    /*  LDA    - INTEGER. */
    /*           On entry, LDA specifies the first dimension of A as declared */
    /*           in the calling (sub) program. LDA must be at least */
    /*           max( 1, m ). */
    /*           Unchanged on exit. */

    /*  X      - REAL             array of DIMENSION at least */
    /*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n' */
    /*           and at least */
    /*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise. */
    /*           Before entry, the incremented array X must contain the */
    /*           vector x. */
    /*           Unchanged on exit. */

    /*  INCX   - INTEGER. */
    /*           On entry, INCX specifies the increment for the elements of */
    /*           X. INCX must not be zero. */
    /*           Unchanged on exit. */

    /*  BETA   - REAL            . */
    /*           On entry, BETA specifies the scalar beta. When BETA is */
    /*           supplied as zero then Y need not be set on input. */
    /*           Unchanged on exit. */

    /*  Y      - REAL             array of DIMENSION at least */
    /*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n' */
    /*           and at least */
    /*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise. */
    /*           Before entry with BETA non-zero, the incremented array Y */
    /*           must contain the vector y. On exit, Y is overwritten by the */
    /*           updated vector y. */

    /*  INCY   - INTEGER. */
    /*           On entry, INCY specifies the increment for the elements of */
    /*           Y. INCY must not be zero. */
    /*           Unchanged on exit. */

    /*  Level 2 Blas routine. */

    /*  -- Written on 22-October-1986. */
    /*     Jack Dongarra, Argonne National Lab. */
    /*     Jeremy Du Croz, Nag Central Office. */
    /*     Sven Hammarling, Nag Central Office. */
    /*     Richard Hanson, Sandia National Labs. */

    /*     .. Parameters .. */
    /*     .. */
    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. External Functions .. */
    /*     .. */
    /*     .. External Subroutines .. */
    /*     .. */
    /*     .. Intrinsic Functions .. */
    /*     .. */

    /*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;
    --y;

    /* Function Body */
    /*info = 0;
     if (! lsame_(trans, "N") && ! lsame_(trans, "T") && ! lsame_(trans, "C")
     ) {
     info = 1;
     } else if (m < 0) {
     info = 2;
     } else if (n < 0) {
     info = 3;
     } else if (lda < max(1,m)) {
     info = 6;
     } else if (incx == 0) {
     info = 8;
     } else if (incy == 0) {
     info = 11;
     }
     if (info != 0) {
     xerbla_("SGEMV ", &info);
     return 0;
     }*/

    /*     Quick return if possible. */

    if (m == 0 || n == 0 || alpha == 0.f && beta == 1.f) {
        return 0;
    }

    /*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
    /*     up the start points in  X  and  Y. */

    //    if (lsame_(trans, "N")) {
    if (1) {
        lenx = n;
        leny = m;
    } else {
        lenx = m;
        leny = n;
    }
    if (incx > 0) {
        kx = 1;
    } else {
        kx = 1 - (lenx - 1) * incx;
    }
    if (incy > 0) {
        ky = 1;
    } else {
        ky = 1 - (leny - 1) * incy;
    }

    /*     Start the operations. In this version the elements of A are */
    /*     accessed sequentially with one pass through A. */


    float32x4_t v4result, v4A, v4X, v4Y, v4AX;


    //	float32x4_t temp_AX, v4bY;
    //	float temp_4AX[4];
    //	float *temp_final_AX;
    const float *ta;
    float *ty;
    //	temp_final_AX = (float *) malloc(m * sizeof(float));

    //	valpha, vbeta,
    //	valpha = vdupq_n_f32(alpha); //load same literal value
    //	vbeta = vdupq_n_f32(beta); //load same literal value

    int i1, i2, remain, do_intger;

    /*     First form  y := beta*y. */

    if (beta != 1.f) {
        if (incy == 1) {
            if (beta == 0.f) {
                i__1 = leny;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    y[i__] = 0.f;
                    /* L10: */
                }
            }
            //			else
            //			{
            //				//Charles: modified
            //				i__1 = leny;
            //#ifdef _NEON1_
            //				//calc bY
            //				remain = leny % 4;
            //				ty = y+1;
            //				for (i__ = 1; i__ <= i__1 - remain; i__ += 4)
            //				{
            //					//calc 4 multiples
            //
            //					v4Y = vld1q_f32(ty);//load 4 Y
            ////					v4Y[0] = *(ty);//load 4 Y
            ////					v4Y[1] = *(ty+1);//load 4 Y
            ////					v4Y[2] = *(ty+2);//load 4 Y
            ////					v4Y[3] = *(ty+3);//load 4 Y
            //
            //					v4bY = vmulq_f32(vbeta, v4Y);// b*Y
            //
            //					vst1q_f32(ty, v4bY);//store 4 values back to y
            //					ty += 4;
            //				}
            //				//calc remain
            //				for (i1 = 0; i1 < remain; ++i1)
            //				{
            //					(*ty) = beta * (*ty);
            //					ty ++;
            //				}
            //#else
            ////				for (i__ = 1; i__ <= i__1; ++i__)
            ////				{
            ////					y[i__] = beta * y[i__];
            ////					/* L20: */
            ////				}
            //#endif
            //			}
        } else {
            iy = ky;
            if (beta == 0.f) {
                i__1 = leny;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    y[iy] = 0.f;
                    iy += incy;
                    /* L30: */
                }
            } else {
                i__1 = leny;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    y[iy] = beta * y[iy];
                    iy += incy;
                    /* L40: */
                }
            }
        }
    }
    if (alpha == 0.f) {
        return 0;
    }
    //    if (lsame_(trans, "N")) {
    if (1) {

        /*        Form  y := alpha*A*x + y. */

        jx = kx;
        if (incy == 1) {
            //Charles: modified

            int row_a[4];
            do_intger = (leny >> 2) + 1;
            remain = leny % 4;
            ta = a;
            ty = y + 1;

            i__1 = leny;
            for (j = 1; j <= leny; j += 4) {
                //calc AX
                //clean 4AX
                v4AX = vdupq_n_f32(0.f);
                i__2 = lenx;

                row_a[0] = (j) * a_dim1;
                row_a[1] = (j + 1) * a_dim1;
                row_a[2] = (j + 2) * a_dim1;
                row_a[3] = (j + 3) * a_dim1;

                for (i__ = 1; i__ <= i__2; ++i__) {
                    v4A[0] = *(ta + i__ + row_a[0]); //row +0
                    v4A[1] = *(ta + i__ + row_a[1]); //row +1
                    v4A[2] = *(ta + i__ + row_a[2]); //row +2
                    v4A[3] = *(ta + i__ + row_a[3]); //row +3

                    //load 4 same x(n)
                    v4X = vdupq_n_f32(x[i__]);

                    //v4AX += v4A * v4X
                    v4AX = vmlaq_f32(v4AX, v4A, v4X);

                }

                //calc y = bY
                v4Y = vld1q_f32(ty); //load 4 Y

                //Vector multiply by scalar,
                /*float32x4_t vmulq_n_f32(float32x4_t a, float32_t b); // VMUL.F32 q0,q0,d0[0]*/
                v4Y = vmulq_n_f32(v4Y, beta);
                //				v4Y = vmulq_f32(vbeta, v4Y); //b*Y

                //calc y = aAX + y
                /*
                 * Vector multiply accumulate with scalar
                 * float32x4_t vmlaq_n_f32(float32x4_t a, float32x4_t b, float32_t c);    // VMLA.F32 q0, q0, d0[0]
                 */
                v4result = vmlaq_n_f32(v4Y, v4AX, alpha);
                //				v4result = vmlaq_f32(v4Y, valpha, v4AX);

                if (j == do_intger) //deal with the remaining
                {
                    for (i1 = 0; i1 < remain; ++i1) {
                        (*ty) = v4result[i1];
                        ty++;
                    }
                } else {
                    vst1q_f32(ty, v4result);
                    ty += 4;
                }
            }

            

        } else {
            i__1 = n;
            for (j = 1; j <= i__1; ++j) {
                if (x[jx] != 0.f) {
                    temp = alpha * x[jx];
                    iy = ky;
                    i__2 = m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        y[iy] += temp * a[i__ + j * a_dim1];
                        iy += incy;
                        /* L70: */
                    }
                }
                jx += incx;
                /* L80: */
            }
        }
    } else {

        /*        Form  y := alpha*A'*x + y. */

        jy = ky;
        if (incx == 1) {
            i__1 = n;
            for (j = 1; j <= i__1; ++j) {
                temp = 0.f;
                i__2 = m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    temp += a[i__ + j * a_dim1] * x[i__];
                    /* L90: */
                }
                y[jy] += alpha * temp;
                jy += incy;
                /* L100: */
            }
        } else {
            i__1 = n;
            for (j = 1; j <= i__1; ++j) {
                temp = 0.f;
                ix = kx;
                i__2 = m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    temp += a[i__ + j * a_dim1] * x[ix];
                    ix += incx;
                    /* L110: */
                }
                y[jy] += alpha * temp;
                jy += incy;
                /* L120: */
            }
        }
    }

    return 0;

    /*     End of SGEMV . */
#else
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    int i__, j, ix, iy, jx, jy, kx, ky, info;
    float temp;
    int lenx, leny;
    //    extern logical lsame_(char *, char *);
    //    extern /* Subroutine */ int xerbla_(char *, int *);

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*  SGEMV  performs one of the matrix-vector operations */

    /*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y, */

    /*  where alpha and beta are scalars, x and y are vectors and A is an */
    /*  m by n matrix. */

    /*  Arguments */
    /*  ========== */

    /*  TRANS  - CHARACTER*1. */
    /*           On entry, TRANS specifies the operation to be performed as */
    /*           follows: */

    /*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y. */

    /*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y. */

    /*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y. */

    /*           Unchanged on exit. */

    /*  M      - INTEGER. */
    /*           On entry, M specifies the number of rows of the matrix A. */
    /*           M must be at least zero. */
    /*           Unchanged on exit. */

    /*  N      - INTEGER. */
    /*           On entry, N specifies the number of columns of the matrix A. */
    /*           N must be at least zero. */
    /*           Unchanged on exit. */

    /*  ALPHA  - REAL            . */
    /*           On entry, ALPHA specifies the scalar alpha. */
    /*           Unchanged on exit. */

    /*  A      - REAL             array of DIMENSION ( LDA, n ). */
    /*           Before entry, the leading m by n part of the array A must */
    /*           contain the matrix of coefficients. */
    /*           Unchanged on exit. */

    /*  LDA    - INTEGER. */
    /*           On entry, LDA specifies the first dimension of A as declared */
    /*           in the calling (sub) program. LDA must be at least */
    /*           max( 1, m ). */
    /*           Unchanged on exit. */

    /*  X      - REAL             array of DIMENSION at least */
    /*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n' */
    /*           and at least */
    /*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise. */
    /*           Before entry, the incremented array X must contain the */
    /*           vector x. */
    /*           Unchanged on exit. */

    /*  INCX   - INTEGER. */
    /*           On entry, INCX specifies the increment for the elements of */
    /*           X. INCX must not be zero. */
    /*           Unchanged on exit. */

    /*  BETA   - REAL            . */
    /*           On entry, BETA specifies the scalar beta. When BETA is */
    /*           supplied as zero then Y need not be set on input. */
    /*           Unchanged on exit. */

    /*  Y      - REAL             array of DIMENSION at least */
    /*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n' */
    /*           and at least */
    /*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise. */
    /*           Before entry with BETA non-zero, the incremented array Y */
    /*           must contain the vector y. On exit, Y is overwritten by the */
    /*           updated vector y. */

    /*  INCY   - INTEGER. */
    /*           On entry, INCY specifies the increment for the elements of */
    /*           Y. INCY must not be zero. */
    /*           Unchanged on exit. */

    /*  Level 2 Blas routine. */

    /*  -- Written on 22-October-1986. */
    /*     Jack Dongarra, Argonne National Lab. */
    /*     Jeremy Du Croz, Nag Central Office. */
    /*     Sven Hammarling, Nag Central Office. */
    /*     Richard Hanson, Sandia National Labs. */

    /*     .. Parameters .. */
    /*     .. */
    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. External Functions .. */
    /*     .. */
    /*     .. External Subroutines .. */
    /*     .. */
    /*     .. Intrinsic Functions .. */
    /*     .. */

    /*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;
    --y;

    /* Function Body */
    /*info = 0;
     if (! lsame_(trans, "N") && ! lsame_(trans, "T") && ! lsame_(trans, "C")
     ) {
     info = 1;
     } else if (m < 0) {
     info = 2;
     } else if (n < 0) {
     info = 3;
     } else if (lda < max(1,m)) {
     info = 6;
     } else if (incx == 0) {
     info = 8;
     } else if (incy == 0) {
     info = 11;
     }
     if (info != 0) {
     xerbla_("SGEMV ", &info);
     return 0;
     }*/

    /*     Quick return if possible. */

    if (m == 0 || n == 0 || alpha == 0.f && beta == 1.f) {
        return 0;
    }

    /*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
    /*     up the start points in  X  and  Y. */

    //    if (lsame_(trans, "N")) {
    if (1) {
        lenx = n;
        leny = m;
    } else {
        lenx = m;
        leny = n;
    }
    if (incx > 0) {
        kx = 1;
    } else {
        kx = 1 - (lenx - 1) * incx;
    }
    if (incy > 0) {
        ky = 1;
    } else {
        ky = 1 - (leny - 1) * incy;
    }

    /*     Start the operations. In this version the elements of A are */
    /*     accessed sequentially with one pass through A. */

    /*     First form  y := beta*y. */

    if (beta != 1.f) {
        if (incy == 1) {
            if (beta == 0.f) {
                i__1 = leny;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    y[i__] = 0.f;
                    /* L10: */
                }
            } else {
                //Charles: modified
                i__1 = leny;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    y[i__] = beta * y[i__];
                    /* L20: */
                }
            }
        } else {
            iy = ky;
            if (beta == 0.f) {
                i__1 = leny;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    y[iy] = 0.f;
                    iy += incy;
                    /* L30: */
                }
            } else {
                i__1 = leny;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    y[iy] = beta * y[iy];
                    iy += incy;
                    /* L40: */
                }
            }
        }
    }
    if (alpha == 0.f) {
        return 0;
    }
    //    if (lsame_(trans, "N")) {
    if (1) {

        /*        Form  y := alpha*A*x + y. */

        jx = kx;
        if (incy == 1) {
            //Charles: modified
            i__1 = leny;
            for (j = 1; j <= i__1; ++j) {
                i__2 = lenx;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    temp = alpha * a[i__ + j * a_dim1];
                    y[j] += temp * x[i__];
                    /* L50: */
                }
                /* L60: */
            }
        } else {
            i__1 = n;
            for (j = 1; j <= i__1; ++j) {
                if (x[jx] != 0.f) {
                    temp = alpha * x[jx];
                    iy = ky;
                    i__2 = m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        y[iy] += temp * a[i__ + j * a_dim1];
                        iy += incy;
                        /* L70: */
                    }
                }
                jx += incx;
                /* L80: */
            }
        }
    } else {

        /*        Form  y := alpha*A'*x + y. */

        jy = ky;
        if (incx == 1) {
            i__1 = n;
            for (j = 1; j <= i__1; ++j) {
                temp = 0.f;
                i__2 = m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    temp += a[i__ + j * a_dim1] * x[i__];
                    /* L90: */
                }
                y[jy] += alpha * temp;
                jy += incy;
                /* L100: */
            }
        } else {
            i__1 = n;
            for (j = 1; j <= i__1; ++j) {
                temp = 0.f;
                ix = kx;
                i__2 = m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    temp += a[i__ + j * a_dim1] * x[ix];
                    ix += incx;
                    /* L110: */
                }
                y[jy] += alpha * temp;
                jy += incy;
                /* L120: */
            }
        }
    }

    return 0;

    /*     End of SGEMV . */
#endif
} /* sgemv_ */