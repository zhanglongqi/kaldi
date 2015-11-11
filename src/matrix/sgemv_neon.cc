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
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    int i__, j, ix, iy, jx, jy, kx, ky, info;
    float temp;
    int lenx, leny;

    //#define x(I) X[(I)-1]
    //#define y(I) Y[(I)-1]
    //#define a(I) A[(I)-1]

    //	extern logical lsame_(char *, char *);
    //	extern /* Subroutine */int xerbla_(char *, int *);

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
    //	a_offset = 1 + a_dim1;
    //	a -= a_offset;
    //	--x;
    //	--y;
    /* Function Body */
    /*	info = 0;
     if (!lsame_(trans, "N") && !lsame_(trans, "T") && !lsame_(trans, "C"))
     {
     info = 1;
     }
     else if (m< 0)
     {
     info = 2;
     }
     else if (n < 0)
     {
     info = 3;
     }
     else if (*lda < max(1,m))
     {
     info = 6;
     }
     else if (incx == 0)
     {
     info = 8;
     }
     else if (incy == 0)
     {
     info = 11;
     }
     if (info != 0)
     {
     xerbla_("SGEMV ", &info);
     return 0;
     }
     */
    /*     Quick return if possible. */

    if (m == 0 || n == 0 || alpha == 0.f && beta == 1.f) {
        return 0;
    }

    /*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
    /*     up the start points in  X  and  Y. */

    if (trans == 'o') {
        lenx = n;
        leny = m;
    } else {
        lenx = m;
        leny = n;
    }
    if (incx > 0) {
        kx = 1;
    } else {
        kx = 1 - (lenx - 1) * incx; //Charles: kx >= 1
    }
    if (incy > 0) {
        ky = 1;
    } else {
        ky = 1 - (leny - 1) * incy; //Charles: ky >= 1
    }
#ifdef _NEON_
    int row, column;
    float32x4_t v4result, v4A, v4X, valpha, vbeta, v4Y, v4AX;

    valpha = vdupq_n_f32(alpha); //load same literal value
    vbeta = vdupq_n_f32(beta); //load same literal value

#endif

    /*     Start the operations. In this version the elements of A are */
    /*     accessed sequentially with one pass through A. */

    /*     First form  y := beta*y. */

    if (beta != 1.f) {
        if (incy == 1) {
            if (beta == 0.f) //Charles: if beta is 0, then clr all data to 0
            {

                i__1 = leny;
                for (i__ = 0; i__ < i__1; ++i__) {
                    y[i__] = 0.f;
                    /* L10: */
                }
            } else //Charles: real operation
            {

                //#ifdef _NEON_
                //				i__1 = leny;
                //				if (i__1 < 4)
                //				{
                //					for (i__ = 0; i__ < i__1; ++i__)
                //					{
                //						y[i__] = beta * y[i__];
                //						/* L20: */
                //					}
                //				}
                //				else //Charles: using neon
                //				{
                //					i1 = leny / 4;
                //					i2 = leny % 4;
                //					y++; //for compensating the --y in the Parameter adjustments
                //					//calc 4 data
                //					for (i__ = 1; i__ <= i1; ++i__)
                //					{
                //						vy = vld1q_f32(y); //load 4 value from y
                //						v4result = vmulq_f32(vbeta, vy); //result = b*y
                //						vst1q_f32(y, v4result); //store 4 values back to vy
                //						y += 4;
                //					}
                //					//calc remain
                //					for (i__ = 1; i__ <= i2; ++i__)
                //					{
                //						*y = beta * (*y);
                //						y++;
                //						/* L20: */
                //					}
                //				}
                //#else
                //				i__1 = leny;
                //				for (i__ = 0; i__ < i__1; ++i__)
                //				{
                //					y[i__] = beta * y[i__];
                //					/* L20: */
                //				}
                //#endif
            }
        } else //for incy > 1
        {
            iy = ky;
            if (beta == 0.f) {
                i__1 = leny;
                for (i__ = 0; i__ < i__1; ++i__) {
                    y[iy] = 0.f;
                    iy += incy;
                    /* L30: */
                }
            } else {
                i__1 = leny;
                for (i__ = 0; i__ < i__1; ++i__) {
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
    if (trans == 'o') {
       
        /*        Form  y := alpha*A*x + y. */

        jx = kx;
        if (incy == 1) {
#ifdef _NEON_
            float32x4_t temp_AX, v4bY;
            float temp_4AX[4];
            float *temp_final_AX;
            float *tx, *ta;
            temp_final_AX = (float*) malloc(m * sizeof (float));

            int i1, i2;
            i1 = n >> 2;
            i2 = n % 4;

            //calc AX
            ta = (float*) a;
            for (row = 0; row < m; ++row) {
                //calc 4 multiples
                v4AX = vdupq_n_f32(0); //initial temp_AX

                int icolumn;
                tx = (float*) x;

                temp_final_AX[row] = 0;
                //calc 4 multiples
                for (icolumn = 0; icolumn < i1; ++icolumn) {

                    v4A = vld1q_f32(ta); //load 4 values
                    v4X = vld1q_f32(tx);
                    v4AX = vmlaq_f32(v4AX, v4A, v4X); //result = a + b*c
                    vst1q_f32(temp_4AX, v4AX);
                    temp_final_AX[row] = temp_4AX[0] + temp_4AX[1] + temp_4AX[2]
                            + temp_4AX[3];
                    ta += 4;
                    tx += 4;
                }
                //calc remain
                for (icolumn = 0; icolumn < i2; ++icolumn) {
                    temp_final_AX[row] += (*ta) * (*tx);
                    ta++;
                    tx++;
                }
            }

            //calc y = aAX + bY;
            //y is m vector
            i1 = m >> 2;
            i2 = m % 4;
            //calc 4 multiples
            for (row = 0; row < i1; ++row) {
                v4Y = vld1q_f32(y); //load 4 Y
                v4AX = vld1q_f32(temp_final_AX); //load 4 AX
                v4bY = vmulq_f32(vbeta, v4Y); // b*Y
                v4result = vmlaq_f32(v4bY, valpha, v4AX); //result = a + b*c
                vst1q_f32(y, v4result); //store 4 values back to y
                y += 4;
                temp_final_AX += 4;
            }
            //calc remain
            for (row = 0; row < i2; ++row) {
                (*y) = alpha * (*temp_final_AX) + beta * (*y);
                y++;
                temp_final_AX++;
            }
#else
            i__1 = n;
            for (j = 0; j < i__1; ++j) {
                if (x[jx] != 0.f) {
                    temp = alpha * x[jx];
                    i__2 = m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        y[i__] += temp * a[i__ + j * a_dim1];
                        /* L50: */
                    }
                }
                jx += incx;
                /* L60: */
            }
#endif
        } else //for incx > 1
        {
            i__1 = n;
            for (j = 0; j < i__1; ++j) {
                if (x[jx] != 0.f) {
                    temp = alpha * x[jx];
                    iy = ky;
                    i__2 = m;
                    for (i__ = 0; i__ < i__2; ++i__) {
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
            for (j = 0; j < i__1; ++j) {
                temp = 0.f;
                i__2 = m;
                for (i__ = 0; i__ < i__2; ++i__) {
                    temp += a[i__ + j * a_dim1] * x[i__];
                    /* L90: */
                }
                y[jy] += alpha * temp;
                jy += incy;
                /* L100: */
            }
        } else {
            i__1 = n;
            for (j = 0; j < i__1; ++j) {
                temp = 0.f;
                ix = kx;
                i__2 = m;
                for (i__ = 0; i__ < i__2; ++i__) {
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

} /* sgemv_ */

