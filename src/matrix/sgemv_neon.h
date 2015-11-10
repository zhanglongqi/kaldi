/* 
 * File:   sgemv_neon.h
 * Author: longqi
 *
 * Created on November 9, 2015, 4:21 PM
 */

#ifndef SGEMV_NEON_H
#define	SGEMV_NEON_H

//#define _NEON_
#define _ACCELERATE_
#include <cstdlib>
#include "f2c.h"
#include "blaswrap.h"
logical lsame_(char *ca, char *cb);

int accelerate_sgemv(char trans, int m, int n,
		float alpha, const float *a, int lda, const float *x, int incx, float beta,
		float *y, int incy);

#endif	/* SGEMV_NEON_H */

