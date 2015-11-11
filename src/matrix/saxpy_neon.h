/* 
 * File:   saxpy_neon.h
 * Author: longqi
 *
 * Created on November 11, 2015, 9:24 PM
 */

#ifndef SAXPY_NEON_H
#define	SAXPY_NEON_H

#define _NEON_
#define _ACCELERATE_
#include <cstdlib>
#include "f2c.h"

int accelerate_saxpy(const int n, const float sa, const float *sx, const int incx,
        float *sy, const int incy);

#endif	/* SAXPY_NEON_H */

