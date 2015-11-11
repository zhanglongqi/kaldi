#include "f2c.h"
#define _NEON_
#define _ACCELERATE_
int accelerate_saxpy(const int n,const float sa,const float *sx,const int incx,
		float *sy,const int incy);