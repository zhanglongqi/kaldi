#include "f2c.h"
#include <arm_neon.h>

int accelerate_saxpy(const int n,const float sa,const float *sx,const int incx,
		float *sy,const int incy);