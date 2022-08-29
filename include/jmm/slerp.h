#pragma once

#include "def.h"

void nlerp3(dbl const p[3][3], dbl const w[3], dbl q[3]);

void slerp2(dbl const p0[3], dbl const p1[3], dbl const w[2], dbl q[3]);
void slerp3(dbl const p[3][3], dbl const w[3], dbl q[3], dbl tol);
void slerp4(dbl const p[4][3], dbl const w[4], dbl q[3], dbl tol);
