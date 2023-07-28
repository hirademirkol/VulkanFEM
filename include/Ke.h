/*************************************************************************************

CUDAFEM

Author: Christian Dick

Copyright (c) Christian Dick

mailto:chr.dick@googlemail.com

This source code is property of Christian Dick. All rights reserved.
Unauthorized use prohibited.

*************************************************************************************/

#ifndef __TUM3D__KE_H__
#define __TUM3D__KE_H__

#include <math.h>
#include <memory.h>

template<typename scalar>
void ComputeKe(scalar a1, scalar a2, scalar a3, scalar E, scalar nu, scalar Ke[24][24]);

template<typename scalar>
void ComputeBavg(scalar a1, scalar a2, scalar a3, scalar Bavg[6][24]);

template<typename scalar>
void ComputeS(scalar a1, scalar a2, scalar a3, scalar E, scalar nu, scalar S[6][24]);


template void ComputeKe<double>(double a1, double a2, double a3, double E, double nu, double Ke[24][24]);
template void ComputeKe<float>(float a1, float a2, float a3, float E, float nu, float Ke[24][24]);

template void ComputeBavg<double>(double a1, double a2, double a3, double Bavg[6][24]);
template void ComputeBavg<float>(float a1, float a2, float a3, float Bavg[6][24]);

template void ComputeS<double>(double a1, double a2, double a3, double E, double nu, double S[6][24]);
template void ComputeS<float>(float a1, float a2, float a3, float E, float nu, float S[6][24]);

#endif // __TUM3D__KE_H__
