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

#endif // __TUM3D__KE_H__
