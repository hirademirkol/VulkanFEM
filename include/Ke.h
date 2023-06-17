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
void ComputeKe(scalar a1, scalar a2, scalar a3, scalar E, scalar nu, scalar Ke[24][24])
{
	memset(Ke, 0, sizeof(scalar)*24*24);

	scalar lambda = E*nu/((1+nu)*(1-2*nu));
	scalar mu = E/(2*(1+nu));

	Ke[0][0] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.2e1 / 0.3e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) - 0.2e1 / 0.3e1 * a1 / a2 / a3 * (scalar) mu) * a3 * a3 / 0.2e1 + (a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) (lambda + 2 * mu)) * pow(a2, 0.3e1) * a3 / 0.3e1 + a1 * a2 / a3 * (scalar) mu + (-0.2e1 / a1 * pow(a2, -0.2e1) * (scalar) mu - 0.2e1 / a1 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1 + (-0.2e1 / 0.3e1 * a1 / a2 * pow(a3, -0.2e1) * (scalar) mu - 0.2e1 / a1 / a2 * (scalar) (lambda + 2 * mu)) * a2 * a2 * a3 / 0.2e1 + 0.1e1 / a1 * (scalar) (lambda + 2 * mu) * a2 * a3 + a1 / a2 * (scalar) mu * a3 + (pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu + pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a3 / 0.3e1;

	Ke[0][1] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) lambda + pow(a3, -0.2e1) * (scalar) mu) * pow(a3, 0.3e1) / 0.3e1 + ((0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + ((-0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 - 0.2e1 / a3 * (scalar) mu - 0.2e1 / a3 * (scalar) lambda) * a3 * a3 / 0.2e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * (scalar) mu - pow(a2, -0.2e1) * (scalar) lambda) * a2 * a2 * a3 / 0.2e1 + (scalar) mu * a3 + (scalar) lambda * a3 + (-pow(a1, -0.2e1) / a2 * (scalar) lambda - pow(a1, -0.2e1) / a2 * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1;

	Ke[0][2] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((-0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) mu * a2 - pow(a3, -0.2e1) * (scalar) lambda * a2) * a3 * a3 / 0.2e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 + ((0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - 0.2e1 / a2 / a3 * (scalar) lambda - 0.2e1 / a2 / a3 * (scalar) mu) * a2 * a2 * a3 / 0.2e1 + (-pow(a1, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1 + a2 * (scalar) mu + a2 * (scalar) lambda;

	Ke[0][3] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (0.2e1 / 0.3e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) - a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 + (a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) (lambda + 2 * mu)) * pow(a2, 0.3e1) * a3 / 0.3e1 - 0.1e1 / a1 * (scalar) (lambda + 2 * mu) * a2 * a3 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu - pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a3 / 0.3e1 + (-a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + 0.2e1 / a1 / a2 * (scalar) (lambda + 2 * mu)) * a2 * a2 * a3 / 0.2e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu + 0.1e1 / a1 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1;

	Ke[0][4] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) mu) * pow(a3, 0.3e1) / 0.3e1 + (((0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 + (-0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + 0.2e1 / a3 * (scalar) mu) * a3 * a3 / 0.2e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * (scalar) mu) * a2 * a2 * a3 / 0.2e1 - (scalar) mu * a3 + (pow(a1, -0.2e1) / a2 * (scalar) lambda + pow(a1, -0.2e1) / a2 * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1;

	Ke[0][5] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) mu * a2) * a3 * a3 / 0.2e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 + ((-0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + 0.2e1 / a2 / a3 * (scalar) mu) * a2 * a2 * a3 / 0.2e1 + (pow(a1, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1 - a2 * (scalar) mu;

	Ke[0][6] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1 + 0.2e1 / 0.3e1 * a1 / a2 / a3 * (scalar) mu) * a3 * a3 / 0.2e1 + (-a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) (lambda + 2 * mu)) * pow(a2, 0.3e1) * a3 / 0.3e1 - a1 / a2 * (scalar) mu * a3 / 0.3e1 + (a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + 0.1e1 / a1 / a2 * (scalar) (lambda + 2 * mu)) * a2 * a2 * a3 / 0.2e1;

	Ke[0][7] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 + 0.1e1 / a3 * (scalar) lambda) * a3 * a3 / 0.2e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * (scalar) mu + pow(a2, -0.2e1) * (scalar) lambda) * a2 * a2 * a3 / 0.2e1 - (scalar) lambda * a3 / 0.2e1;

	Ke[0][8] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 + ((-pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda - pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + 0.1e1 / a2 / a3 * (scalar) lambda + 0.1e1 / a2 / a3 * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[0][9] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 + (-a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) (lambda + 2 * mu)) * pow(a2, 0.3e1) * a3 / 0.3e1 - a1 / a2 * (scalar) mu * a3 / 0.6e1 + (a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - 0.1e1 / a1 / a2 * (scalar) (lambda + 2 * mu)) * a2 * a2 * a3 / 0.2e1;

	Ke[0][10] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((-0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 + 0.1e1 / a3 * (scalar) lambda) * a3 * a3 / 0.2e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * (scalar) mu) * a2 * a2 * a3 / 0.2e1 - (scalar) lambda * a3 / 0.2e1;

	Ke[0][11] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 + ((pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda + pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 / a3 * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[0][12] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.9e1;

	Ke[0][13] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) lambda - pow(a3, -0.2e1) * (scalar) mu) * pow(a3, 0.3e1) / 0.3e1 + (((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda - pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + 0.1e1 / a3 * (scalar) lambda + 0.1e1 / a3 * (scalar) mu) * a3 * a3 / 0.2e1;

	Ke[0][14] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + (-pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) mu * a2 + ((0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) lambda * a2) * a3 * a3 / 0.2e1 - a2 * (scalar) lambda / 0.6e1;

	Ke[0][15] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[0][16] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) mu) * pow(a3, 0.3e1) / 0.3e1 + (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda + pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - 0.1e1 / a3 * (scalar) mu) * a3 * a3 / 0.2e1;

	Ke[0][17] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((-0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) mu * a2) * a3 * a3 / 0.2e1 - a2 * (scalar) lambda / 0.6e1;

	Ke[0][18] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[0][19] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 - 0.1e1 / a3 * (scalar) lambda / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[0][20] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 - a2 * (scalar) lambda / 0.12e2;

	Ke[0][21] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.36e2;

	Ke[0][22] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 - 0.1e1 / a3 * (scalar) lambda / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[0][23] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 - a2 * (scalar) lambda / 0.12e2;

	Ke[1][0] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) lambda + pow(a3, -0.2e1) * (scalar) mu) * pow(a3, 0.3e1) / 0.3e1 + ((0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + ((-0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 - 0.2e1 / a3 * (scalar) mu - 0.2e1 / a3 * (scalar) lambda) * a3 * a3 / 0.2e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * (scalar) mu - pow(a2, -0.2e1) * (scalar) lambda) * a2 * a2 * a3 / 0.2e1 + (scalar) mu * a3 + (scalar) lambda * a3 + (-pow(a1, -0.2e1) / a2 * (scalar) lambda - pow(a1, -0.2e1) / a2 * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1;

	Ke[1][1] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.2e1 / 0.3e1 / a1 * a2 / a3 * (scalar) mu - 0.2e1 / 0.3e1 * a1 / a2 / a3 * (scalar) (lambda + 2 * mu)) * a3 * a3 / 0.2e1 + (a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 + a1 / a2 * (scalar) (lambda + 2 * mu) * a3 + a1 * a2 / a3 * (scalar) mu + (-0.2e1 / 0.3e1 * a1 / a2 * pow(a3, -0.2e1) * (scalar) mu - 0.2e1 / a1 / a2 * (scalar) mu) * a2 * a2 * a3 / 0.2e1 + (-0.2e1 / a1 * pow(a2, -0.2e1) * (scalar) (lambda + 2 * mu) - 0.2e1 / a1 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1 + (pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) (lambda + 2 * mu) + pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a3 / 0.3e1 + 0.1e1 / a1 * (scalar) mu * a2 * a3;

	Ke[1][2] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1 + (0.2e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - a1 * pow(a3, -0.2e1) * (scalar) lambda - a1 * pow(a3, -0.2e1) * (scalar) mu) * a3 * a3 / 0.2e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.2e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.2e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) lambda * a1 - a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 * a3 / 0.2e1 + (pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda + pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * pow(a1, 0.3e1) * a2 * a3 / 0.3e1 + (-0.2e1 / a1 / a2 / a3 * (scalar) lambda - 0.2e1 / a1 / a2 / a3 * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1 + a1 * (scalar) lambda + a1 * (scalar) mu;

	Ke[1][3] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) lambda) * pow(a3, 0.3e1) / 0.3e1 + (((0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) lambda) * a2 * a2 / 0.2e1 + (-0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + 0.2e1 / a3 * (scalar) lambda) * a3 * a3 / 0.2e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * (scalar) lambda) * a2 * a2 * a3 / 0.2e1 - (scalar) lambda * a3 + (pow(a1, -0.2e1) / a2 * (scalar) lambda + pow(a1, -0.2e1) / a2 * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1;

	Ke[1][4] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (0.2e1 / 0.3e1 / a1 * a2 / a3 * (scalar) mu - a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1) * a3 * a3 / 0.2e1 + (a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 - 0.1e1 / a1 * (scalar) mu * a2 * a3 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) (lambda + 2 * mu) - pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a3 / 0.3e1 + (-a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + 0.2e1 / a1 / a2 * (scalar) mu) * a2 * a2 * a3 / 0.2e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) (lambda + 2 * mu) + 0.1e1 / a1 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1;

	Ke[1][5] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1 + (-0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 / 0.2e1 + (-pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda - pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * pow(a1, 0.3e1) * a2 * a3 / 0.3e1 + (0.1e1 / a1 / a2 / a3 * (scalar) lambda + 0.1e1 / a1 / a2 / a3 * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1;

	Ke[1][6] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 + 0.1e1 / a3 * (scalar) mu) * a3 * a3 / 0.2e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * (scalar) mu + pow(a2, -0.2e1) * (scalar) lambda) * a2 * a2 * a3 / 0.2e1 - (scalar) mu * a3 / 0.2e1;

	Ke[1][7] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + 0.2e1 / 0.3e1 * a1 / a2 / a3 * (scalar) (lambda + 2 * mu)) * a3 * a3 / 0.2e1 + (-a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 - a1 / a2 * (scalar) (lambda + 2 * mu) * a3 / 0.3e1 + (a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + 0.1e1 / a1 / a2 * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[1][8] = (((0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 - a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + a1 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 + ((-0.2e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.2e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) lambda * a1 + (pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 * a3 / 0.2e1 - a1 * (scalar) mu / 0.3e1;

	Ke[1][9] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((-0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) lambda) * a2 * a2 / 0.2e1 + 0.1e1 / a3 * (scalar) mu) * a3 * a3 / 0.2e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * (scalar) lambda) * a2 * a2 * a3 / 0.2e1 - (scalar) mu * a3 / 0.2e1;

	Ke[1][10] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1) * a3 * a3 / 0.2e1 + (-a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 - a1 / a2 * (scalar) (lambda + 2 * mu) * a3 / 0.6e1 + (a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - 0.1e1 / a1 / a2 * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[1][11] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + a1 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 / 0.2e1 - a1 * (scalar) mu / 0.6e1;

	Ke[1][12] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) lambda - pow(a3, -0.2e1) * (scalar) mu) * pow(a3, 0.3e1) / 0.3e1 + (((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda - pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + 0.1e1 / a3 * (scalar) lambda + 0.1e1 / a3 * (scalar) mu) * a3 * a3 / 0.2e1;

	Ke[1][13] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.9e1;

	Ke[1][14] = (((0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 - a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1 + (-0.2e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + a1 * pow(a3, -0.2e1) * (scalar) lambda + a1 * pow(a3, -0.2e1) * (scalar) mu) * a3 * a3 / 0.2e1 - a1 * (scalar) lambda / 0.6e1;

	Ke[1][15] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) lambda) * pow(a3, 0.3e1) / 0.3e1 + (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) lambda) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda + pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - 0.1e1 / a3 * (scalar) lambda) * a3 * a3 / 0.2e1;

	Ke[1][16] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[1][17] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1 + (0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1 - a1 * (scalar) lambda / 0.12e2;

	Ke[1][18] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 - 0.1e1 / a3 * (scalar) mu / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[1][19] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[1][20] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 - a1 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * (scalar) lambda / 0.6e1;

	Ke[1][21] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) lambda) * a2 * a2 / 0.2e1 - 0.1e1 / a3 * (scalar) mu / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[1][22] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.36e2;

	Ke[1][23] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 - a1 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * (scalar) lambda / 0.12e2;

	Ke[2][0] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((-0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) mu * a2 - pow(a3, -0.2e1) * (scalar) lambda * a2) * a3 * a3 / 0.2e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 + ((0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - 0.2e1 / a2 / a3 * (scalar) lambda - 0.2e1 / a2 / a3 * (scalar) mu) * a2 * a2 * a3 / 0.2e1 + (-pow(a1, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1 + a2 * (scalar) mu + a2 * (scalar) lambda;

	Ke[2][1] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1 + (0.2e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - a1 * pow(a3, -0.2e1) * (scalar) lambda - a1 * pow(a3, -0.2e1) * (scalar) mu) * a3 * a3 / 0.2e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.2e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.2e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) lambda * a1 - a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 * a3 / 0.2e1 + (pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda + pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * pow(a1, 0.3e1) * a2 * a3 / 0.3e1 + (-0.2e1 / a1 / a2 / a3 * (scalar) lambda - 0.2e1 / a1 / a2 / a3 * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1 + a1 * (scalar) lambda + a1 * (scalar) mu;

	Ke[2][2] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.2e1 / 0.3e1 / a1 * a2 / a3 * (scalar) mu - 0.2e1 / 0.3e1 * a1 / a2 / a3 * (scalar) mu) * a3 * a3 / 0.2e1 + (a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 + 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 + 0.1e1 / a1 * (scalar) mu * a2 * a3 + (pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) + pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a3 / 0.3e1 + (-0.2e1 / 0.3e1 * a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) - 0.2e1 / a1 / a2 * (scalar) mu) * a2 * a2 * a3 / 0.2e1 + a1 * a2 / a3 * (scalar) (lambda + 2 * mu) + (-0.2e1 / a1 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) - 0.2e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1 + a1 / a2 * (scalar) mu * a3;

	Ke[2][3] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) / 0.3e1 + ((0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) lambda * a2) * a3 * a3 / 0.2e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) lambda) * pow(a2, 0.3e1) * a3 / 0.3e1 + ((-0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + 0.2e1 / a2 / a3 * (scalar) lambda) * a2 * a2 * a3 / 0.2e1 + (pow(a1, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1 - a2 * (scalar) lambda;

	Ke[2][4] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1 + (-0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 / 0.2e1 + (-pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda - pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * pow(a1, 0.3e1) * a2 * a3 / 0.3e1 + (0.1e1 / a1 / a2 / a3 * (scalar) lambda + 0.1e1 / a1 / a2 / a3 * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1;

	Ke[2][5] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (0.2e1 / 0.3e1 / a1 * a2 / a3 * (scalar) mu - a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 + (a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 - 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 - 0.1e1 / a1 * (scalar) mu * a2 * a3 + (-pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) - pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a3 / 0.3e1 + (-a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 + 0.2e1 / a1 / a2 * (scalar) mu) * a2 * a2 * a3 / 0.2e1 + (0.1e1 / a1 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) + 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1;

	Ke[2][6] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 + ((-pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda - pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + 0.1e1 / a2 / a3 * (scalar) lambda + 0.1e1 / a2 / a3 * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[2][7] = (((0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 - a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.3e1) * a3 * a3 / 0.2e1 + ((-0.2e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.2e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) lambda * a1 + (pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 * a3 / 0.2e1 - a1 * (scalar) lambda / 0.3e1;

	Ke[2][8] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + 0.2e1 / 0.3e1 * a1 / a2 / a3 * (scalar) mu) * a3 * a3 / 0.2e1 + (-a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 - 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 - a1 / a2 * (scalar) mu * a3 / 0.3e1 + (a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 + 0.1e1 / a1 / a2 * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[2][9] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) lambda) * pow(a2, 0.3e1) * a3 / 0.3e1 + ((pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda + pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 / a3 * (scalar) lambda) * a2 * a2 * a3 / 0.2e1;

	Ke[2][10] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.6e1) * a3 * a3 / 0.2e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 / 0.2e1 - a1 * (scalar) lambda / 0.6e1;

	Ke[2][11] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 + (-a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 + 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 - a1 / a2 * (scalar) mu * a3 / 0.6e1 + (a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 - 0.1e1 / a1 / a2 * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[2][12] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + (-pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) mu * a2 + ((0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) lambda * a2) * a3 * a3 / 0.2e1 - a2 * (scalar) mu / 0.6e1;

	Ke[2][13] = (((0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 - a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1 + (-0.2e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + a1 * pow(a3, -0.2e1) * (scalar) lambda + a1 * pow(a3, -0.2e1) * (scalar) mu) * a3 * a3 / 0.2e1 - a1 * (scalar) mu / 0.6e1;

	Ke[2][14] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.9e1;

	Ke[2][15] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) / 0.3e1 + ((-0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) lambda * a2) * a3 * a3 / 0.2e1 - a2 * (scalar) mu / 0.6e1;

	Ke[2][16] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1 + (0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1 - a1 * (scalar) mu / 0.12e2;

	Ke[2][17] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.18e2;

	Ke[2][18] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 - a2 * (scalar) mu / 0.12e2;

	Ke[2][19] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 - a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.3e1) * a3 * a3 / 0.2e1 - a1 * (scalar) mu / 0.6e1;

	Ke[2][20] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.18e2;

	Ke[2][21] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) / 0.3e1 + ((pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 - a2 * (scalar) mu / 0.12e2;

	Ke[2][22] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 - a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.6e1) * a3 * a3 / 0.2e1 - a1 * (scalar) mu / 0.12e2;

	Ke[2][23] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.36e2;

	Ke[3][0] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (0.2e1 / 0.3e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) - a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 + (a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) (lambda + 2 * mu)) * pow(a2, 0.3e1) * a3 / 0.3e1 - 0.1e1 / a1 * (scalar) (lambda + 2 * mu) * a2 * a3 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu - pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a3 / 0.3e1 + (-a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + 0.2e1 / a1 / a2 * (scalar) (lambda + 2 * mu)) * a2 * a2 * a3 / 0.2e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu + 0.1e1 / a1 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1;

	Ke[3][1] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) lambda) * pow(a3, 0.3e1) / 0.3e1 + (((0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) lambda) * a2 * a2 / 0.2e1 + (-0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + 0.2e1 / a3 * (scalar) lambda) * a3 * a3 / 0.2e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * (scalar) lambda) * a2 * a2 * a3 / 0.2e1 - (scalar) lambda * a3 + (pow(a1, -0.2e1) / a2 * (scalar) lambda + pow(a1, -0.2e1) / a2 * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1;

	Ke[3][2] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) / 0.3e1 + ((0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) lambda * a2) * a3 * a3 / 0.2e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) lambda) * pow(a2, 0.3e1) * a3 / 0.3e1 + ((-0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + 0.2e1 / a2 / a3 * (scalar) lambda) * a2 * a2 * a3 / 0.2e1 + (pow(a1, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1 - a2 * (scalar) lambda;

	Ke[3][3] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.2e1 / 0.3e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) - 0.2e1 / 0.3e1 * a1 / a2 / a3 * (scalar) mu) * a3 * a3 / 0.2e1 + (a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) (lambda + 2 * mu)) * pow(a2, 0.3e1) * a3 / 0.3e1 + 0.1e1 / a1 * (scalar) (lambda + 2 * mu) * a2 * a3 + (pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu + pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a3 / 0.3e1 + (-0.2e1 / 0.3e1 * a1 / a2 * pow(a3, -0.2e1) * (scalar) mu - 0.2e1 / a1 / a2 * (scalar) (lambda + 2 * mu)) * a2 * a2 * a3 / 0.2e1;

	Ke[3][4] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + ((-0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + (0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + (pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 * a3 / 0.4e1 + (-pow(a1, -0.2e1) / a2 * (scalar) lambda - pow(a1, -0.2e1) / a2 * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1;

	Ke[3][5] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) / 0.6e1 + (-0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + (pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) * a3 / 0.6e1 + (0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 * a2 * a3 / 0.4e1 + (-pow(a1, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1;

	Ke[3][6] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 + (-a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) (lambda + 2 * mu)) * pow(a2, 0.3e1) * a3 / 0.3e1 - a1 / a2 * (scalar) mu * a3 / 0.6e1 + (a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - 0.1e1 / a1 / a2 * (scalar) (lambda + 2 * mu)) * a2 * a2 * a3 / 0.2e1;

	Ke[3][7] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((-0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) lambda) * a2 * a2 / 0.2e1 - 0.1e1 / a3 * (scalar) lambda) * a3 * a3 / 0.2e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * (scalar) lambda) * a2 * a2 * a3 / 0.2e1 + (scalar) lambda * a3 / 0.2e1;

	Ke[3][8] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) lambda) * pow(a2, 0.3e1) * a3 / 0.3e1 + ((pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda + pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 / a3 * (scalar) lambda) * a2 * a2 * a3 / 0.2e1;

	Ke[3][9] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1 + 0.2e1 / 0.3e1 * a1 / a2 / a3 * (scalar) mu) * a3 * a3 / 0.2e1 + (-a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) (lambda + 2 * mu)) * pow(a2, 0.3e1) * a3 / 0.3e1 - a1 / a2 * (scalar) mu * a3 / 0.3e1 + (a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + 0.1e1 / a1 / a2 * (scalar) (lambda + 2 * mu)) * a2 * a2 * a3 / 0.2e1;

	Ke[3][10] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + ((0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 - 0.1e1 / a3 * (scalar) lambda) * a3 * a3 / 0.2e1 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 * a3 / 0.4e1 + (scalar) lambda * a3 / 0.2e1;

	Ke[3][11] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) / 0.6e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1) * a3 * a3 / 0.2e1 + (pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) * a3 / 0.6e1 + (-pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda - pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 * a2 * a3 / 0.4e1;

	Ke[3][12] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[3][13] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) lambda) * pow(a3, 0.3e1) / 0.3e1 + (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) lambda) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda + pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - 0.1e1 / a3 * (scalar) lambda) * a3 * a3 / 0.2e1;

	Ke[3][14] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) / 0.3e1 + ((-0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) lambda * a2) * a3 * a3 / 0.2e1 + a2 * (scalar) lambda / 0.6e1;

	Ke[3][15] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.9e1;

	Ke[3][16] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + (-pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda - pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[3][17] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) / 0.6e1 + (0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + (-pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + a2 * (scalar) lambda / 0.6e1;

	Ke[3][18] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.36e2;

	Ke[3][19] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) lambda) * a2 * a2 / 0.2e1 + 0.1e1 / a3 * (scalar) lambda / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[3][20] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) / 0.3e1 + ((pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + a2 * (scalar) lambda / 0.12e2;

	Ke[3][21] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[3][22] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 - pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + 0.1e1 / a3 * (scalar) lambda / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[3][23] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) / 0.6e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1) * a3 * a3 / 0.2e1 + a2 * (scalar) lambda / 0.12e2;

	Ke[4][0] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) mu) * pow(a3, 0.3e1) / 0.3e1 + (((0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 + (-0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + 0.2e1 / a3 * (scalar) mu) * a3 * a3 / 0.2e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * (scalar) mu) * a2 * a2 * a3 / 0.2e1 - (scalar) mu * a3 + (pow(a1, -0.2e1) / a2 * (scalar) lambda + pow(a1, -0.2e1) / a2 * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1;

	Ke[4][1] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (0.2e1 / 0.3e1 / a1 * a2 / a3 * (scalar) mu - a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1) * a3 * a3 / 0.2e1 + (a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 - 0.1e1 / a1 * (scalar) mu * a2 * a3 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) (lambda + 2 * mu) - pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a3 / 0.3e1 + (-a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + 0.2e1 / a1 / a2 * (scalar) mu) * a2 * a2 * a3 / 0.2e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) (lambda + 2 * mu) + 0.1e1 / a1 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1;

	Ke[4][2] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1 + (-0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 / 0.2e1 + (-pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda - pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * pow(a1, 0.3e1) * a2 * a3 / 0.3e1 + (0.1e1 / a1 / a2 / a3 * (scalar) lambda + 0.1e1 / a1 / a2 / a3 * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1;

	Ke[4][3] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + ((-0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + (0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + (pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 * a3 / 0.4e1 + (-pow(a1, -0.2e1) / a2 * (scalar) lambda - pow(a1, -0.2e1) / a2 * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1;

	Ke[4][4] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.2e1 / 0.3e1 / a1 * a2 / a3 * (scalar) mu - 0.2e1 / 0.3e1 * a1 / a2 / a3 * (scalar) (lambda + 2 * mu)) * a3 * a3 / 0.2e1 + (a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 + 0.1e1 / a1 * (scalar) mu * a2 * a3 + (pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) (lambda + 2 * mu) + pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a3 / 0.3e1 + (-0.2e1 / 0.3e1 * a1 / a2 * pow(a3, -0.2e1) * (scalar) mu - 0.2e1 / a1 / a2 * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[4][5] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 / 0.6e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1) * a3 * a3 / 0.2e1 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 * a3 / 0.6e1 + (pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda + pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * pow(a1, 0.3e1) * a2 * a3 / 0.3e1;

	Ke[4][6] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((-0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 - 0.1e1 / a3 * (scalar) mu) * a3 * a3 / 0.2e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * (scalar) mu) * a2 * a2 * a3 / 0.2e1 + (scalar) mu * a3 / 0.2e1;

	Ke[4][7] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1) * a3 * a3 / 0.2e1 + (-a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 - a1 / a2 * (scalar) (lambda + 2 * mu) * a3 / 0.6e1 + (a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - 0.1e1 / a1 / a2 * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[4][8] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + a1 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 / 0.2e1 - a1 * (scalar) mu / 0.6e1;

	Ke[4][9] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + ((0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 - 0.1e1 / a3 * (scalar) mu) * a3 * a3 / 0.2e1 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 * a3 / 0.4e1 + (scalar) mu * a3 / 0.2e1;

	Ke[4][10] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + 0.2e1 / 0.3e1 * a1 / a2 / a3 * (scalar) (lambda + 2 * mu)) * a3 * a3 / 0.2e1 + (-a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 - a1 / a2 * (scalar) (lambda + 2 * mu) * a3 / 0.3e1 + (a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + 0.1e1 / a1 / a2 * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[4][11] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 / 0.6e1 + a1 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 + (pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 * a3 / 0.6e1 - a1 * (scalar) mu / 0.3e1;

	Ke[4][12] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) mu) * pow(a3, 0.3e1) / 0.3e1 + (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda + pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - 0.1e1 / a3 * (scalar) mu) * a3 * a3 / 0.2e1;

	Ke[4][13] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[4][14] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1 + (0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1 - a1 * (scalar) lambda / 0.12e2;

	Ke[4][15] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + (-pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda - pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[4][16] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.9e1;

	Ke[4][17] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 / 0.6e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1) * a3 * a3 / 0.2e1 - a1 * (scalar) lambda / 0.6e1;

	Ke[4][18] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 + 0.1e1 / a3 * (scalar) mu / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[4][19] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.36e2;

	Ke[4][20] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 - a1 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * (scalar) lambda / 0.12e2;

	Ke[4][21] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 - pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + 0.1e1 / a3 * (scalar) mu / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[4][22] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[4][23] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 / 0.6e1 - a1 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * (scalar) lambda / 0.6e1;

	Ke[5][0] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) mu * a2) * a3 * a3 / 0.2e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 + ((-0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + 0.2e1 / a2 / a3 * (scalar) mu) * a2 * a2 * a3 / 0.2e1 + (pow(a1, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1 - a2 * (scalar) mu;

	Ke[5][1] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1 + (-0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 / 0.2e1 + (-pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda - pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * pow(a1, 0.3e1) * a2 * a3 / 0.3e1 + (0.1e1 / a1 / a2 / a3 * (scalar) lambda + 0.1e1 / a1 / a2 / a3 * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1;

	Ke[5][2] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (0.2e1 / 0.3e1 / a1 * a2 / a3 * (scalar) mu - a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 + (a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 - 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 - 0.1e1 / a1 * (scalar) mu * a2 * a3 + (-pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) - pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a3 / 0.3e1 + (-a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 + 0.2e1 / a1 / a2 * (scalar) mu) * a2 * a2 * a3 / 0.2e1 + (0.1e1 / a1 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) + 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1;

	Ke[5][3] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) / 0.6e1 + (-0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + (pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) * a3 / 0.6e1 + (0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 * a2 * a3 / 0.4e1 + (-pow(a1, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a3 / 0.2e1;

	Ke[5][4] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 / 0.6e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1) * a3 * a3 / 0.2e1 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 * a3 / 0.6e1 + (pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda + pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * pow(a1, 0.3e1) * a2 * a3 / 0.3e1;

	Ke[5][5] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.2e1 / 0.3e1 / a1 * a2 / a3 * (scalar) mu - 0.2e1 / 0.3e1 * a1 / a2 / a3 * (scalar) mu) * a3 * a3 / 0.2e1 + (a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 + 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 + 0.1e1 / a1 * (scalar) mu * a2 * a3 + (pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) + pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a3 / 0.3e1 + (-0.2e1 / 0.3e1 * a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) - 0.2e1 / a1 / a2 * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[5][6] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 + ((pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda + pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 / a3 * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[5][7] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.6e1) * a3 * a3 / 0.2e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 / 0.2e1 - a1 * (scalar) lambda / 0.6e1;

	Ke[5][8] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 + (-a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 + 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 - a1 / a2 * (scalar) mu * a3 / 0.6e1 + (a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 - 0.1e1 / a1 / a2 * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[5][9] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) / 0.6e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1) * a3 * a3 / 0.2e1 + (pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) * a3 / 0.6e1 + (-pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda - pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 * a2 * a3 / 0.4e1;

	Ke[5][10] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 / 0.6e1 + a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.3e1) * a3 * a3 / 0.2e1 + (pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 * a3 / 0.6e1 - a1 * (scalar) lambda / 0.3e1;

	Ke[5][11] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + 0.2e1 / 0.3e1 * a1 / a2 / a3 * (scalar) mu) * a3 * a3 / 0.2e1 + (-a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 - 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 - a1 / a2 * (scalar) mu * a3 / 0.3e1 + (a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 + 0.1e1 / a1 / a2 * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[5][12] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((-0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) mu * a2) * a3 * a3 / 0.2e1 + a2 * (scalar) mu / 0.6e1;

	Ke[5][13] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1 + (0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1 - a1 * (scalar) mu / 0.12e2;

	Ke[5][14] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.18e2;

	Ke[5][15] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) / 0.6e1 + (0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + (-pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + a2 * (scalar) mu / 0.6e1;

	Ke[5][16] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 / 0.6e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1) * a3 * a3 / 0.2e1 - a1 * (scalar) mu / 0.6e1;

	Ke[5][17] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.9e1;

	Ke[5][18] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + a2 * (scalar) mu / 0.12e2;

	Ke[5][19] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 - a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.6e1) * a3 * a3 / 0.2e1 - a1 * (scalar) mu / 0.12e2;

	Ke[5][20] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.36e2;

	Ke[5][21] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) / 0.6e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1) * a3 * a3 / 0.2e1 + a2 * (scalar) mu / 0.12e2;

	Ke[5][22] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 / 0.6e1 - a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.3e1) * a3 * a3 / 0.2e1 - a1 * (scalar) mu / 0.6e1;

	Ke[5][23] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.18e2;

	Ke[6][0] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1 + 0.2e1 / 0.3e1 * a1 / a2 / a3 * (scalar) mu) * a3 * a3 / 0.2e1 + (-a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) (lambda + 2 * mu)) * pow(a2, 0.3e1) * a3 / 0.3e1 - a1 / a2 * (scalar) mu * a3 / 0.3e1 + (a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + 0.1e1 / a1 / a2 * (scalar) (lambda + 2 * mu)) * a2 * a2 * a3 / 0.2e1;

	Ke[6][1] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 + 0.1e1 / a3 * (scalar) mu) * a3 * a3 / 0.2e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * (scalar) mu + pow(a2, -0.2e1) * (scalar) lambda) * a2 * a2 * a3 / 0.2e1 - (scalar) mu * a3 / 0.2e1;

	Ke[6][2] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 + ((-pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda - pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + 0.1e1 / a2 / a3 * (scalar) lambda + 0.1e1 / a2 / a3 * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[6][3] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 + (-a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) (lambda + 2 * mu)) * pow(a2, 0.3e1) * a3 / 0.3e1 - a1 / a2 * (scalar) mu * a3 / 0.6e1 + (a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - 0.1e1 / a1 / a2 * (scalar) (lambda + 2 * mu)) * a2 * a2 * a3 / 0.2e1;

	Ke[6][4] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((-0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 - 0.1e1 / a3 * (scalar) mu) * a3 * a3 / 0.2e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * (scalar) mu) * a2 * a2 * a3 / 0.2e1 + (scalar) mu * a3 / 0.2e1;

	Ke[6][5] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 + ((pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda + pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 / a3 * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[6][6] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.2e1 / 0.3e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) - 0.2e1 / 0.3e1 * a1 / a2 / a3 * (scalar) mu) * a3 * a3 / 0.2e1 + (a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) (lambda + 2 * mu)) * pow(a2, 0.3e1) * a3 / 0.3e1 + a1 / a2 * (scalar) mu * a3 / 0.3e1;

	Ke[6][7] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 * pow(a3, 0.3e1) / 0.6e1 + ((-0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 * a3 * a3 / 0.4e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * (scalar) mu - pow(a2, -0.2e1) * (scalar) lambda) * a2 * a2 * a3 / 0.2e1;

	Ke[6][8] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 * a3 / 0.6e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1;

	Ke[6][9] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (0.2e1 / 0.3e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) - a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 + (a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) (lambda + 2 * mu)) * pow(a2, 0.3e1) * a3 / 0.3e1 + a1 / a2 * (scalar) mu * a3 / 0.6e1;

	Ke[6][10] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 * pow(a3, 0.3e1) / 0.6e1 + ((0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 * a3 * a3 / 0.4e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[6][11] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 * a3 / 0.6e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1;

	Ke[6][12] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[6][13] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 - 0.1e1 / a3 * (scalar) mu / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[6][14] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 - a2 * (scalar) lambda / 0.12e2;

	Ke[6][15] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.36e2;

	Ke[6][16] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 + 0.1e1 / a3 * (scalar) mu / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[6][17] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 - a2 * (scalar) lambda / 0.12e2;

	Ke[6][18] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.9e1;

	Ke[6][19] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 * pow(a3, 0.3e1) / 0.6e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 * a3 * a3 / 0.4e1;

	Ke[6][20] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 * a3 / 0.6e1 - a2 * (scalar) lambda / 0.6e1;

	Ke[6][21] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[6][22] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 * pow(a3, 0.3e1) / 0.6e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 * a3 * a3 / 0.4e1;

	Ke[6][23] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 * a3 / 0.6e1 - a2 * (scalar) lambda / 0.6e1;

	Ke[7][0] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 + 0.1e1 / a3 * (scalar) lambda) * a3 * a3 / 0.2e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * (scalar) mu + pow(a2, -0.2e1) * (scalar) lambda) * a2 * a2 * a3 / 0.2e1 - (scalar) lambda * a3 / 0.2e1;

	Ke[7][1] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + 0.2e1 / 0.3e1 * a1 / a2 / a3 * (scalar) (lambda + 2 * mu)) * a3 * a3 / 0.2e1 + (-a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 - a1 / a2 * (scalar) (lambda + 2 * mu) * a3 / 0.3e1 + (a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + 0.1e1 / a1 / a2 * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[7][2] = (((0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 - a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.3e1) * a3 * a3 / 0.2e1 + ((-0.2e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.2e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) lambda * a1 + (pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 * a3 / 0.2e1 - a1 * (scalar) lambda / 0.3e1;

	Ke[7][3] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((-0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) lambda) * a2 * a2 / 0.2e1 - 0.1e1 / a3 * (scalar) lambda) * a3 * a3 / 0.2e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * (scalar) lambda) * a2 * a2 * a3 / 0.2e1 + (scalar) lambda * a3 / 0.2e1;

	Ke[7][4] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1) * a3 * a3 / 0.2e1 + (-a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 - a1 / a2 * (scalar) (lambda + 2 * mu) * a3 / 0.6e1 + (a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - 0.1e1 / a1 / a2 * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[7][5] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.6e1) * a3 * a3 / 0.2e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 / 0.2e1 - a1 * (scalar) lambda / 0.6e1;

	Ke[7][6] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 * pow(a3, 0.3e1) / 0.6e1 + ((-0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 * a3 * a3 / 0.4e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * (scalar) mu - pow(a2, -0.2e1) * (scalar) lambda) * a2 * a2 * a3 / 0.2e1;

	Ke[7][7] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.2e1 / 0.3e1 / a1 * a2 / a3 * (scalar) mu - 0.2e1 / 0.3e1 * a1 / a2 / a3 * (scalar) (lambda + 2 * mu)) * a3 * a3 / 0.2e1 + (a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 + a1 / a2 * (scalar) (lambda + 2 * mu) * a3 / 0.3e1;

	Ke[7][8] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 * a3 * a3 / 0.4e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.2e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.2e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) lambda * a1 - a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[7][9] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 * pow(a3, 0.3e1) / 0.6e1 + ((0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) lambda) * a2 * a2 * a3 * a3 / 0.4e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * (scalar) lambda) * a2 * a2 * a3 / 0.2e1;

	Ke[7][10] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (0.2e1 / 0.3e1 / a1 * a2 / a3 * (scalar) mu - a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1) * a3 * a3 / 0.2e1 + (a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 + a1 / a2 * (scalar) (lambda + 2 * mu) * a3 / 0.6e1;

	Ke[7][11] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 * a3 / 0.4e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 / 0.2e1;

	Ke[7][12] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 - 0.1e1 / a3 * (scalar) lambda / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[7][13] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[7][14] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 - a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.3e1) * a3 * a3 / 0.2e1 + a1 * (scalar) lambda / 0.6e1;

	Ke[7][15] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) lambda) * a2 * a2 / 0.2e1 + 0.1e1 / a3 * (scalar) lambda / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[7][16] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.36e2;

	Ke[7][17] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 - a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.6e1) * a3 * a3 / 0.2e1 + a1 * (scalar) lambda / 0.12e2;

	Ke[7][18] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 * pow(a3, 0.3e1) / 0.6e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 * a3 * a3 / 0.4e1;

	Ke[7][19] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.9e1;

	Ke[7][20] = ((0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 - a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 * a3 * a3 / 0.4e1 + a1 * (scalar) lambda / 0.6e1;

	Ke[7][21] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 * pow(a3, 0.3e1) / 0.6e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) lambda) * a2 * a2 * a3 * a3 / 0.4e1;

	Ke[7][22] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[7][23] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 * a3 / 0.4e1 + a1 * (scalar) lambda / 0.12e2;

	Ke[8][0] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 + ((-pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda - pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + 0.1e1 / a2 / a3 * (scalar) lambda + 0.1e1 / a2 / a3 * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[8][1] = (((0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 - a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + a1 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 + ((-0.2e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.2e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) lambda * a1 + (pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 * a3 / 0.2e1 - a1 * (scalar) mu / 0.3e1;

	Ke[8][2] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + 0.2e1 / 0.3e1 * a1 / a2 / a3 * (scalar) mu) * a3 * a3 / 0.2e1 + (-a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 - 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 - a1 / a2 * (scalar) mu * a3 / 0.3e1 + (a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 + 0.1e1 / a1 / a2 * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[8][3] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) lambda) * pow(a2, 0.3e1) * a3 / 0.3e1 + ((pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda + pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 / a3 * (scalar) lambda) * a2 * a2 * a3 / 0.2e1;

	Ke[8][4] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + a1 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 / 0.2e1 - a1 * (scalar) mu / 0.6e1;

	Ke[8][5] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 + (-a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 + 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 - a1 / a2 * (scalar) mu * a3 / 0.6e1 + (a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 - 0.1e1 / a1 / a2 * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[8][6] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 * a3 / 0.6e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1;

	Ke[8][7] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 * a3 * a3 / 0.4e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.2e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.2e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) lambda * a1 - a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[8][8] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.2e1 / 0.3e1 / a1 * a2 / a3 * (scalar) mu - 0.2e1 / 0.3e1 * a1 / a2 / a3 * (scalar) mu) * a3 * a3 / 0.2e1 + (a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 + 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 + a1 / a2 * (scalar) mu * a3 / 0.3e1;

	Ke[8][9] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) * a3 * a3 / 0.6e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) lambda) * pow(a2, 0.3e1) * a3 / 0.3e1;

	Ke[8][10] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 * a3 / 0.4e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 / 0.2e1;

	Ke[8][11] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (0.2e1 / 0.3e1 / a1 * a2 / a3 * (scalar) mu - a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 + (a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 - 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 + a1 / a2 * (scalar) mu * a3 / 0.6e1;

	Ke[8][12] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 - a2 * (scalar) mu / 0.12e2;

	Ke[8][13] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 - a1 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 + a1 * (scalar) mu / 0.6e1;

	Ke[8][14] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.18e2;

	Ke[8][15] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) / 0.3e1 + ((pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 - a2 * (scalar) mu / 0.12e2;

	Ke[8][16] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 - a1 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 + a1 * (scalar) mu / 0.12e2;

	Ke[8][17] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.36e2;

	Ke[8][18] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 * a3 / 0.6e1 - a2 * (scalar) mu / 0.6e1;

	Ke[8][19] = ((0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 - a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 * a3 * a3 / 0.4e1 + a1 * (scalar) mu / 0.6e1;

	Ke[8][20] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.9e1;

	Ke[8][21] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) * a3 * a3 / 0.6e1 - a2 * (scalar) mu / 0.6e1;

	Ke[8][22] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 * a3 / 0.4e1 + a1 * (scalar) mu / 0.12e2;

	Ke[8][23] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.18e2;

	Ke[9][0] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 + (-a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) (lambda + 2 * mu)) * pow(a2, 0.3e1) * a3 / 0.3e1 - a1 / a2 * (scalar) mu * a3 / 0.6e1 + (a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - 0.1e1 / a1 / a2 * (scalar) (lambda + 2 * mu)) * a2 * a2 * a3 / 0.2e1;

	Ke[9][1] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((-0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) lambda) * a2 * a2 / 0.2e1 + 0.1e1 / a3 * (scalar) mu) * a3 * a3 / 0.2e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * (scalar) lambda) * a2 * a2 * a3 / 0.2e1 - (scalar) mu * a3 / 0.2e1;

	Ke[9][2] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) lambda) * pow(a2, 0.3e1) * a3 / 0.3e1 + ((pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda + pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 / a3 * (scalar) lambda) * a2 * a2 * a3 / 0.2e1;

	Ke[9][3] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1 + 0.2e1 / 0.3e1 * a1 / a2 / a3 * (scalar) mu) * a3 * a3 / 0.2e1 + (-a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) (lambda + 2 * mu)) * pow(a2, 0.3e1) * a3 / 0.3e1 - a1 / a2 * (scalar) mu * a3 / 0.3e1 + (a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + 0.1e1 / a1 / a2 * (scalar) (lambda + 2 * mu)) * a2 * a2 * a3 / 0.2e1;

	Ke[9][4] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + ((0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 - 0.1e1 / a3 * (scalar) mu) * a3 * a3 / 0.2e1 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 * a3 / 0.4e1 + (scalar) mu * a3 / 0.2e1;

	Ke[9][5] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) / 0.6e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1) * a3 * a3 / 0.2e1 + (pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) * a3 / 0.6e1 + (-pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda - pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 * a2 * a3 / 0.4e1;

	Ke[9][6] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (0.2e1 / 0.3e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) - a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 + (a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) (lambda + 2 * mu)) * pow(a2, 0.3e1) * a3 / 0.3e1 + a1 / a2 * (scalar) mu * a3 / 0.6e1;

	Ke[9][7] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 * pow(a3, 0.3e1) / 0.6e1 + ((0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) lambda) * a2 * a2 * a3 * a3 / 0.4e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * (scalar) lambda) * a2 * a2 * a3 / 0.2e1;

	Ke[9][8] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) * a3 * a3 / 0.6e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) lambda) * pow(a2, 0.3e1) * a3 / 0.3e1;

	Ke[9][9] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.2e1 / 0.3e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) - 0.2e1 / 0.3e1 * a1 / a2 / a3 * (scalar) mu) * a3 * a3 / 0.2e1 + (a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) (lambda + 2 * mu)) * pow(a2, 0.3e1) * a3 / 0.3e1 + a1 / a2 * (scalar) mu * a3 / 0.3e1;

	Ke[9][10] = (pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 * pow(a3, 0.3e1) / 0.12e2 + (-0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a2 * a3 * a3 / 0.8e1 + (pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 * a3 / 0.4e1;

	Ke[9][11] = (pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) * a3 * a3 / 0.12e2 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) * a3 / 0.6e1;

	Ke[9][12] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.36e2;

	Ke[9][13] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) lambda) * a2 * a2 / 0.2e1 - 0.1e1 / a3 * (scalar) mu / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[9][14] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) / 0.3e1 + ((pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + a2 * (scalar) lambda / 0.12e2;

	Ke[9][15] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[9][16] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 - pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + 0.1e1 / a3 * (scalar) mu / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[9][17] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) / 0.6e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1) * a3 * a3 / 0.2e1 + a2 * (scalar) lambda / 0.12e2;

	Ke[9][18] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[9][19] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 * pow(a3, 0.3e1) / 0.6e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) lambda) * a2 * a2 * a3 * a3 / 0.4e1;

	Ke[9][20] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) * a3 * a3 / 0.6e1 + a2 * (scalar) lambda / 0.6e1;

	Ke[9][21] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.9e1;

	Ke[9][22] = (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 * pow(a3, 0.3e1) / 0.12e2 + (pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a2 * a3 * a3 / 0.8e1;

	Ke[9][23] = (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) * a3 * a3 / 0.12e2 + a2 * (scalar) lambda / 0.6e1;

	Ke[10][0] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((-0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 + 0.1e1 / a3 * (scalar) lambda) * a3 * a3 / 0.2e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * (scalar) mu) * a2 * a2 * a3 / 0.2e1 - (scalar) lambda * a3 / 0.2e1;

	Ke[10][1] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1) * a3 * a3 / 0.2e1 + (-a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 - a1 / a2 * (scalar) (lambda + 2 * mu) * a3 / 0.6e1 + (a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - 0.1e1 / a1 / a2 * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[10][2] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.6e1) * a3 * a3 / 0.2e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 / 0.2e1 - a1 * (scalar) lambda / 0.6e1;

	Ke[10][3] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + ((0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 - 0.1e1 / a3 * (scalar) lambda) * a3 * a3 / 0.2e1 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 * a3 / 0.4e1 + (scalar) lambda * a3 / 0.2e1;

	Ke[10][4] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + 0.2e1 / 0.3e1 * a1 / a2 / a3 * (scalar) (lambda + 2 * mu)) * a3 * a3 / 0.2e1 + (-a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 - a1 / a2 * (scalar) (lambda + 2 * mu) * a3 / 0.3e1 + (a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + 0.1e1 / a1 / a2 * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[10][5] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 / 0.6e1 + a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.3e1) * a3 * a3 / 0.2e1 + (pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 * a3 / 0.6e1 - a1 * (scalar) lambda / 0.3e1;

	Ke[10][6] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 * pow(a3, 0.3e1) / 0.6e1 + ((0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - 0.2e1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 * a3 * a3 / 0.4e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[10][7] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (0.2e1 / 0.3e1 / a1 * a2 / a3 * (scalar) mu - a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1) * a3 * a3 / 0.2e1 + (a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 + a1 / a2 * (scalar) (lambda + 2 * mu) * a3 / 0.6e1;

	Ke[10][8] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 * a3 / 0.4e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 / 0.2e1;

	Ke[10][9] = (pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 * pow(a3, 0.3e1) / 0.12e2 + (-0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a2 * a3 * a3 / 0.8e1 + (pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 * a3 / 0.4e1;

	Ke[10][10] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.2e1 / 0.3e1 / a1 * a2 / a3 * (scalar) mu - 0.2e1 / 0.3e1 * a1 / a2 / a3 * (scalar) (lambda + 2 * mu)) * a3 * a3 / 0.2e1 + (a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 + a1 / a2 * (scalar) (lambda + 2 * mu) * a3 / 0.3e1;

	Ke[10][11] = (pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 * a3 * a3 / 0.12e2 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 * a3 / 0.6e1;

	Ke[10][12] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 - 0.1e1 / a3 * (scalar) lambda / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[10][13] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.36e2;

	Ke[10][14] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 - a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.6e1) * a3 * a3 / 0.2e1 + a1 * (scalar) lambda / 0.12e2;

	Ke[10][15] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 - pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + 0.1e1 / a3 * (scalar) lambda / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[10][16] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[10][17] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 / 0.6e1 - a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.3e1) * a3 * a3 / 0.2e1 + a1 * (scalar) lambda / 0.6e1;

	Ke[10][18] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 * pow(a3, 0.3e1) / 0.6e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 * a3 * a3 / 0.4e1;

	Ke[10][19] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[10][20] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 * a3 / 0.4e1 + a1 * (scalar) lambda / 0.12e2;

	Ke[10][21] = (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 * pow(a3, 0.3e1) / 0.12e2 + (pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a2 * a3 * a3 / 0.8e1;

	Ke[10][22] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.9e1;

	Ke[10][23] = (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 * a3 * a3 / 0.12e2 + a1 * (scalar) lambda / 0.6e1;

	Ke[11][0] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 + ((pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda + pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 / a3 * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[11][1] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + a1 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 / 0.2e1 - a1 * (scalar) mu / 0.6e1;

	Ke[11][2] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 + (-a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 + 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 - a1 / a2 * (scalar) mu * a3 / 0.6e1 + (a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 - 0.1e1 / a1 / a2 * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[11][3] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) / 0.6e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1) * a3 * a3 / 0.2e1 + (pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) * a3 / 0.6e1 + (-pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda - pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 * a2 * a3 / 0.4e1;

	Ke[11][4] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 / 0.6e1 + a1 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 + (pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 * a3 / 0.6e1 - a1 * (scalar) mu / 0.3e1;

	Ke[11][5] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + 0.2e1 / 0.3e1 * a1 / a2 / a3 * (scalar) mu) * a3 * a3 / 0.2e1 + (-a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 - 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 - a1 / a2 * (scalar) mu * a3 / 0.3e1 + (a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 + 0.1e1 / a1 / a2 * (scalar) mu) * a2 * a2 * a3 / 0.2e1;

	Ke[11][6] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 * a3 / 0.6e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1;

	Ke[11][7] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 * a3 / 0.4e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 / 0.2e1;

	Ke[11][8] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (0.2e1 / 0.3e1 / a1 * a2 / a3 * (scalar) mu - a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 + (a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 - 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 + a1 / a2 * (scalar) mu * a3 / 0.6e1;

	Ke[11][9] = (pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) * a3 * a3 / 0.12e2 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) * a3 / 0.6e1;

	Ke[11][10] = (pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 * a3 * a3 / 0.12e2 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 * a3 / 0.6e1;

	Ke[11][11] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.2e1 / 0.3e1 / a1 * a2 / a3 * (scalar) mu - 0.2e1 / 0.3e1 * a1 / a2 / a3 * (scalar) mu) * a3 * a3 / 0.2e1 + (a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 + 0.1e1 / a1 * pow(a2, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 / 0.3e1 + a1 / a2 * (scalar) mu * a3 / 0.3e1;

	Ke[11][12] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + a2 * (scalar) mu / 0.12e2;

	Ke[11][13] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 - a1 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 + a1 * (scalar) mu / 0.12e2;

	Ke[11][14] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.36e2;

	Ke[11][15] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) / 0.6e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1) * a3 * a3 / 0.2e1 + a2 * (scalar) mu / 0.12e2;

	Ke[11][16] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 / 0.6e1 - a1 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 + a1 * (scalar) mu / 0.6e1;

	Ke[11][17] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.18e2;

	Ke[11][18] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 * a3 / 0.6e1 + a2 * (scalar) mu / 0.6e1;

	Ke[11][19] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 * a3 / 0.4e1 + a1 * (scalar) mu / 0.12e2;

	Ke[11][20] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.18e2;

	Ke[11][21] = (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) * a3 * a3 / 0.12e2 + a2 * (scalar) mu / 0.6e1;

	Ke[11][22] = (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 * a3 * a3 / 0.12e2 + a1 * (scalar) mu / 0.6e1;

	Ke[11][23] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.9e1;

	Ke[12][0] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.9e1;

	Ke[12][1] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) lambda - pow(a3, -0.2e1) * (scalar) mu) * pow(a3, 0.3e1) / 0.3e1 + (((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda - pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + 0.1e1 / a3 * (scalar) lambda + 0.1e1 / a3 * (scalar) mu) * a3 * a3 / 0.2e1;

	Ke[12][2] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + (-pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) mu * a2 + ((0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) lambda * a2) * a3 * a3 / 0.2e1 - a2 * (scalar) mu / 0.6e1;

	Ke[12][3] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[12][4] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) mu) * pow(a3, 0.3e1) / 0.3e1 + (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda + pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - 0.1e1 / a3 * (scalar) mu) * a3 * a3 / 0.2e1;

	Ke[12][5] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((-0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) mu * a2) * a3 * a3 / 0.2e1 + a2 * (scalar) mu / 0.6e1;

	Ke[12][6] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[12][7] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 - 0.1e1 / a3 * (scalar) lambda / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[12][8] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 - a2 * (scalar) mu / 0.12e2;

	Ke[12][9] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.36e2;

	Ke[12][10] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 - 0.1e1 / a3 * (scalar) lambda / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[12][11] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + a2 * (scalar) mu / 0.12e2;

	Ke[12][12] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.9e1;

	Ke[12][13] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) lambda + pow(a3, -0.2e1) * (scalar) mu) * pow(a3, 0.3e1) / 0.3e1;

	Ke[12][14] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((-0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) mu * a2 - pow(a3, -0.2e1) * (scalar) lambda * a2) * a3 * a3 / 0.2e1;

	Ke[12][15] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[12][16] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) mu) * pow(a3, 0.3e1) / 0.3e1;

	Ke[12][17] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) mu * a2) * a3 * a3 / 0.2e1;

	Ke[12][18] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[12][19] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1;

	Ke[12][20] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[12][21] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.36e2;

	Ke[12][22] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1;

	Ke[12][23] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[13][0] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) lambda - pow(a3, -0.2e1) * (scalar) mu) * pow(a3, 0.3e1) / 0.3e1 + (((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda - pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + 0.1e1 / a3 * (scalar) lambda + 0.1e1 / a3 * (scalar) mu) * a3 * a3 / 0.2e1;

	Ke[13][1] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.9e1;

	Ke[13][2] = (((0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 - a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1 + (-0.2e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + a1 * pow(a3, -0.2e1) * (scalar) lambda + a1 * pow(a3, -0.2e1) * (scalar) mu) * a3 * a3 / 0.2e1 - a1 * (scalar) mu / 0.6e1;

	Ke[13][3] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) lambda) * pow(a3, 0.3e1) / 0.3e1 + (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) lambda) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda + pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - 0.1e1 / a3 * (scalar) lambda) * a3 * a3 / 0.2e1;

	Ke[13][4] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[13][5] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1 + (0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1 - a1 * (scalar) mu / 0.12e2;

	Ke[13][6] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 - 0.1e1 / a3 * (scalar) mu / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[13][7] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[13][8] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 - a1 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 + a1 * (scalar) mu / 0.6e1;

	Ke[13][9] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) lambda) * a2 * a2 / 0.2e1 - 0.1e1 / a3 * (scalar) mu / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[13][10] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.36e2;

	Ke[13][11] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 - a1 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 + a1 * (scalar) mu / 0.12e2;

	Ke[13][12] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) lambda + pow(a3, -0.2e1) * (scalar) mu) * pow(a3, 0.3e1) / 0.3e1;

	Ke[13][13] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.9e1;

	Ke[13][14] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1 + (0.2e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - a1 * pow(a3, -0.2e1) * (scalar) lambda - a1 * pow(a3, -0.2e1) * (scalar) mu) * a3 * a3 / 0.2e1;

	Ke[13][15] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) lambda) * pow(a3, 0.3e1) / 0.3e1;

	Ke[13][16] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[13][17] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1 + (-0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[13][18] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1;

	Ke[13][19] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[13][20] = (((0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 - a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + a1 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1;

	Ke[13][21] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1;

	Ke[13][22] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.36e2;

	Ke[13][23] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + a1 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1;

	Ke[14][0] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + (-pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) mu * a2 + ((0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) lambda * a2) * a3 * a3 / 0.2e1 - a2 * (scalar) lambda / 0.6e1;

	Ke[14][1] = (((0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 - a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1 + (-0.2e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + a1 * pow(a3, -0.2e1) * (scalar) lambda + a1 * pow(a3, -0.2e1) * (scalar) mu) * a3 * a3 / 0.2e1 - a1 * (scalar) lambda / 0.6e1;

	Ke[14][2] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.9e1;

	Ke[14][3] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) / 0.3e1 + ((-0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) lambda * a2) * a3 * a3 / 0.2e1 + a2 * (scalar) lambda / 0.6e1;

	Ke[14][4] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1 + (0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1 - a1 * (scalar) lambda / 0.12e2;

	Ke[14][5] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.18e2;

	Ke[14][6] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 - a2 * (scalar) lambda / 0.12e2;

	Ke[14][7] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 - a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.3e1) * a3 * a3 / 0.2e1 + a1 * (scalar) lambda / 0.6e1;

	Ke[14][8] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.18e2;

	Ke[14][9] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) / 0.3e1 + ((pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + a2 * (scalar) lambda / 0.12e2;

	Ke[14][10] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 - a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.6e1) * a3 * a3 / 0.2e1 + a1 * (scalar) lambda / 0.12e2;

	Ke[14][11] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.36e2;

	Ke[14][12] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((-0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) mu * a2 - pow(a3, -0.2e1) * (scalar) lambda * a2) * a3 * a3 / 0.2e1;

	Ke[14][13] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1 + (0.2e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - a1 * pow(a3, -0.2e1) * (scalar) lambda - a1 * pow(a3, -0.2e1) * (scalar) mu) * a3 * a3 / 0.2e1;

	Ke[14][14] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.9e1;

	Ke[14][15] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) / 0.3e1 + ((0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) lambda * a2) * a3 * a3 / 0.2e1;

	Ke[14][16] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1 + (-0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[14][17] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.18e2;

	Ke[14][18] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[14][19] = (((0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 - a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.3e1) * a3 * a3 / 0.2e1;

	Ke[14][20] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.18e2;

	Ke[14][21] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[14][22] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.6e1) * a3 * a3 / 0.2e1;

	Ke[14][23] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.36e2;

	Ke[15][0] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[15][1] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) lambda) * pow(a3, 0.3e1) / 0.3e1 + (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) lambda) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda + pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - 0.1e1 / a3 * (scalar) lambda) * a3 * a3 / 0.2e1;

	Ke[15][2] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) / 0.3e1 + ((-0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) lambda * a2) * a3 * a3 / 0.2e1 - a2 * (scalar) mu / 0.6e1;

	Ke[15][3] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.9e1;

	Ke[15][4] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + (-pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda - pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[15][5] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) / 0.6e1 + (0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + (-pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + a2 * (scalar) mu / 0.6e1;

	Ke[15][6] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.36e2;

	Ke[15][7] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) lambda) * a2 * a2 / 0.2e1 + 0.1e1 / a3 * (scalar) lambda / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[15][8] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) / 0.3e1 + ((pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 - a2 * (scalar) mu / 0.12e2;

	Ke[15][9] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[15][10] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 - pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + 0.1e1 / a3 * (scalar) lambda / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[15][11] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) / 0.6e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1) * a3 * a3 / 0.2e1 + a2 * (scalar) mu / 0.12e2;

	Ke[15][12] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[15][13] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) lambda) * pow(a3, 0.3e1) / 0.3e1;

	Ke[15][14] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) / 0.3e1 + ((0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) lambda * a2) * a3 * a3 / 0.2e1;

	Ke[15][15] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.9e1;

	Ke[15][16] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * pow(a3, 0.3e1) / 0.3e1;

	Ke[15][17] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) / 0.6e1 + (-0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + (pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[15][18] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.36e2;

	Ke[15][19] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1;

	Ke[15][20] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[15][21] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[15][22] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1;

	Ke[15][23] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) / 0.6e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1) * a3 * a3 / 0.2e1;

	Ke[16][0] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) mu) * pow(a3, 0.3e1) / 0.3e1 + (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda + pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - 0.1e1 / a3 * (scalar) mu) * a3 * a3 / 0.2e1;

	Ke[16][1] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[16][2] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1 + (0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1 - a1 * (scalar) mu / 0.12e2;

	Ke[16][3] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + (-pow(a1, -0.2e1) / a2 / a3 * (scalar) lambda - pow(a1, -0.2e1) / a2 / a3 * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[16][4] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.9e1;

	Ke[16][5] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 / 0.6e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1) * a3 * a3 / 0.2e1 - a1 * (scalar) mu / 0.6e1;

	Ke[16][6] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 + 0.1e1 / a3 * (scalar) mu / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[16][7] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.36e2;

	Ke[16][8] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 - a1 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 + a1 * (scalar) mu / 0.12e2;

	Ke[16][9] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 - pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + 0.1e1 / a3 * (scalar) mu / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[16][10] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[16][11] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 / 0.6e1 - a1 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 + a1 * (scalar) mu / 0.6e1;

	Ke[16][12] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) mu) * pow(a3, 0.3e1) / 0.3e1;

	Ke[16][13] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[16][14] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1 + (-0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[16][15] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * pow(a3, 0.3e1) / 0.3e1;

	Ke[16][16] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.9e1;

	Ke[16][17] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 / 0.6e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1) * a3 * a3 / 0.2e1;

	Ke[16][18] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1;

	Ke[16][19] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.36e2;

	Ke[16][20] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + a1 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1;

	Ke[16][21] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1;

	Ke[16][22] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[16][23] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 / 0.6e1 + a1 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1;

	Ke[17][0] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((-0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) mu * a2) * a3 * a3 / 0.2e1 - a2 * (scalar) lambda / 0.6e1;

	Ke[17][1] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1 + (0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1 - a1 * (scalar) lambda / 0.12e2;

	Ke[17][2] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.18e2;

	Ke[17][3] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) / 0.6e1 + (0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + (-pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + a2 * (scalar) lambda / 0.6e1;

	Ke[17][4] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 / 0.6e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1) * a3 * a3 / 0.2e1 - a1 * (scalar) lambda / 0.6e1;

	Ke[17][5] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.9e1;

	Ke[17][6] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 - a2 * (scalar) lambda / 0.12e2;

	Ke[17][7] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 - a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.6e1) * a3 * a3 / 0.2e1 + a1 * (scalar) lambda / 0.12e2;

	Ke[17][8] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.36e2;

	Ke[17][9] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) / 0.6e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1) * a3 * a3 / 0.2e1 + a2 * (scalar) lambda / 0.12e2;

	Ke[17][10] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 / 0.6e1 - a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.3e1) * a3 * a3 / 0.2e1 + a1 * (scalar) lambda / 0.6e1;

	Ke[17][11] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.18e2;

	Ke[17][12] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.2e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + (-pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) mu * a2) * a3 * a3 / 0.2e1;

	Ke[17][13] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1 + (-0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[17][14] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.18e2;

	Ke[17][15] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) / 0.6e1 + (-0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 * pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + (pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[17][16] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 / 0.6e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 / 0.3e1) * a3 * a3 / 0.2e1;

	Ke[17][17] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.9e1;

	Ke[17][18] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[17][19] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.6e1) * a3 * a3 / 0.2e1;

	Ke[17][20] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.36e2;

	Ke[17][21] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) / 0.6e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1) * a3 * a3 / 0.2e1;

	Ke[17][22] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 / 0.6e1 + a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.3e1) * a3 * a3 / 0.2e1;

	Ke[17][23] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.18e2;

	Ke[18][0] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[18][1] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 - 0.1e1 / a3 * (scalar) mu / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[18][2] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 - a2 * (scalar) mu / 0.12e2;

	Ke[18][3] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.36e2;

	Ke[18][4] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 + 0.1e1 / a3 * (scalar) mu / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[18][5] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + a2 * (scalar) mu / 0.12e2;

	Ke[18][6] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.9e1;

	Ke[18][7] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 * pow(a3, 0.3e1) / 0.6e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 * a3 * a3 / 0.4e1;

	Ke[18][8] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 * a3 / 0.6e1 - a2 * (scalar) mu / 0.6e1;

	Ke[18][9] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[18][10] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 * pow(a3, 0.3e1) / 0.6e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 * a3 * a3 / 0.4e1;

	Ke[18][11] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 * a3 / 0.6e1 + a2 * (scalar) mu / 0.6e1;

	Ke[18][12] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[18][13] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1;

	Ke[18][14] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[18][15] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.36e2;

	Ke[18][16] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1;

	Ke[18][17] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[18][18] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.9e1;

	Ke[18][19] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 * pow(a3, 0.3e1) / 0.6e1;

	Ke[18][20] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 * a3 / 0.6e1;

	Ke[18][21] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[18][22] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 * pow(a3, 0.3e1) / 0.6e1;

	Ke[18][23] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 * a3 / 0.6e1;

	Ke[19][0] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 - 0.1e1 / a3 * (scalar) lambda / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[19][1] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[19][2] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 - a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.3e1) * a3 * a3 / 0.2e1 - a1 * (scalar) mu / 0.6e1;

	Ke[19][3] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) lambda) * a2 * a2 / 0.2e1 + 0.1e1 / a3 * (scalar) lambda / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[19][4] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.36e2;

	Ke[19][5] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 - a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.6e1) * a3 * a3 / 0.2e1 - a1 * (scalar) mu / 0.12e2;

	Ke[19][6] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 * pow(a3, 0.3e1) / 0.6e1 + ((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 * a3 * a3 / 0.4e1;

	Ke[19][7] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.9e1;

	Ke[19][8] = ((0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 - a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 * a3 * a3 / 0.4e1 + a1 * (scalar) mu / 0.6e1;

	Ke[19][9] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 * pow(a3, 0.3e1) / 0.6e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) lambda) * a2 * a2 * a3 * a3 / 0.4e1;

	Ke[19][10] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[19][11] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 * a3 / 0.4e1 + a1 * (scalar) mu / 0.12e2;

	Ke[19][12] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1;

	Ke[19][13] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[19][14] = (((0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 - a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.3e1) * a3 * a3 / 0.2e1;

	Ke[19][15] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1;

	Ke[19][16] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.36e2;

	Ke[19][17] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.6e1) * a3 * a3 / 0.2e1;

	Ke[19][18] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 * pow(a3, 0.3e1) / 0.6e1;

	Ke[19][19] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.9e1;

	Ke[19][20] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 * a3 * a3 / 0.4e1;

	Ke[19][21] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 * pow(a3, 0.3e1) / 0.6e1;

	Ke[19][22] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[19][23] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 * a3 / 0.4e1;

	Ke[20][0] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 - a2 * (scalar) lambda / 0.12e2;

	Ke[20][1] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 - a1 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * (scalar) lambda / 0.6e1;

	Ke[20][2] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.18e2;

	Ke[20][3] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) / 0.3e1 + ((pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 + a2 * (scalar) lambda / 0.12e2;

	Ke[20][4] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 - a1 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * (scalar) lambda / 0.12e2;

	Ke[20][5] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.36e2;

	Ke[20][6] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 * a3 / 0.6e1 - a2 * (scalar) lambda / 0.6e1;

	Ke[20][7] = ((0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 - a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 * a3 * a3 / 0.4e1 + a1 * (scalar) lambda / 0.6e1;

	Ke[20][8] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.9e1;

	Ke[20][9] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) * a3 * a3 / 0.6e1 + a2 * (scalar) lambda / 0.6e1;

	Ke[20][10] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 * a3 / 0.4e1 + a1 * (scalar) lambda / 0.12e2;

	Ke[20][11] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.18e2;

	Ke[20][12] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[20][13] = (((0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 - a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + a1 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1;

	Ke[20][14] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.18e2;

	Ke[20][15] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[20][16] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + a1 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1;

	Ke[20][17] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.36e2;

	Ke[20][18] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 * a3 / 0.6e1;

	Ke[20][19] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.2e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda * a1 + a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 * a3 * a3 / 0.4e1;

	Ke[20][20] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.9e1;

	Ke[20][21] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) * a3 * a3 / 0.6e1;

	Ke[20][22] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 * a3 / 0.4e1;

	Ke[20][23] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.18e2;

	Ke[21][0] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.36e2;

	Ke[21][1] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) lambda) * a2 * a2 / 0.2e1 - 0.1e1 / a3 * (scalar) mu / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[21][2] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) / 0.3e1 + ((pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 - a2 * (scalar) mu / 0.12e2;

	Ke[21][3] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[21][4] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 - pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + 0.1e1 / a3 * (scalar) mu / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[21][5] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) / 0.6e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1) * a3 * a3 / 0.2e1 + a2 * (scalar) mu / 0.12e2;

	Ke[21][6] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[21][7] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 * pow(a3, 0.3e1) / 0.6e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) lambda) * a2 * a2 * a3 * a3 / 0.4e1;

	Ke[21][8] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) * a3 * a3 / 0.6e1 - a2 * (scalar) mu / 0.6e1;

	Ke[21][9] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.9e1;

	Ke[21][10] = (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 * pow(a3, 0.3e1) / 0.12e2 + (pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a2 * a3 * a3 / 0.8e1;

	Ke[21][11] = (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) * a3 * a3 / 0.12e2 + a2 * (scalar) mu / 0.6e1;

	Ke[21][12] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.36e2;

	Ke[21][13] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1;

	Ke[21][14] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[21][15] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[21][16] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + pow(a3, -0.2e1) * (scalar) mu / 0.2e1) * pow(a3, 0.3e1) / 0.3e1;

	Ke[21][17] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) / 0.6e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1) * a3 * a3 / 0.2e1;

	Ke[21][18] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[21][19] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * a2 * a2 * pow(a3, 0.3e1) / 0.6e1;

	Ke[21][20] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda) * pow(a2, 0.3e1) * a3 * a3 / 0.6e1;

	Ke[21][21] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.9e1;

	Ke[21][22] = (pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 * pow(a3, 0.3e1) / 0.12e2;

	Ke[21][23] = (pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) * a3 * a3 / 0.12e2;

	Ke[22][0] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 + pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + (((pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 / 0.2e1 - 0.1e1 / a3 * (scalar) lambda / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[22][1] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.36e2;

	Ke[22][2] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 - a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.6e1) * a3 * a3 / 0.2e1 - a1 * (scalar) mu / 0.12e2;

	Ke[22][3] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 - pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + 0.1e1 / a3 * (scalar) lambda / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[22][4] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[22][5] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 / 0.6e1 - a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.3e1) * a3 * a3 / 0.2e1 - a1 * (scalar) mu / 0.6e1;

	Ke[22][6] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 * pow(a3, 0.3e1) / 0.6e1 + ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) / a3 * (scalar) mu) * a2 * a2 * a3 * a3 / 0.4e1;

	Ke[22][7] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[22][8] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 * a3 / 0.4e1 + a1 * (scalar) mu / 0.12e2;

	Ke[22][9] = (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 * pow(a3, 0.3e1) / 0.12e2 + (pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) / a3 * (scalar) mu) * a1 * a1 * a2 * a2 * a3 * a3 / 0.8e1;

	Ke[22][10] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) (lambda + 2 * mu) / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) mu / 0.9e1;

	Ke[22][11] = (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 * a3 * a3 / 0.12e2 + a1 * (scalar) mu / 0.6e1;

	Ke[22][12] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1 - pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1;

	Ke[22][13] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.36e2;

	Ke[22][14] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.6e1) * a3 * a3 / 0.2e1;

	Ke[22][15] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1 + pow(a3, -0.2e1) * (scalar) lambda / 0.2e1) * pow(a3, 0.3e1) / 0.3e1;

	Ke[22][16] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[22][17] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 / 0.6e1 + a1 * pow(a3, -0.2e1) * (scalar) lambda / 0.3e1) * a3 * a3 / 0.2e1;

	Ke[22][18] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 * pow(a3, 0.3e1) / 0.6e1;

	Ke[22][19] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.18e2;

	Ke[22][20] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 * a3 / 0.4e1;

	Ke[22][21] = (pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 * pow(a3, 0.3e1) / 0.12e2;

	Ke[22][22] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) (lambda + 2 * mu) / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) mu / 0.9e1;

	Ke[22][23] = (pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 * a3 * a3 / 0.12e2;

	Ke[23][0] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1 - a2 * (scalar) lambda / 0.12e2;

	Ke[23][1] = (((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 - a1 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * (scalar) lambda / 0.12e2;

	Ke[23][2] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.36e2;

	Ke[23][3] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) / 0.6e1 + (-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1) * a3 * a3 / 0.2e1 + a2 * (scalar) lambda / 0.12e2;

	Ke[23][4] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 / 0.6e1 - a1 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * (scalar) lambda / 0.6e1;

	Ke[23][5] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.6e1 - a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.18e2;

	Ke[23][6] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 * a3 / 0.6e1 - a2 * (scalar) lambda / 0.6e1;

	Ke[23][7] = ((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 * a3 / 0.4e1 + a1 * (scalar) lambda / 0.12e2;

	Ke[23][8] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.18e2;

	Ke[23][9] = (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) * a3 * a3 / 0.12e2 + a2 * (scalar) lambda / 0.6e1;

	Ke[23][10] = (-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 * a3 * a3 / 0.12e2 + a1 * (scalar) lambda / 0.6e1;

	Ke[23][11] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + (0.1e1 / a1 * a2 / a3 * (scalar) mu / 0.3e1 + a1 / a2 / a3 * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1 - a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.9e1;

	Ke[23][12] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 - pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) / 0.3e1 + ((-pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + 0.1e1 / a2 * pow(a3, -0.2e1) * (scalar) mu) * a2 * a2 / 0.2e1) * a3 * a3 / 0.2e1;

	Ke[23][13] = (((pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (-0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 / 0.2e1 + a1 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * a3 * a3 / 0.2e1;

	Ke[23][14] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.36e2;

	Ke[23][15] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) / 0.6e1 + (pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) / a2 * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * a2 * a2 / 0.4e1) * a3 * a3 / 0.2e1;

	Ke[23][16] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 / 0.6e1 + a1 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * a3 * a3 / 0.2e1;

	Ke[23][17] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1 - a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.18e2;

	Ke[23][18] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1 + pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a2, 0.3e1) * a3 * a3 / 0.6e1;

	Ke[23][19] = ((-pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda - pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) / 0.3e1 + (0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + 0.1e1 / a1 * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 / 0.2e1) * a2 * a2 * a3 * a3 / 0.4e1;

	Ke[23][20] = (-0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.6e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.18e2;

	Ke[23][21] = (pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * a1 * a1 * pow(a2, 0.3e1) * a3 * a3 / 0.12e2;

	Ke[23][22] = (pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) lambda + pow(a1, -0.2e1) * pow(a2, -0.2e1) * pow(a3, -0.2e1) * (scalar) mu) * pow(a1, 0.3e1) * a2 * a2 * a3 * a3 / 0.12e2;

	Ke[23][23] = (0.1e1 / a1 * a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1 + a1 / a2 * pow(a3, -0.2e1) * (scalar) mu / 0.3e1) * pow(a3, 0.3e1) / 0.3e1 + a1 * a2 / a3 * (scalar) (lambda + 2 * mu) / 0.9e1;

}

template<typename scalar>
void ComputeBavg(scalar a1, scalar a2, scalar a3, scalar Bavg[6][24])
{
	memset(Bavg, 0, sizeof(scalar)*6*24);
	Bavg[0][0] = -0.1e1 / a1 / 0.4e1;
	Bavg[0][1] = 0;
	Bavg[0][2] = 0;
	Bavg[0][3] = 0.1e1 / a1 / 0.4e1;
	Bavg[0][4] = 0;
	Bavg[0][5] = 0;
	Bavg[0][6] = -0.1e1 / a1 / 0.4e1;
	Bavg[0][7] = 0;
	Bavg[0][8] = 0;
	Bavg[0][9] = 0.1e1 / a1 / 0.4e1;
	Bavg[0][10] = 0;
	Bavg[0][11] = 0;
	Bavg[0][12] = -0.1e1 / a1 / 0.4e1;
	Bavg[0][13] = 0;
	Bavg[0][14] = 0;
	Bavg[0][15] = 0.1e1 / a1 / 0.4e1;
	Bavg[0][16] = 0;
	Bavg[0][17] = 0;
	Bavg[0][18] = -0.1e1 / a1 / 0.4e1;
	Bavg[0][19] = 0;
	Bavg[0][20] = 0;
	Bavg[0][21] = 0.1e1 / a1 / 0.4e1;
	Bavg[0][22] = 0;
	Bavg[0][23] = 0;
	Bavg[1][0] = 0;
	Bavg[1][1] = -0.1e1 / a2 / 0.4e1;
	Bavg[1][2] = 0;
	Bavg[1][3] = 0;
	Bavg[1][4] = -0.1e1 / a2 / 0.4e1;
	Bavg[1][5] = 0;
	Bavg[1][6] = 0;
	Bavg[1][7] = 0.1e1 / a2 / 0.4e1;
	Bavg[1][8] = 0;
	Bavg[1][9] = 0;
	Bavg[1][10] = 0.1e1 / a2 / 0.4e1;
	Bavg[1][11] = 0;
	Bavg[1][12] = 0;
	Bavg[1][13] = -0.1e1 / a2 / 0.4e1;
	Bavg[1][14] = 0;
	Bavg[1][15] = 0;
	Bavg[1][16] = -0.1e1 / a2 / 0.4e1;
	Bavg[1][17] = 0;
	Bavg[1][18] = 0;
	Bavg[1][19] = 0.1e1 / a2 / 0.4e1;
	Bavg[1][20] = 0;
	Bavg[1][21] = 0;
	Bavg[1][22] = 0.1e1 / a2 / 0.4e1;
	Bavg[1][23] = 0;
	Bavg[2][0] = 0;
	Bavg[2][1] = 0;
	Bavg[2][2] = -0.1e1 / a3 / 0.4e1;
	Bavg[2][3] = 0;
	Bavg[2][4] = 0;
	Bavg[2][5] = -0.1e1 / a3 / 0.4e1;
	Bavg[2][6] = 0;
	Bavg[2][7] = 0;
	Bavg[2][8] = -0.1e1 / a3 / 0.4e1;
	Bavg[2][9] = 0;
	Bavg[2][10] = 0;
	Bavg[2][11] = -0.1e1 / a3 / 0.4e1;
	Bavg[2][12] = 0;
	Bavg[2][13] = 0;
	Bavg[2][14] = 0.1e1 / a3 / 0.4e1;
	Bavg[2][15] = 0;
	Bavg[2][16] = 0;
	Bavg[2][17] = 0.1e1 / a3 / 0.4e1;
	Bavg[2][18] = 0;
	Bavg[2][19] = 0;
	Bavg[2][20] = 0.1e1 / a3 / 0.4e1;
	Bavg[2][21] = 0;
	Bavg[2][22] = 0;
	Bavg[2][23] = 0.1e1 / a3 / 0.4e1;
	Bavg[3][0] = -0.1e1 / a2 / 0.8e1;
	Bavg[3][1] = -0.1e1 / a1 / 0.8e1;
	Bavg[3][2] = 0;
	Bavg[3][3] = -0.1e1 / a2 / 0.8e1;
	Bavg[3][4] = 0.1e1 / a1 / 0.8e1;
	Bavg[3][5] = 0;
	Bavg[3][6] = 0.1e1 / a2 / 0.8e1;
	Bavg[3][7] = -0.1e1 / a1 / 0.8e1;
	Bavg[3][8] = 0;
	Bavg[3][9] = 0.1e1 / a2 / 0.8e1;
	Bavg[3][10] = 0.1e1 / a1 / 0.8e1;
	Bavg[3][11] = 0;
	Bavg[3][12] = -0.1e1 / a2 / 0.8e1;
	Bavg[3][13] = -0.1e1 / a1 / 0.8e1;
	Bavg[3][14] = 0;
	Bavg[3][15] = -0.1e1 / a2 / 0.8e1;
	Bavg[3][16] = 0.1e1 / a1 / 0.8e1;
	Bavg[3][17] = 0;
	Bavg[3][18] = 0.1e1 / a2 / 0.8e1;
	Bavg[3][19] = -0.1e1 / a1 / 0.8e1;
	Bavg[3][20] = 0;
	Bavg[3][21] = 0.1e1 / a2 / 0.8e1;
	Bavg[3][22] = 0.1e1 / a1 / 0.8e1;
	Bavg[3][23] = 0;
	Bavg[4][0] = -0.1e1 / a3 / 0.8e1;
	Bavg[4][1] = 0;
	Bavg[4][2] = -0.1e1 / a1 / 0.8e1;
	Bavg[4][3] = -0.1e1 / a3 / 0.8e1;
	Bavg[4][4] = 0;
	Bavg[4][5] = 0.1e1 / a1 / 0.8e1;
	Bavg[4][6] = -0.1e1 / a3 / 0.8e1;
	Bavg[4][7] = 0;
	Bavg[4][8] = -0.1e1 / a1 / 0.8e1;
	Bavg[4][9] = -0.1e1 / a3 / 0.8e1;
	Bavg[4][10] = 0;
	Bavg[4][11] = 0.1e1 / a1 / 0.8e1;
	Bavg[4][12] = 0.1e1 / a3 / 0.8e1;
	Bavg[4][13] = 0;
	Bavg[4][14] = -0.1e1 / a1 / 0.8e1;
	Bavg[4][15] = 0.1e1 / a3 / 0.8e1;
	Bavg[4][16] = 0;
	Bavg[4][17] = 0.1e1 / a1 / 0.8e1;
	Bavg[4][18] = 0.1e1 / a3 / 0.8e1;
	Bavg[4][19] = 0;
	Bavg[4][20] = -0.1e1 / a1 / 0.8e1;
	Bavg[4][21] = 0.1e1 / a3 / 0.8e1;
	Bavg[4][22] = 0;
	Bavg[4][23] = 0.1e1 / a1 / 0.8e1;
	Bavg[5][0] = 0;
	Bavg[5][1] = -0.1e1 / a3 / 0.8e1;
	Bavg[5][2] = -0.1e1 / a2 / 0.8e1;
	Bavg[5][3] = 0;
	Bavg[5][4] = -0.1e1 / a3 / 0.8e1;
	Bavg[5][5] = -0.1e1 / a2 / 0.8e1;
	Bavg[5][6] = 0;
	Bavg[5][7] = -0.1e1 / a3 / 0.8e1;
	Bavg[5][8] = 0.1e1 / a2 / 0.8e1;
	Bavg[5][9] = 0;
	Bavg[5][10] = -0.1e1 / a3 / 0.8e1;
	Bavg[5][11] = 0.1e1 / a2 / 0.8e1;
	Bavg[5][12] = 0;
	Bavg[5][13] = 0.1e1 / a3 / 0.8e1;
	Bavg[5][14] = -0.1e1 / a2 / 0.8e1;
	Bavg[5][15] = 0;
	Bavg[5][16] = 0.1e1 / a3 / 0.8e1;
	Bavg[5][17] = -0.1e1 / a2 / 0.8e1;
	Bavg[5][18] = 0;
	Bavg[5][19] = 0.1e1 / a3 / 0.8e1;
	Bavg[5][20] = 0.1e1 / a2 / 0.8e1;
	Bavg[5][21] = 0;
	Bavg[5][22] = 0.1e1 / a3 / 0.8e1;
	Bavg[5][23] = 0.1e1 / a2 / 0.8e1;
}

template<typename scalar>
void ComputeS(scalar a1, scalar a2, scalar a3, scalar E, scalar nu, scalar S[6][24])
{
	memset(S, 0, sizeof(scalar)*6*24);

	scalar lambda = E*nu/((1+nu)*(1-2*nu));
	scalar mu = E/(2*(1+nu));

	S[0][0] = -(scalar) ((lambda + 2 * mu) / a1) / 0.4e1;

	S[0][1] = -(scalar) (lambda / a2) / 0.4e1;

	S[0][2] = -(scalar) (1 / a3 * lambda) / 0.4e1;

	S[0][3] = (scalar) ((lambda + 2 * mu) / a1) / 0.4e1;

	S[0][4] = -(scalar) (lambda / a2) / 0.4e1;

	S[0][5] = -(scalar) (1 / a3 * lambda) / 0.4e1;

	S[0][6] = -(scalar) ((lambda + 2 * mu) / a1) / 0.4e1;

	S[0][7] = (scalar) (lambda / a2) / 0.4e1;

	S[0][8] = -(scalar) (1 / a3 * lambda) / 0.4e1;

	S[0][9] = (scalar) ((lambda + 2 * mu) / a1) / 0.4e1;

	S[0][10] = (scalar) (lambda / a2) / 0.4e1;

	S[0][11] = -(scalar) (1 / a3 * lambda) / 0.4e1;

	S[0][12] = -(scalar) ((lambda + 2 * mu) / a1) / 0.4e1;

	S[0][13] = -(scalar) (lambda / a2) / 0.4e1;

	S[0][14] = (scalar) (1 / a3 * lambda) / 0.4e1;

	S[0][15] = (scalar) ((lambda + 2 * mu) / a1) / 0.4e1;

	S[0][16] = -(scalar) (lambda / a2) / 0.4e1;

	S[0][17] = (scalar) (1 / a3 * lambda) / 0.4e1;

	S[0][18] = -(scalar) ((lambda + 2 * mu) / a1) / 0.4e1;

	S[0][19] = (scalar) (lambda / a2) / 0.4e1;

	S[0][20] = (scalar) (1 / a3 * lambda) / 0.4e1;

	S[0][21] = (scalar) ((lambda + 2 * mu) / a1) / 0.4e1;

	S[0][22] = (scalar) (lambda / a2) / 0.4e1;

	S[0][23] = (scalar) (1 / a3 * lambda) / 0.4e1;

	S[1][0] = -(scalar) (lambda / a1) / 0.4e1;

	S[1][1] = -(scalar) ((lambda + 2 * mu) / a2) / 0.4e1;

	S[1][2] = -(scalar) (1 / a3 * lambda) / 0.4e1;

	S[1][3] = (scalar) (lambda / a1) / 0.4e1;

	S[1][4] = -(scalar) ((lambda + 2 * mu) / a2) / 0.4e1;

	S[1][5] = -(scalar) (1 / a3 * lambda) / 0.4e1;

	S[1][6] = -(scalar) (lambda / a1) / 0.4e1;

	S[1][7] = (scalar) ((lambda + 2 * mu) / a2) / 0.4e1;

	S[1][8] = -(scalar) (1 / a3 * lambda) / 0.4e1;

	S[1][9] = (scalar) (lambda / a1) / 0.4e1;

	S[1][10] = (scalar) ((lambda + 2 * mu) / a2) / 0.4e1;

	S[1][11] = -(scalar) (1 / a3 * lambda) / 0.4e1;

	S[1][12] = -(scalar) (lambda / a1) / 0.4e1;

	S[1][13] = -(scalar) ((lambda + 2 * mu) / a2) / 0.4e1;

	S[1][14] = (scalar) (1 / a3 * lambda) / 0.4e1;

	S[1][15] = (scalar) (lambda / a1) / 0.4e1;

	S[1][16] = -(scalar) ((lambda + 2 * mu) / a2) / 0.4e1;

	S[1][17] = (scalar) (1 / a3 * lambda) / 0.4e1;

	S[1][18] = -(scalar) (lambda / a1) / 0.4e1;

	S[1][19] = (scalar) ((lambda + 2 * mu) / a2) / 0.4e1;

	S[1][20] = (scalar) (1 / a3 * lambda) / 0.4e1;

	S[1][21] = (scalar) (lambda / a1) / 0.4e1;

	S[1][22] = (scalar) ((lambda + 2 * mu) / a2) / 0.4e1;

	S[1][23] = (scalar) (1 / a3 * lambda) / 0.4e1;

	S[2][0] = -(scalar) (lambda / a1) / 0.4e1;

	S[2][1] = -(scalar) (lambda / a2) / 0.4e1;

	S[2][2] = -(scalar) ((lambda + 2 * mu) / a3) / 0.4e1;

	S[2][3] = (scalar) (lambda / a1) / 0.4e1;

	S[2][4] = -(scalar) (lambda / a2) / 0.4e1;

	S[2][5] = -(scalar) ((lambda + 2 * mu) / a3) / 0.4e1;

	S[2][6] = -(scalar) (lambda / a1) / 0.4e1;

	S[2][7] = (scalar) (lambda / a2) / 0.4e1;

	S[2][8] = -(scalar) ((lambda + 2 * mu) / a3) / 0.4e1;

	S[2][9] = (scalar) (lambda / a1) / 0.4e1;

	S[2][10] = (scalar) (lambda / a2) / 0.4e1;

	S[2][11] = -(scalar) ((lambda + 2 * mu) / a3) / 0.4e1;

	S[2][12] = -(scalar) (lambda / a1) / 0.4e1;

	S[2][13] = -(scalar) (lambda / a2) / 0.4e1;

	S[2][14] = (scalar) ((lambda + 2 * mu) / a3) / 0.4e1;

	S[2][15] = (scalar) (lambda / a1) / 0.4e1;

	S[2][16] = -(scalar) (lambda / a2) / 0.4e1;

	S[2][17] = (scalar) ((lambda + 2 * mu) / a3) / 0.4e1;

	S[2][18] = -(scalar) (lambda / a1) / 0.4e1;

	S[2][19] = (scalar) (lambda / a2) / 0.4e1;

	S[2][20] = (scalar) ((lambda + 2 * mu) / a3) / 0.4e1;

	S[2][21] = (scalar) (lambda / a1) / 0.4e1;

	S[2][22] = (scalar) (lambda / a2) / 0.4e1;

	S[2][23] = (scalar) ((lambda + 2 * mu) / a3) / 0.4e1;

	S[3][0] = -(scalar) (mu / a2) / 0.2e1;

	S[3][1] = -(scalar) (mu / a1) / 0.2e1;

	S[3][3] = -(scalar) (mu / a2) / 0.2e1;

	S[3][4] = (scalar) (mu / a1) / 0.2e1;

	S[3][6] = (scalar) (mu / a2) / 0.2e1;

	S[3][7] = -(scalar) (mu / a1) / 0.2e1;

	S[3][9] = (scalar) (mu / a2) / 0.2e1;

	S[3][10] = (scalar) (mu / a1) / 0.2e1;

	S[3][12] = -(scalar) (mu / a2) / 0.2e1;

	S[3][13] = -(scalar) (mu / a1) / 0.2e1;

	S[3][15] = -(scalar) (mu / a2) / 0.2e1;

	S[3][16] = (scalar) (mu / a1) / 0.2e1;

	S[3][18] = (scalar) (mu / a2) / 0.2e1;

	S[3][19] = -(scalar) (mu / a1) / 0.2e1;

	S[3][21] = (scalar) (mu / a2) / 0.2e1;

	S[3][22] = (scalar) (mu / a1) / 0.2e1;

	S[4][0] = -(scalar) (1 / a3 * mu) / 0.2e1;

	S[4][2] = -(scalar) (mu / a1) / 0.2e1;

	S[4][3] = -(scalar) (1 / a3 * mu) / 0.2e1;

	S[4][5] = (scalar) (mu / a1) / 0.2e1;

	S[4][6] = -(scalar) (1 / a3 * mu) / 0.2e1;

	S[4][8] = -(scalar) (mu / a1) / 0.2e1;

	S[4][9] = -(scalar) (1 / a3 * mu) / 0.2e1;

	S[4][11] = (scalar) (mu / a1) / 0.2e1;

	S[4][12] = (scalar) (1 / a3 * mu) / 0.2e1;

	S[4][14] = -(scalar) (mu / a1) / 0.2e1;

	S[4][15] = (scalar) (1 / a3 * mu) / 0.2e1;

	S[4][17] = (scalar) (mu / a1) / 0.2e1;

	S[4][18] = (scalar) (1 / a3 * mu) / 0.2e1;

	S[4][20] = -(scalar) (mu / a1) / 0.2e1;

	S[4][21] = (scalar) (1 / a3 * mu) / 0.2e1;

	S[4][23] = (scalar) (mu / a1) / 0.2e1;

	S[5][1] = -(scalar) (1 / a3 * mu) / 0.2e1;

	S[5][2] = -(scalar) (mu / a2) / 0.2e1;

	S[5][4] = -(scalar) (1 / a3 * mu) / 0.2e1;

	S[5][5] = -(scalar) (mu / a2) / 0.2e1;

	S[5][7] = -(scalar) (1 / a3 * mu) / 0.2e1;

	S[5][8] = (scalar) (mu / a2) / 0.2e1;

	S[5][10] = -(scalar) (1 / a3 * mu) / 0.2e1;

	S[5][11] = (scalar) (mu / a2) / 0.2e1;

	S[5][13] = (scalar) (1 / a3 * mu) / 0.2e1;

	S[5][14] = -(scalar) (mu / a2) / 0.2e1;

	S[5][16] = (scalar) (1 / a3 * mu) / 0.2e1;

	S[5][17] = -(scalar) (mu / a2) / 0.2e1;

	S[5][19] = (scalar) (1 / a3 * mu) / 0.2e1;

	S[5][20] = (scalar) (mu / a2) / 0.2e1;

	S[5][22] = (scalar) (1 / a3 * mu) / 0.2e1;

	S[5][23] = (scalar) (mu / a2) / 0.2e1;

}

template void ComputeKe<double>(double a1, double a2, double a3, double E, double nu, double Ke[24][24]);
template void ComputeKe<float>(float a1, float a2, float a3, float E, float nu, float Ke[24][24]);

template void ComputeBavg<double>(double a1, double a2, double a3, double Bavg[6][24]);
template void ComputeBavg<float>(float a1, float a2, float a3, float Bavg[6][24]);

template void ComputeS<double>(double a1, double a2, double a3, double E, double nu, double S[6][24]);
template void ComputeS<float>(float a1, float a2, float a3, float E, float nu, float S[6][24]);

#endif // __TUM3D__KE_H__
