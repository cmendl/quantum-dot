/// \file density_int.c
/// \brief Calculate integral Ne(r) of the density rho(r), as well as the inverse function of Ne.
//
//	Copyright (c) 2014, Christian B. Mendl
//	All rights reserved.
//	http://christian.mendl.net
//
//	This program is free software; you can redistribute it and/or
//	modify it under the terms of the Simplified BSD License
//	http://www.opensource.org/licenses/bsd-license.php
//
//	Reference:
//	  Christian B. Mendl, Francesc Malet, Paola Gori-Giorgi
//	  Wigner localization in quantum dots from Kohn-Sham density functional theory without symmetry breaking
//	  Physical Review B 89, 125106 (2014)
//	  (preprint http://arxiv.org/abs/1311.6011)
//_______________________________________________________________________________________________________________________
//

#include "density_int.h"
#include "util.h"
#include "error_desc.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>


//_______________________________________________________________________________________________________________________
///
/// \brief Evaluate N_e^{-1}(r)
///
static inline double EvaluateInvNe(const double dr, const int ngrid, const double *Ne, const int *invtab, const double d)
{
	int i;

	if (d < 0)
	{
		return 0;
	}
	if (d >= 1)
	{
		return *(double *)&DBL_INF;
	}

	// index of uniform d grid
	i = (int)floor(d * ngrid);
	assert(i < ngrid);

	i = invtab[i];
	while (Ne[i+1] <= d)
	{
		i++;
		assert(i < ngrid-1);
	}
	assert(Ne[i] <= d && d < Ne[i+1]);

	// simple linear interpolation
	return (i + (d - Ne[i])/(Ne[i+1] - Ne[i]))*dr;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Use simple trapezoidal rule to integrate density
///
/// \param nelec number of electrons
/// \param dr delta r (uniform grid spacing)
/// \param ngrid number of uniform grid points
/// \param rho discretized density defined on uniform grid points i*dr, i = 0, ..., ngrid-1
/// \param nshpt number of (non-uniform) points per shell, defined by inverse function of density integral
/// \param data output structure
///
int IntegrateDensity(const int nelec, const double dr, const int ngrid, const double *rho, const int nshpt, dataNe_t *data)
{
	int i, j;
	int ret = 0;

	data->nelec = nelec;
	data->nshpt = nshpt;

	// TODO: use higher-level interpolation

	// data values of density integral
	double *Ne = malloc(ngrid*sizeof(double));
	Ne[0] = 0;
	for (i = 1; i < ngrid; i++)
	{
		// Jacobi element of polar coordinates: 2*pi*r, here r = i*dr
		Ne[i] = Ne[i-1] + M_PI*dr*dr*((i-1)*rho[i-1] + i*rho[i])/nelec;

		// consistency check
		if (Ne[i-1] >= Ne[i]) {
			fprintf(stderr, "Warning: integral of density 'Ne' not strictly increasing\n");
			ret |= E_DENSITY_INT;
		}
	}

	// rescale such that last point is exactly 1
	double scale = 1.0 / Ne[ngrid-1];
	for (i = 0; i < ngrid-1; i++)
	{
		Ne[i] *= scale;
	}
	Ne[ngrid-1] = 1;	// explicitly set to 1

	// create lookup table for inverse function
	int *invtab = malloc(ngrid*sizeof(int));
	j = 0;
	for (i = 0; i < ngrid; i++)
	{
		// look for last data point which is still smaller or equal i / ngrid
		double d = (double)i / ngrid;
		while (Ne[j+1] < d)
		{
			j++;
			assert(j < ngrid);
		}
		invtab[i] = j;
	}

	// calculate non-uniform grid defined by inverse function of density integral (includes shell radii at 'j*nshpt')
	data->rgrid = malloc((nelec*nshpt+1)*sizeof(double));
	for (i = 0; i < nelec*nshpt; i++)
	{
		data->rgrid[i] = EvaluateInvNe(dr, ngrid, Ne, invtab, (double)i/(nelec*nshpt));
	}
	// outermost point is at infinity
	data->rgrid[nelec*nshpt] = *(double *)&DBL_INF;

	// clean up
	free(invtab);
	free(Ne);

	return ret;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Free memory of density integral structure
///
void FreeDensityInt(dataNe_t *data)
{
	free(data->rgrid);
	data->nshpt = 0;
}
