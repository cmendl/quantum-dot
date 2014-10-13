/// \file interpolation.c
/// \brief Piecewise cubic polynomial interpolation.
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

#include "interpolation.h"
#include <math.h>
#include <stdlib.h>
#include <assert.h>


void AllocateInterpolation(const int n, interpolation_t *ip)
{
	ip->n = n;
	ip->x = (double *)malloc(n*sizeof(double));
	ip->f = (double *)malloc(n*sizeof(double));
}


void FreeInterpolation(interpolation_t *ip)
{
	free(ip->f);
	free(ip->x);
	ip->n = 0;
}


//_______________________________________________________________________________________________________________________
//


/// \brief Polynomial interpolation using the Lagrange form
static double EvaluateLagrancePoly(const double *x, const double *f, const int order, const double pt)
{
	int i, j;
	double ret = 0;

	for (j = 0; j <= order; j++)
	{
		double ell = 1;
		for (i = 0; i < j; i++)
		{
			ell *= (pt - x[i])/(x[j] - x[i]);
		}
		for (i = j+1; i <= order; i++)
		{
			ell *= (pt - x[i])/(x[j] - x[i]);
		}

		ret += f[j]*ell;
	}

	return ret;
}


double EvaluateInterpolationCubic(const interpolation_t *ip, const double pt)
{
	// check if 'pt' is within bounds
	// might have to use extrapolation

	if (pt < ip->x[1])
	{
		return EvaluateLagrancePoly(&ip->x[0], &ip->f[0], 3, pt);
	}

	if (pt >= ip->x[ip->n-2])
	{
		return EvaluateLagrancePoly(&ip->x[ip->n-4], &ip->f[ip->n-4], 3, pt);
	}

	// determine appropriate interval by binary search
	int imin = 1;	// start from second point
	int imax = ip->n-2;
	while (imin < imax - 1)
	{
		// calculate midpoint to cut set in half
		int imid = (imin + imax)/2;
 
		if (pt < ip->x[imid])
		{
			// pt is in lower half
			assert(imax > imid);
			imax = imid;
		}
		else
		{
			// pt is in upper half
			assert(imin < imid);
			imin = imid;
		}
	}
	assert(ip->x[imin] <= pt && pt < ip->x[imin+1]);

	return EvaluateLagrancePoly(&ip->x[imin-1], &ip->f[imin-1], 3, pt);
}
