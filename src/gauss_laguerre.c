/// \file gauss_laguerre.c
/// \brief Generalized Gauss-Laguerre quadrature and evaluation of Laguerre polynomials.
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

#include "gauss_laguerre.h"
#include "util.h"
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>


#include "../data/gauss_laguerre_data.c"


void GenerateGaussLaguerre(const int m, const int npt, const int dim, gauss_laguerre_t *gl)
{
	int i, j;
	double sqmfac = sqrt((double)factorial(m));

	gl->m = m;
	gl->npt = npt;
	gl->dim = dim;

	// other cases not supported yet
	assert(npt == 64 || npt == 96);
	assert(dim == 16 || dim == 32);

	gl->x = (double *)malloc(npt*sizeof(double));
	gl->w = (double *)malloc(npt*sizeof(double));
	gl->v = (double *)malloc(npt*dim*sizeof(double));
	gl->u = (double *)malloc(npt*dim*sizeof(double));

	// TODO: universal code
	if (m == 0)
	{
		if (npt == 64)
		{
			memcpy(gl->x, x64m0, npt*sizeof(double));
			memcpy(gl->w, w64m0, npt*sizeof(double));
			if (dim == 16)
			{
				memcpy(gl->v, v64_16_m0, npt*dim*sizeof(double));
			}
			else if (dim == 32)
			{
				memcpy(gl->v, v64_32_m0, npt*dim*sizeof(double));
			}
			else
			{
				assert(false);
			}
		}
		else if (npt == 96)
		{
			memcpy(gl->x, x96m0, npt*sizeof(double));
			memcpy(gl->w, w96m0, npt*sizeof(double));
			if (dim == 16)
			{
				memcpy(gl->v, v96_16_m0, npt*dim*sizeof(double));
			}
			else if (dim == 32)
			{
				memcpy(gl->v, v96_32_m0, npt*dim*sizeof(double));
			}
			else
			{
				assert(false);
			}
		}
	}
	else if (m == 1)
	{
		if (npt == 64)
		{
			memcpy(gl->x, x64m1, npt*sizeof(double));
			memcpy(gl->w, w64m1, npt*sizeof(double));
			if (dim == 16)
			{
				memcpy(gl->v, v64_16_m1, npt*dim*sizeof(double));
			}
			else if (dim == 32)
			{
				memcpy(gl->v, v64_32_m1, npt*dim*sizeof(double));
			}
			else
			{
				assert(false);
			}
		}
		else if (npt == 96)
		{
			memcpy(gl->x, x96m1, npt*sizeof(double));
			memcpy(gl->w, w96m1, npt*sizeof(double));
			if (dim == 16)
			{
				memcpy(gl->v, v96_16_m1, npt*dim*sizeof(double));
			}
			else if (dim == 32)
			{
				memcpy(gl->v, v96_32_m1, npt*dim*sizeof(double));
			}
			else
			{
				assert(false);
			}
		}
		else
		{
			assert(false);
		}
	}
	else if (m == 2)
	{
		if (npt == 64)
		{
			memcpy(gl->x, x64m2, npt*sizeof(double));
			memcpy(gl->w, w64m2, npt*sizeof(double));
			if (dim == 16)
			{
				memcpy(gl->v, v64_16_m2, npt*dim*sizeof(double));
			}
			else if (dim == 32)
			{
				memcpy(gl->v, v64_32_m2, npt*dim*sizeof(double));
			}
			else
			{
				assert(false);
			}
		}
		else if (npt == 96)
		{
			memcpy(gl->x, x96m2, npt*sizeof(double));
			memcpy(gl->w, w96m2, npt*sizeof(double));
			if (dim == 16)
			{
				memcpy(gl->v, v96_16_m2, npt*dim*sizeof(double));
			}
			else if (dim == 32)
			{
				memcpy(gl->v, v96_32_m2, npt*dim*sizeof(double));
			}
			else
			{
				assert(false);
			}
		}
		else
		{
			assert(false);
		}
	}
	else if (m == 3)
	{
		if (npt == 64)
		{
			memcpy(gl->x, x64m3, npt*sizeof(double));
			memcpy(gl->w, w64m3, npt*sizeof(double));
			if (dim == 16)
			{
				memcpy(gl->v, v64_16_m3, npt*dim*sizeof(double));
			}
			else if (dim == 32)
			{
				memcpy(gl->v, v64_32_m3, npt*dim*sizeof(double));
			}
			else
			{
				assert(false);
			}
		}
		else if (npt == 96)
		{
			memcpy(gl->x, x96m3, npt*sizeof(double));
			memcpy(gl->w, w96m3, npt*sizeof(double));
			if (dim == 16)
			{
				memcpy(gl->v, v96_16_m3, npt*dim*sizeof(double));
			}
			else if (dim == 32)
			{
				memcpy(gl->v, v96_32_m3, npt*dim*sizeof(double));
			}
			else
			{
				assert(false);
			}
		}
		else
		{
			assert(false);
		}
	}
	else if (m == 4)
	{
		if (npt == 64)
		{
			memcpy(gl->x, x64m4, npt*sizeof(double));
			memcpy(gl->w, w64m4, npt*sizeof(double));
			if (dim == 16)
			{
				memcpy(gl->v, v64_16_m4, npt*dim*sizeof(double));
			}
			else if (dim == 32)
			{
				memcpy(gl->v, v64_32_m4, npt*dim*sizeof(double));
			}
			else
			{
				assert(false);
			}
		}
		else if (npt == 96)
		{
			memcpy(gl->x, x96m4, npt*sizeof(double));
			memcpy(gl->w, w96m4, npt*sizeof(double));
			if (dim == 16)
			{
				memcpy(gl->v, v96_16_m4, npt*dim*sizeof(double));
			}
			else if (dim == 32)
			{
				memcpy(gl->v, v96_32_m4, npt*dim*sizeof(double));
			}
			else
			{
				assert(false);
			}
		}
		else
		{
			assert(false);
		}
	}
	else if (m == 5)
	{
		if (npt == 64)
		{
			memcpy(gl->x, x64m5, npt*sizeof(double));
			memcpy(gl->w, w64m5, npt*sizeof(double));
			if (dim == 16)
			{
				memcpy(gl->v, v64_16_m5, npt*dim*sizeof(double));
			}
			else if (dim == 32)
			{
				memcpy(gl->v, v64_32_m5, npt*dim*sizeof(double));
			}
			else
			{
				assert(false);
			}
		}
		else if (npt == 96)
		{
			memcpy(gl->x, x96m5, npt*sizeof(double));
			memcpy(gl->w, w96m5, npt*sizeof(double));
			if (dim == 16)
			{
				memcpy(gl->v, v96_16_m5, npt*dim*sizeof(double));
			}
			else if (dim == 32)
			{
				memcpy(gl->v, v96_32_m5, npt*dim*sizeof(double));
			}
			else
			{
				assert(false);
			}
		}
		else
		{
			assert(false);
		}
	}
	else
	{
		assert(false);	// not supported yet
	}

	// rescale polynomials by sqrt(x^m / m!) exp(-x/2)
	for (j = 0; j < npt; j++)
	{
		double t = pow(gl->x[j], 0.5*m)/sqmfac * exp(-0.5*gl->x[j]);
		for (i = 0; i < dim; i++)
		{
			gl->u[i+j*dim] = t * gl->v[i+j*dim];
		}
	}
}


void FreeGaussLaguerre(gauss_laguerre_t *gl)
{
	free(gl->u);
	free(gl->v);
	free(gl->w);
	free(gl->x);
	gl->npt = 0;
}
