/// \file comotion.c
/// \brief Evaluate the radial co-motion functions f_i.
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

#include "comotion.h"
#include "util.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>



//_______________________________________________________________________________________________________________________
///
/// \brief Retrieve the grid index of electron 'k', given the grid index 'ir' of the first electron
///
int CoMotionIndexMap(const int nelec, const int nshpt, const int ir, const int k)
{
	assert(0 <= ir && ir <= nelec*nshpt);

	if (k % 2 == 1)	// k odd
	{
		return abs((k+1)*nshpt - ir);
	}
	else			// k even
	{
		return nelec*nshpt - abs(nelec*nshpt - k*nshpt - ir);
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Evaluate co-motion functions f at the radius which is inside shell 'shellindex' at the non-uniform grid point 'pointindex'
///
void EvaluateCoMotion(const dataNe_t *data, const int shellindex, const int pointindex, double *f)
{
	int k;

	assert(0 <= pointindex && pointindex < data->nshpt);
	assert(shellindex < data->nelec || (shellindex == data->nelec && pointindex == 0));		// allow outmost point at infinity

	// effective index of current r
	const int ir = pointindex + shellindex*data->nshpt;

	// for pointindex == 0, effectively evaluate co-motion functions at r = a_shellindex, and one entry in 'f' can become 'inf'

	for (k = 0; k < data->nelec; k++)
	{
		f[k] = data->rgrid[CoMotionIndexMap(data->nelec, data->nshpt, ir, k)];
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Return the shell of electron 'k' for first electron within shell 'shellindex'
///
int GetCoMotionShell(const int nelec, const int shellindex, const int k)
{
	assert(0 <= k);
	assert(0 <= shellindex && shellindex < nelec);

	if (k % 2 == 1)	// k odd
	{
		return mini(abs(k + 1 - shellindex), abs(k - shellindex));
	}
	else			// k even
	{
		return mini(nelec - abs(nelec - k - shellindex), nelec - abs(nelec - k - shellindex - 1));
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Return the index of the electron at infinity, or -1 if all electron radii are finite
///
int GetInfiniteCoMotionShell(const int nelec, const int shellindex)
{
	int i;

	assert(0 <= shellindex && shellindex <= nelec);

	// special case
	if (nelec % 2 == 0 && shellindex == 0)
	{
		return nelec-1;
	}

	i = nelec - shellindex;
	if (i % 2 == 0)
	{
		return i;
	}
	else
	{
		return -1;
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Return the index of the electron at the origin, or -1 if all electron radii are larger than zero
///
int GetZeroCoMotionShell(const int shellindex)
{
	assert(0 <= shellindex);

	// special case
	if (shellindex == 0) {
		return 0;
	}

	if (shellindex % 2 == 0)
	{
		return shellindex - 1;
	}
	else
	{
		return -1;
	}
}
