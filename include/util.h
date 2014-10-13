//	Utility functions.
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

#ifndef UTIL_H
#define UTIL_H

#include <stdbool.h>


#ifdef _MSC_VER

#define strcasecmp _strcmpi

#endif


//_______________________________________________________________________________________________________________________
//

/// \brief Binary bit pattern representing INF (double format)
static const unsigned long long DBL_INF = 0x7ff0000000000000;


//_______________________________________________________________________________________________________________________
//

static inline double square(double x)
{
	return x*x;
}


static inline int factorial(int n)
{
	int f;

	if (n < 0) {
		return -1;
	}

	// handle the case n == 0
	f = 1;
	while (n > 1)
	{
		f *= n;
		--n;
	}

	return f;
}


//_______________________________________________________________________________________________________________________
//

static inline int mini(int x, int y)
{
	if (x <= y)
	{
		return x;
	}
	else
	{
		return y;
	}
}

static inline int maxi(int x, int y)
{
	if (x >= y)
	{
		return x;
	}
	else
	{
		return y;
	}
}


//_______________________________________________________________________________________________________________________
//

static inline double minf(double x, double y)
{
	if (x <= y)
	{
		return x;
	}
	else
	{
		return y;
	}
}

static inline double maxf(double x, double y)
{
	if (x >= y)
	{
		return x;
	}
	else
	{
		return y;
	}
}


//_______________________________________________________________________________________________________________________
//

double Norm(const double *x, int n);


//_______________________________________________________________________________________________________________________
//

int ReadDoubleData(const char *filename, double *data, int n);

int WriteDoubleData(const char *filename, double *data, int n, bool append);



#endif
