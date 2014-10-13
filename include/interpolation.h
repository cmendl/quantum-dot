//	Piecewise cubic polynomial interpolation.
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

#ifndef INTERPOLATION_H
#define INTERPOLATION_H


/// \brief Function evaluated at discrete points,
/// for constructing a piecewise cubic polynomial interpolation
typedef struct
{
	int n;			//!< number of points
	double *x;		//!< points (sorted in ascending order)
	double *f;		//!< corresponding function values at points
}
interpolation_t;


void AllocateInterpolation(const int n, interpolation_t *ip);

void FreeInterpolation(interpolation_t *ip);


//_______________________________________________________________________________________________________________________
//


// evaluate at point 'pt'
double EvaluateInterpolationCubic(const interpolation_t *ip, const double pt);



#endif
