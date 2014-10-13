//	Generalized Gauss-Laguerre quadrature and evaluation of Laguerre polynomials.
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

#ifndef GAUSS_LAGUERRE_H
#define GAUSS_LAGUERRE_H


/// \brief Gauss-Laguerre quadrature and Laguerre polynomials evaluated at quadrature points
typedef struct
{
	int m;			//!< angular momentum quantum number
	int npt;		//!< number of points
	int dim;		//!< number of basis functions used
	double *x;		//!< points
	double *w;		//!< weights
	double *v;		//!< first 'dim' Laguerre polynomials evaluated at points 'x', matrix of dimension 'dim x n'
	double *u;		//!< first 'dim' Laguerre polynomials multiplied by sqrt of weight function, sqrt(x^m / m!) exp(-x/2)
}
gauss_laguerre_t;


void GenerateGaussLaguerre(const int m, const int npt, const int dim, gauss_laguerre_t *gl);

void FreeGaussLaguerre(gauss_laguerre_t *gl);



#endif
