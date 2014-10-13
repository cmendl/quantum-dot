//	Calculate SCE potential v, using Newton iteration to optimize the angles.
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

#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "density_int.h"
#include "gauss_laguerre.h"
#include "error_desc.h"


//_______________________________________________________________________________________________________________________
///
/// \brief SCE potential
///
/// radial grid points Ne^{-1}(i + j/nshpt), j = 0, 1, ..., nshpt-1 for each shell
///
typedef struct
{
	int nelec;			//!< number of electrons
	int nshpt;			//!< number of (non-uniform) radial grid points for each shell
	double *rgrid;		//!< radial grid points Ne^{-1}(i + j/nshpt), j = 0, 1, ..., nshpt-1 for each shell i
	double *dv;			//!< derivative of SCE potential evaluated at grid points
	double *v;			//!< SCE potential (integral of derivative) evaluated at grid points
	double *f;			//!< radial co-motion functions evaluated at grid points
	double *phi;		//!< corresponding optimized angles
	double *g;			//!< 'g' function evaluated at grid points
}
potential_t;


error_desc_t GeneratePotential(const dataNe_t *data, const double *phiStart, potential_t *pot);

void FreePotential(potential_t *pot);


//_______________________________________________________________________________________________________________________
//


// assuming that 'vdiscr' points to an array of length gl->dim^2
void DiscretizePotential(const potential_t *pot, const gauss_laguerre_t *gl, const double omega, double *vdiscr);



#endif
