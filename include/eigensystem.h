//	Solve the Kohn-Sham matrix equation (in the Laguerre basis of radial eigenfunctions of the 2D harmonic oscillator).
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

#ifndef EIGENSYSTEM_H
#define EIGENSYSTEM_H

#include "gauss_laguerre.h"
#include "potential.h"


/// \brief eigenvalues (energies) and eigenfunctions of the Kohn-Sham Hamiltonian
/// discretized on a Gauss-Laguerre grid
typedef struct
{
	int m;				//!< angular quantum number
	int npt;			//!< number of Gauss-Laguerre quadrature points
	int dim;			//!< dimension: number of basis functions
	double omega;		//!< omega parameter of harmonic potential
	double *en;			//!< energies (eigenvalues)
	double *psi;		//!< eigenfunctions (coefficients of Laguerre basis), matrix of size 'dim x dim'
	double *rho;		//!< corresponding density evaluated on Gauss-Laguerre grid, matrix of size 'npt x dim'
}
eigensystem_t;


void CalculateEigensystem(const potential_t *pot, const gauss_laguerre_t *gl, const double omega, eigensystem_t *eigensystem);

void FreeEigensystem(eigensystem_t *eigensystem);



#endif
