/// \file eigensystem.c
/// \brief Solve the Kohn-Sham matrix equation (in the Laguerre basis of radial eigenfunctions of the 2D harmonic oscillator).
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

#include "eigensystem.h"
#include "util.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>
#include <acml.h>



void CalculateEigensystem(const potential_t *pot, const gauss_laguerre_t *gl, const double omega, eigensystem_t *eigensystem)
{
	int i;
	int info;

	eigensystem->m = gl->m;
	eigensystem->npt = gl->npt;
	eigensystem->dim = gl->dim;
	eigensystem->omega = omega;

	// discretize potential
	eigensystem->psi = (double *)malloc(gl->dim*gl->dim*sizeof(double));
	DiscretizePotential(pot, gl, omega, eigensystem->psi);

	// add eigenvalues of Gauss-Laguerre eigenfunctions to diagonal
	for (i = 0; i < gl->dim; i++)
	{
		eigensystem->psi[i+i*eigensystem->dim] += 2*omega*(i+0.5 + 0.5*eigensystem->m);
	}

	eigensystem->en = (double *)malloc(gl->dim*sizeof(double));

	// find eigenvalues and correspoding eigenvectors of Hamiltonian;
	// 'dsyev' stores eigenvectors in 'eigensystem->psi'
	dsyev('V', 'U', gl->dim, eigensystem->psi, gl->dim, eigensystem->en, &info);
	if (info != 0)
	{
		fprintf(stderr, "dsyev failed, info: %i\n", info);
	}

	// calculate density
	eigensystem->rho = (double *)malloc(gl->npt*eigensystem->dim*sizeof(double));
	dgemm('T', 'N', gl->npt, eigensystem->dim, eigensystem->dim, 1.0, gl->u, gl->dim, eigensystem->psi, eigensystem->dim, 0, eigensystem->rho, eigensystem->npt);
	// take pointwise square
	for (i = 0; i < eigensystem->dim*gl->npt; i++)
	{
		eigensystem->rho[i] = omega/M_PI * square(eigensystem->rho[i]);
	}
}


void FreeEigensystem(eigensystem_t *eigensystem)
{
	free(eigensystem->rho);
	free(eigensystem->psi);
	free(eigensystem->en);
	eigensystem->dim = 0;
}
