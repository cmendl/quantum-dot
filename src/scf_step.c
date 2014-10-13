/// \file scf_step.c
/// \brief One self-consistent field (SCF) iteration step.
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

#include "scf_step.h"
#include "eigensystem.h"
#include "interpolation.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <acml.h>



//_______________________________________________________________________________________________________________________
///
/// \brief Calculate a self-consistent field (SCF) iteration step
///
/// \param nelec number of electrons
/// \param omega spring constant of harmonic oscillator
/// \param rho radial density discretized on uniform grid
/// \param dr delta r of uniform radial grid for density
/// \param ngrid number of radial grid points for density
/// \param nshpt number of (uniform) radial grid points within each shell, used to discretize co-motion functions
/// \param spinpol spin-polarized or not
/// \param kBT k_B T (Boltzmann constant times temperature)
/// \param phiStart initial starting angles for optimization
/// \param pot data structure for Kantorovich potential (output)
/// \param occlist energy level occupancy (output)
/// \param rho_next new density, discretized on uniform radial grid (output)
///
/// \return detailed error description, if any
///
error_desc_t SelfConsistentFieldStep(const int nelec, const double omega, const double *rho, const double dr, const int ngrid, const int nshpt,
	const bool spinpol, const double kBT, const double *phiStart, potential_t *pot, en_occlist_t *occlist, double *rho_next)
{
	const int m_max = (nelec <= 10 ? 3 : 5);	// maximum angular quantum number
	const int npt = 64;		// number of Gauss-Laguerre points
	const int dim = 16;		// number of basis functions

	int i, m;
	error_desc_t err = { 0 };

	dataNe_t data;
	err.codes |= IntegrateDensity(nelec, dr, ngrid, rho, nshpt, &data);

	error_desc_t cur_err = GeneratePotential(&data, phiStart, pot);
	AccumulateError(&err, &cur_err);

	// at most 'dim' levels
	en_levels_t levels;
	AllocateEnergyLevels(m_max, mini(dim, spinpol ? nelec/2 : (nelec+1)/4), &levels);

	gauss_laguerre_t *gl  = malloc((m_max+1)*sizeof(gauss_laguerre_t));
	eigensystem_t *eigsys = malloc((m_max+1)*sizeof(eigensystem_t));
	// potential is evaluated on Gauss-Laguerre quadrature points
	for (m = 0; m <= m_max; m++)
	{
		// Gauss-Laguerre quadrature
		GenerateGaussLaguerre(m, npt, dim, &gl[m]);

		CalculateEigensystem(pot, &gl[m], omega, &eigsys[m]);

		// copy eigenvalues
		memcpy(&levels.en[m*levels.num], eigsys[m].en, levels.num*sizeof(double));
	}

	CalculateOccupancy(&levels, nelec, spinpol, kBT, occlist);

	// new density for each 'm', evaluated on Gauss-Laguerre grid
	interpolation_t *ip_rho1 = malloc((m_max+1)*sizeof(interpolation_t));
	// accumulate density for each 'm' separately (evaluated on Gauss-Laguerre points)
	for (m = 0; m <= m_max; m++)
	{
		AllocateInterpolation(npt, &ip_rho1[m]);

		// r = sqrt(x/omega)
		for (i = 0; i < npt; i++)
		{
			ip_rho1[m].x[i] = sqrt(gl[m].x[i]/omega);
		}

		// set density values initially to zero (for accumulating density)
		memset(ip_rho1[m].f, 0, npt*sizeof(double));
	}
	for (i = 0; i < occlist->length; i++)
	{
		m = occlist->occupancy[i].m;
		// rho1[m] += occ * rho_level
		daxpy(npt, occlist->occupancy[i].occ, &eigsys[m].rho[occlist->occupancy[i].n*npt], 1, ip_rho1[m].f, 1);
	}

	// use interpolation to evaluate density on original 'rho' grid
	memset(rho_next, 0, ngrid*sizeof(double));
	for (i = 0; i < ngrid; i++)
	{
		for (m = 0; m <= m_max; m++)
		{
			rho_next[i] += EvaluateInterpolationCubic(&ip_rho1[m], i*dr);
		}
	}

	// clean up
	FreeEnergyLevels(&levels);
	for (m = 0; m <= m_max; m++)
	{
		FreeInterpolation(&ip_rho1[m]);
		FreeEigensystem(&eigsys[m]);
		FreeGaussLaguerre(&gl[m]);
	}
	free(ip_rho1);
	free(eigsys);
	free(gl);
	FreeDensityInt(&data);

	return err;
}
