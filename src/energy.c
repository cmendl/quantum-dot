/// \file energy.c
/// \brief Sort and fill energy levels to obtain occupancy of each level.
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

#include "energy.h"
#include <stdlib.h>
#include <assert.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>



/// \brief Allocate memory for energy levels structure
void AllocateEnergyLevels(const int m_max, const int num, en_levels_t *levels)
{
	levels->m_max = m_max;
	levels->num   = num;

	levels->en = malloc((m_max+1)*num * sizeof(double));
}


/// \brief Free memory of energy levels structure
void FreeEnergyLevels(en_levels_t *levels)
{
	free(levels->en);

	levels->num   = 0;
	levels->m_max = 0;
}


//_______________________________________________________________________________________________________________________
//


/// \brief Comparison function used for sorting
static int en_occupancy_cmp(const void *e1, const void *e2)
{
	const en_occupancy_t *a = (en_occupancy_t *)e1;
	const en_occupancy_t *b = (en_occupancy_t *)e2;

	// sort according to energy
	if (a->en < b->en)
	{
		return -1;
	}
	else if (a->en > b->en)
	{
		return 1;
	}

	// should actually never reach this point

	// compare 'm' quantum numbers
	if (a->m < b->m)
	{
		return -1;
	}
	else if (a->m > b->m)
	{
		return 1;
	}

	// compare energy quantum numbers
	if (a->n < b->n)
	{
		return -1;
	}
	else if (a->n > b->n)
	{
		return 1;
	}

	return 0;
}


/// \brief Get the maximum occupancy of an energy level, depending on 'm' quantum number and whether state is spin-polarized
static inline int GetMaxOccupancy(const int m, const bool spinpol)
{
	int occ = 1;

	if (m != 0)
	{
		occ *= 2;
	}
	if (!spinpol)
	{
		occ *= 2;
	}

	return occ;
}


//_______________________________________________________________________________________________________________________
//


/// \brief Fermi-Dirac distribution function
static double FermiDirac(const double kBT, const double mu, const double en)
{
	return 1 / (exp((en - mu)/kBT) + 1);
}


/// \brief GSL parameters for NumParticleDifference()
typedef struct
{
	int nelec;			//!< physical number of electrons
	int nlevels;		//!< number of energy levels
	double *en;			//!< energies
	int *occ;			//!< maximum occupancy of each energy level
	double kBT;			//!< k_B T (Boltzmann constant times temperature)
}
gsl_FD_params_t;


/// \brief Calculate difference in number of particles given a Fermi-Dirac distribution with chemical potential 'mu'
static double NumParticleDifference(double mu, void *p)
{
	int i;

	const gsl_FD_params_t *params = (gsl_FD_params_t *)p;

	double N = 0;
	for (i = 0; i < params->nlevels; i++)
	{
		N += params->occ[i] * FermiDirac(params->kBT, mu, params->en[i]);
	}

	return N - params->nelec;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Calculate chemical potential 'mu' such that FermiDirac distribution gives exactly 'nelec' electrons
///
static double CalculateChemicalPotential(const en_levels_t *levels, const int nelec, const bool spinpol, const double kBT)
{
	int m, n, i;

	const int maxiter = 128;

	gsl_FD_params_t params;
	params.nelec = nelec;
	params.nlevels = (levels->m_max+1)*levels->num;
	params.en  = malloc(params.nlevels * sizeof(double));
	params.occ = malloc(params.nlevels * sizeof(int));
	params.kBT = kBT;
	// fill energy levels and occupancies
	for (m = 0; m <= levels->m_max; m++)
	{
		for (n = 0; n < levels->num; n++)
		{
			params.en [n + m*levels->num] = levels->en[n + m*levels->num];
			params.occ[n + m*levels->num] = GetMaxOccupancy(m, spinpol);
		}
	}

	// first interval used in interval bisection method
	double mu_lo = 0.0;
	double mu_hi = levels->en[(levels->m_max+1)*levels->num-1];		// use "last" energy level as upper bound

	gsl_function F;
	F.function = &NumParticleDifference;
	F.params = &params;

	const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
	gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s, &F, mu_lo, mu_hi);

	double mu = 0;
	int status = GSL_CONTINUE;
	for (i = 0; i < maxiter && status == GSL_CONTINUE; i++)
	{
		status = gsl_root_fsolver_iterate(s);
		if (status != GSL_SUCCESS) {
			fprintf(stderr, "CalculateOccupancy() warning: gsl_root_fsolver_iterate() returned with status code %i\n", status);
		}

		mu = gsl_root_fsolver_root(s);
		mu_lo  = gsl_root_fsolver_x_lower(s);
		mu_hi  = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval(mu_lo, mu_hi, 0, 1e-8);
	}

	if (status != GSL_SUCCESS) {
		fprintf(stderr, "CalculateOccupancy() warning: not converged after %i iterations, status: %i\n", maxiter, status);
	}

	// clean up
	free(params.occ);
	free(params.en);
	gsl_root_fsolver_free(s);

	return mu;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Calculate occupancy of energy levels
///
void CalculateOccupancy(const en_levels_t *levels, const int nelec, const bool spinpol, const double kBT, en_occlist_t *occlist)
{
	int i, m, n;

	// accumulate all energies into one list
	occlist->occupancy = malloc((levels->m_max+1)*levels->num * sizeof(en_occupancy_t));

	if (kBT == 0)	// zero temperature
	{
		for (m = 0; m <= levels->m_max; m++)
		{
			for (n = 0; n < levels->num; n++)
			{
				en_occupancy_t *encur = &occlist->occupancy[n + m*levels->num];
				encur->en  = levels->en[n + m*levels->num];
				encur->m   = m;
				encur->n   = n;
				encur->occ = 0;		// initially set to zero
			}
		}

		// sort energies
		qsort(occlist->occupancy, (levels->m_max+1)*levels->num, sizeof(en_occupancy_t), en_occupancy_cmp);

		// fill up energy levels
		en_occupancy_t *encur = occlist->occupancy;
		occlist->length = 1;
		for (i = 0; i < nelec; i++)
		{
			if (encur->occ >= GetMaxOccupancy(encur->m, spinpol)) {
				encur++;	// use next level
				occlist->length++;
			}

			encur->occ++;
		}
	}
	else	// k_B T > 0
	{
		assert(kBT > 0);

		const double mu = CalculateChemicalPotential(levels, nelec, spinpol, kBT);

		// all energy levels have nonzero occupancy since Fermi-Dirac function is strictly positive
		occlist->length = (levels->m_max+1)*levels->num;
		for (m = 0; m <= levels->m_max; m++)
		{
			for (n = 0; n < levels->num; n++)
			{
				en_occupancy_t *encur = &occlist->occupancy[n + m*levels->num];
				encur->en  = levels->en[n + m*levels->num];
				encur->m   = m;
				encur->n   = n;
				encur->occ = GetMaxOccupancy(m, spinpol) * FermiDirac(kBT, mu, encur->en);
			}
		}

		// sort energies
		qsort(occlist->occupancy, (levels->m_max+1)*levels->num, sizeof(en_occupancy_t), en_occupancy_cmp);
	}
}



/// \brief Free memory of energy level occupancy list
void FreeOccupancy(en_occlist_t *occlist)
{
	free(occlist->occupancy);
	occlist->length = 0;
}
