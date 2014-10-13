/// \file potential.c
/// \brief Calculate SCE potential v, using Newton iteration to optimize the angles.
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

#include "potential.h"
#include "coulomb.h"
#include "comotion.h"
#include "angles.h"
#include "util.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include <acml.h>


//_______________________________________________________________________________________________________________________
///
/// \brief Remove entry at position 'i' from array
///
static double ArrayDrop(const int n, double *data, const int idrop)
{
	int i;
	double tmp = data[idrop];

	// shift-copy entries
	for (i = idrop; i < n-1; i++)
	{
		data[i] = data[i+1];
	}

	return tmp;
}

//_______________________________________________________________________________________________________________________
///
/// \brief Insert entry at position 'i', assuming that allocated memory is sufficient;
/// 'n' is the array length before insertion
///
static void ArrayInsert(const int n, double *data, const double val, const int iins)
{
	int i;

	// shift-copy entries
	for (i = n; i > iins; i--)
	{
		data[i] = data[i-1];
	}

	// copy value back
	data[iins] = val;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Generate potential structure based on non-uniform grid defined by inverse density function
///
error_desc_t GeneratePotential(const dataNe_t *data, const double *phiStart, potential_t *pot)
{
	int i, j, k;
	error_desc_t err = { 0 };

	pot->nelec = data->nelec;
	pot->nshpt = data->nshpt;

	// copy radial non-uniform grid
	pot->rgrid = malloc((pot->nelec*pot->nshpt+1) * sizeof(double));
	memcpy(pot->rgrid, data->rgrid, (pot->nelec*pot->nshpt+1) * sizeof(double));

	// allocate memory for 'dv' and 'v'
	// additional point for 'r = infinity'
	pot->dv  = malloc((pot->nelec*pot->nshpt+1) * sizeof(double));
	pot->v   = malloc((pot->nelec*pot->nshpt+1) * sizeof(double));
	pot->f   = malloc( pot->nelec*pot->nelec*pot->nshpt * sizeof(double));
	pot->phi = malloc( pot->nelec*pot->nelec*pot->nshpt * sizeof(double));
	pot->g   = malloc( pot->nelec*pot->nshpt * sizeof(double));

	// current 'f' and 'phi'
	double *f   = malloc(pot->nelec*sizeof(double));	// radial part of co-motion functions
	double *phi = malloc(pot->nelec*sizeof(double));	// polar angles

	if (pot->nelec == 2)
	{
		phi[0] = 0;
		phi[1] = M_PI;		// second electron always at opposite side

		// for each shell...
		for (i = 0; i < pot->nelec; i++)
		{
			// first case: electron radii coinciding with shell radii
			{
				EvaluateCoMotion(data, i, 0, f);

				if (i == 0)
				{
					assert(f[0] == 0 && isinf(f[1]));
					pot->dv[i*pot->nshpt] = 0;
				}
				else
				{
					// calculate derivative of potential via generalized Eq. (39) in Seidl, Gori-Giorgi and Savin, Phys. Rev. A 75, 042511 (2007)
					pot->dv[i*pot->nshpt] = VeeDr(pot->nelec, f, phi);
				}

				// save values
				memcpy(  &pot->f[i*pot->nshpt*pot->nelec],   f, pot->nelec*sizeof(double));
				memcpy(&pot->phi[i*pot->nshpt*pot->nelec], phi, pot->nelec*sizeof(double));
			}

			// next, electron radii different from shell radii
			for (j = 1; j < pot->nshpt; j++)
			{
				EvaluateCoMotion(data, i, j, f);

				// calculate derivative of potential via generalized Eq. (39) in Seidl, Gori-Giorgi and Savin, Phys. Rev. A 75, 042511 (2007)
				pot->dv[i*pot->nshpt + j] = VeeDr(pot->nelec, f, phi);

				// save values
				memcpy(&pot->f  [(i*pot->nshpt+j)*pot->nelec],   f, pot->nelec*sizeof(double));
				memcpy(&pot->phi[(i*pot->nshpt+j)*pot->nelec], phi, pot->nelec*sizeof(double));
			}
		}
	}
	else	// pot->nelec > 2
	{
		assert(pot->nelec > 2);

		error_desc_t cur_err;

		// temporary storage for angles
		double *chi	= (double *)malloc(pot->nelec*sizeof(double));

		// one round of angle optimization to refine starting values for actual optimization
		EvaluateCoMotion(data, 0, pot->nshpt/2, f);
		memcpy(phi, phiStart, pot->nelec*sizeof(double));
		assert(phiStart[0] == 0);

		OptimizeAngles(pot->nelec, f, phi);
		// skip error checking in preliminary angle optimization routine

		// optimize angles for first electron in first shell (for electron radii different from shell radii),
		// remaining co-motion function radii and angles result from group property of co-motion functions...
		// postpone angle optimization for electron radii coinciding with shell radii;
		// start from 'nshpt/2' outwards
		for (j = pot->nshpt/2; j < pot->nshpt; j++)
		{
			EvaluateCoMotion(data, 0, j, f);	// first electron in first shell

			cur_err = OptimizeAngles(pot->nelec, f, phi);
			AccumulateError(&err, &cur_err);

			//if (j == 1)		// first point inside shell
			//{
			//	// try disturbing angles a bit to find minimum value
			//	double optV = Vee(pot->nelec, f, phi);
			//	for (l = 0; l < 4; l++)
			//	{
			//		for (k = 1; k < pot->nelec; k++)	// angle of first electron always zero
			//		{
			//			double curV;

			//			memcpy(chi, phi, pot->nelec*sizeof(double));
			//			chi[k] += 0.01*M_PI;
			//			cur_err = OptimizeAngles(pot->nelec, f, chi);
			//			AccumulateError(&err, &cur_err);

			//			curV = Vee(pot->nelec, f, chi);
			//			if (curV < optV)
			//			{
			//				// save values
			//				memcpy(phi, chi, pot->nelec*sizeof(double));
			//				optV = curV;
			//			}
			//		}
			//	}
			//}

			// calculate derivative of potential via generalized Eq. (39) in Seidl, Gori-Giorgi and Savin, Phys. Rev. A 75, 042511 (2007)
			pot->dv[j] = VeeDr(pot->nelec, f, phi);

			// save values
			memcpy(&pot->f  [j*pot->nelec],   f, pot->nelec*sizeof(double));
			memcpy(&pot->phi[j*pot->nelec], phi, pot->nelec*sizeof(double));
		}

		// restore optimal angles from 'nshpt/2'
		memcpy(phi, &pot->phi[pot->nshpt/2 * pot->nelec], pot->nelec*sizeof(double));

		// proceed from 'nshpt/2 - 1' inwards
		for (j = pot->nshpt/2 - 1; j >= 1; j--)
		{
			EvaluateCoMotion(data, 0, j, f);	// first electron in first shell

			cur_err = OptimizeAngles(pot->nelec, f, phi);
			AccumulateError(&err, &cur_err);

			// calculate derivative of potential via generalized Eq. (39) in Seidl, Gori-Giorgi and Savin, Phys. Rev. A 75, 042511 (2007)
			pot->dv[j] = VeeDr(pot->nelec, f, phi);

			// save values
			memcpy(&pot->f  [j*pot->nelec],   f, pot->nelec*sizeof(double));
			memcpy(&pot->phi[j*pot->nelec], phi, pot->nelec*sizeof(double));
		}

		// remaining co-motion function radii and angles result from group property of co-motion functions
		for (i = 1; i < pot->nelec; i++)
		{
			// electron radii different from shell radii
			for (j = 1; j < pot->nshpt; j++)
			{
				// corresponding j index for first shell
				int j_base = (i % 2 == 0 ? j : pot->nshpt-j);

				for (k = 0; k < data->nelec; k++)
				{
					// shell index of current electron 'k'
					int k_base = GetCoMotionShell(pot->nelec, i, k);

					// for the first electron in first shell, each remaining electron 'k' is in shell 'k'
					f  [k] = pot->f  [pot->nelec*j_base + k_base];
					phi[k] = pot->phi[pot->nelec*j_base + k_base];
				}

				// shift angles such that first electron has angle zero
				double phi_shift = phi[0];
				for (k = 0; k < data->nelec; k++)
				{
					phi[k] -= phi_shift;
					if (phi[k] < -M_PI)		// using "if" instead of "while" to avoid infinite loop if (phi[k]) is +-1.#INF
					{
						phi[k] += 2*M_PI;
					}
					else if (phi[k] > M_PI)
					{
						phi[k] -= 2*M_PI;
					}
				}

				// calculate derivative of potential via generalized Eq. (39) in Seidl, Gori-Giorgi and Savin, Phys. Rev. A 75, 042511 (2007)
				pot->dv[i*pot->nshpt + j] = VeeDr(pot->nelec, f, phi);

				// save values
				memcpy(&pot->f  [(i*pot->nshpt+j)*pot->nelec],   f, pot->nelec*sizeof(double));
				memcpy(&pot->phi[(i*pot->nshpt+j)*pot->nelec], phi, pot->nelec*sizeof(double));
			}
		}

		// now optimize angles for electron radii coinciding with shell radii, using above optimized 'phi' values as starting point
		for (i = 0; i < pot->nelec; i++)
		{
			// for i==0 shell, first electron at origin, thus angles have rotational degeneracy and might jump when r > 0

			// index of electron at infinity, if any, otherwise -1
			int j_inf  = GetInfiniteCoMotionShell(pot->nelec, i);
			int j_zero = GetZeroCoMotionShell(i);
			// corresponding angles, saved for later
			double phi_inf = 0, phi_zero = 0;

			// effective number of electrons used in optimization
			int neff = pot->nelec;

			EvaluateCoMotion(data, i, 0, f);
			// check
			assert(j_inf  == -1 || isinf(f[j_inf]));
			assert(j_zero == -1 || f[j_zero] == 0);

			// starting point for optimization
			memcpy(phi, &pot->phi[(i*pot->nshpt+1)*pot->nelec], pot->nelec*sizeof(double));

			if (j_inf >= 0)
			{
				// remove electron at infinity for angular optimization
				phi_inf = ArrayDrop(neff, phi, j_inf);
				ArrayDrop(neff, f, j_inf);
				neff--;
				// account for shifted indices
				if (j_inf < j_zero) {
					j_zero--;
				}
			}
			if (j_zero >= 0)
			{
				// remove electron at origin for angular optimization
				phi_zero = ArrayDrop(neff, phi, j_zero);
				ArrayDrop(neff, f, j_zero);
				neff--;
			}

			// perform optimization
			cur_err = OptimizeAngles(neff, f, phi);
			AccumulateError(&err, &cur_err);

			// restore electron at origin, if any
			if (j_zero >= 0)
			{
				ArrayInsert(neff, phi, phi_zero, j_zero);
				ArrayInsert(neff, f,   0,        j_zero);
				neff++;

				assert(f[j_zero] == 0);
			}

			// calculate derivative of potential via generalized Eq. (39) in Seidl, Gori-Giorgi and Savin, Phys. Rev. A 75, 042511 (2007)
			pot->dv[i*pot->nshpt] = VeeDr(neff, f, phi);

			// restore 'phi' value of electron at infinity, if any
			if (j_inf >= 0)
			{
				// account for shifted indices
				if (j_inf <= j_zero) {
					j_zero++;
				}

				ArrayInsert(neff, phi, phi_inf,             j_inf);
				ArrayInsert(neff, f,   *(double *)&DBL_INF, j_inf);
				neff++;
			}

			// check
			assert(neff == pot->nelec);
			assert(j_inf  == -1 || isinf(f[j_inf]));
			assert(j_zero == -1 || f[j_zero] == 0);

			// save values
			memcpy(&pot->f  [i*pot->nshpt*pot->nelec], f,   pot->nelec*sizeof(double));
			memcpy(&pot->phi[i*pot->nshpt*pot->nelec], phi, pot->nelec*sizeof(double));
		}

		free(chi);
	}

	free(phi);
	free(f);


	// potential and its derivative are zero at infinity
	pot->dv[pot->nelec*pot->nshpt] = 0;
	pot->v [pot->nelec*pot->nshpt] = 0;

	// integrate derivative of potential to obtain actual potential;
	// integration direction is inwards to adhere to asymptotic (N-1)/r law

	// use asymptotic (N-1)/r law for outmost point
	pot->v[pot->nelec*pot->nshpt-1] = (pot->nelec - 1)/data->rgrid[pot->nelec*pot->nshpt-1];

	for (i = pot->nelec*pot->nshpt-2; i >= 0; i--)
	{
		// simple trapezoidal rule
		double dr = data->rgrid[i+1] - data->rgrid[i];
		pot->v[i] = pot->v[i+1] - 0.5*dr*(pot->dv[i] + pot->dv[i+1]);
	}


	// evaluate 'g' function

	for (i = 0; i < pot->nelec*pot->nshpt; i++)
	{
		pot->g[i] = Vee(pot->nelec, &pot->f[i*pot->nelec], &pot->phi[i*pot->nelec]);
		for (k = 0; k < pot->nelec; k++)
		{
			pot->g[i] -= pot->v[CoMotionIndexMap(pot->nelec, pot->nshpt, i, k)];
		}
	}


	return err;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Free memory of potential structure
///
void FreePotential(potential_t *pot)
{
	free(pot->g);
	free(pot->phi);
	free(pot->f);
	free(pot->v);
	free(pot->dv);
	free(pot->rgrid);
	pot->nshpt = 0;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Evaluate Kantorovich potential function v(r) at radius r
///
static inline double EvaluatePotential(const potential_t *pot, const double r)
{
	int i;

	assert(r >= 0);

	if (r >= pot->rgrid[pot->nelec*pot->nshpt-1])
	{
		// asymptotic (N-1)/r law
		return (pot->nelec - 1)/r;
	}

	// find appropriate grid point
	i = 0;
	int inext = pot->nelec*pot->nshpt-1;
	while (i+1 < inext)
	{
		// middle index
		int mid = (i + inext)/2;
		if (pot->rgrid[mid] <= r)
		{
			i = mid;
		}
		else
		{
			inext = mid;
		}
		assert(i < inext);
	}
	assert(pot->rgrid[i] <= r && r < pot->rgrid[i+1]);

	// simple linear interpolation
	double t = (r - pot->rgrid[i]) / (pot->rgrid[i+1] - pot->rgrid[i]);
	return pot->v[i] + t*(pot->v[i+1] - pot->v[i]);
}


//_______________________________________________________________________________________________________________________
///
/// \brief Discretize potential on Gauss-Laguerre grid
///
void DiscretizePotential(const potential_t *pot, const gauss_laguerre_t *gl, const double omega, double *vdiscr)
{
	int i, j;
	double *t   = (double *)malloc(gl->npt*sizeof(double));
	double *psi = (double *)malloc(gl->dim*gl->npt*sizeof(double));

	for (j = 0; j < gl->npt; j++)
	{
		double r = sqrt(gl->x[j]/omega);
		t[j] = EvaluatePotential(pot, r)*gl->w[j];
	}

	// effectively multiply by diagonal matrix with entries 't'
	for (j = 0; j < gl->npt; j++)
	{
		for (i = 0; i < gl->dim; i++)
		{
			psi[i+j*gl->dim] = gl->v[i+j*gl->dim]*t[j];
		}
	}

	// multiply matrices
	dgemm('N', 'T', gl->dim, gl->dim, gl->npt, 1.0, gl->v, gl->dim, psi, gl->dim, 0, vdiscr, gl->dim);

	free(psi);
	free(t);
}
