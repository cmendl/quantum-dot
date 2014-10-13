/// \file angles.c
/// \brief Optimize polar angles with respect to Coulomb repulsion, given the radii.
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

#include "angles.h"
#include "coulomb.h"
#include "util.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <float.h>
#include <assert.h>
#include <stdio.h>
// LAPACK functions
#include <acml.h>
// GSL library
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_siman.h>


static const int maxiter		= 64;			//!< maximum number of iterations for local angle optimization
static const double step_size	= 0.01;			//!< step size for angle optimization


//_______________________________________________________________________________________________________________________
//
// GSL wrappers of Vee() for local minimization

/// \brief GSL parameters for angle optimization
typedef struct
{
	double phi0;
	const double *r;
}
gsl_angle_params_t;


/// \brief GSL wrapper for Vee()
static double gslVee(const gsl_vector *phiRest, void *p)
{
	int i;
	const gsl_angle_params_t *params = (gsl_angle_params_t *)p;
	double v_val;

	// first angle is held fixed
	const int nelec = (int)phiRest->size + 1;

	double *phi = (double *)malloc(nelec*sizeof(double));
	phi[0] = params->phi0;
	for (i = 1; i < nelec; i++) {
		phi[i] = gsl_vector_get(phiRest, i-1);
	}

	v_val = Vee(nelec, params->r, phi);

	free(phi);

	return v_val;
}

/// \brief GSL wrapper for VeeGradPhi()
static void gslVeeGradPhi(const gsl_vector *phiRest, void *p, gsl_vector *df)
{
	int i;
	const gsl_angle_params_t *params = (gsl_angle_params_t *)p;
	double *grad;

	// first angle is held fixed
	const int nelec = (int)phiRest->size + 1;

	double *phi = (double *)malloc(nelec*sizeof(double));
	phi[0] = params->phi0;
	for (i = 1; i < nelec; i++) {
		phi[i] = gsl_vector_get(phiRest, i-1);
	}

	grad = (double *)malloc(nelec*sizeof(double));

	// calculate gradient
	VeeGradPhi(nelec, params->r, phi, grad);

	// ignore derivative with respect to first angle (held fixed)
	for (i = 1; i < nelec; i++) {
		gsl_vector_set(df, i-1, grad[i]);
	}

	free(grad);
	free(phi);
}

/// \brief GSL wrapper for both 'Vee()' and 'VeeGradPhi()' together
static void gslVeeValGradPhi(const gsl_vector *phiRest, void *p, double *f, gsl_vector *df)
{
	int i;
	const gsl_angle_params_t *params = (gsl_angle_params_t *)p;
	double *phi;
	double *grad;

	// first angle is held fixed
	const int nelec = (int)phiRest->size + 1;

	phi = (double *)malloc(nelec*sizeof(double));
	phi[0] = params->phi0;
	for (i = 1; i < nelec; i++) {
		phi[i] = gsl_vector_get(phiRest, i-1);
	}

	*f = Vee(nelec, params->r, phi);

	grad = (double *)malloc(nelec*sizeof(double));

	// calculate gradient
	VeeGradPhi(nelec, params->r, phi, grad);

	// ignore derivative with respect to first angle (which is held fixed)
	for (i = 1; i < nelec; i++) {
		gsl_vector_set(df, i-1, grad[i]);
	}

	free(grad);
	free(phi);
}


//_______________________________________________________________________________________________________________________
///
/// \brief Optimize polar angles with respect to Coulomb repulsion, given the radii
///
/// \param nelec number of electrons
/// \param r radii, array of length 'nelec'
/// \param phi serves as input (starting values) and output (optimized angles); array of length 'nelec';
///        phi[0] is not changed due to invariance under global rotation
///
/// \return detailed error description, if any
///
error_desc_t OptimizeAngles(const int nelec, const double *r, double *phi)
{
	int i, j, k;
	error_desc_t err = { 0 };
	int status;

	// tolerance for convergence test based on gradient
	const double tol = nelec < 7 ? 1e-12 : 1e-10;

	const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_vector_bfgs2;	// alternative: gsl_multimin_fdfminimizer_conjugate_fr
	gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc(T, nelec-1);

	gsl_angle_params_t params;
	params.phi0 = phi[0];
	params.r = r;

	gsl_multimin_function_fdf my_func;
	my_func.n = nelec-1;
	my_func.f = gslVee;
	my_func.df = gslVeeGradPhi;
	my_func.fdf = gslVeeValGradPhi;
	my_func.params = &params;

	// starting point
	// omit first angle phi[0]
	gsl_vector *x = gsl_vector_alloc(nelec-1);
	for (i = 1; i < nelec; i++) {
		gsl_vector_set(x, i-1, phi[i]);
	}

	gsl_multimin_fdfminimizer_set(s, &my_func, x, step_size, 0.1);

	for (k = 0; k < maxiter; k++)
	{
		status = gsl_multimin_fdfminimizer_iterate(s);

		// interval [-pi, pi]
		for (i = 0; i < nelec-1; i++)
		{
			double curphi = gsl_vector_get(s->x, i);
			assert(!isinf(curphi));

			if (curphi < -M_PI)		// using "if" instead of "while" to avoid infinite loop if curphi is +-1.#INF
			{
				curphi += 2*M_PI;
			}
			else if (curphi > M_PI)
			{
				curphi -= 2*M_PI;
			}

			gsl_vector_set(s->x, i, curphi);
		}

		// check status
		if (status != 0) {
			err.codes |= E_ITERATION;
			err.msg = gsl_strerror(status);
			fprintf(stderr, "Warning: 'gsl_multimin_fdfminimizer_iterate' failed at iteration %i, error: %s.\n", k, err.msg);
			break;
		}

		// convergence test
		status = gsl_multimin_test_gradient(s->gradient, tol);
		if (status == GSL_SUCCESS) {
			printf("'gsl_multimin_test_gradient()' successful.\n");
			break;
		}
	}

	// norm of gradient
	double delta = 0;
	for (i = 0; i < nelec-1; i++)
	{
		delta += square(gsl_vector_get(s->gradient, i));
	}
	// 2-norm
	delta = sqrt(delta);
	printf("final delta: %g\n", delta);

	if (k == maxiter)
	{
		fprintf(stderr, "Warning: 'OptimizeAngles' did not converge within %i iterations, final delta: %g.\n", maxiter, delta);
		err.codes |= E_CONVERGENCE;
		err.conv_delta = delta;
	}

	// copy values back
	for (i = 1; i < nelec; i++) {
		phi[i] = gsl_vector_get(s->x, i-1);
	}

	gsl_multimin_fdfminimizer_free(s);
	gsl_vector_free(x);


	// consistency check: all eigenvalues of Hessian matrix should be non-negative
	{
		double *hessfull = (double *)malloc(nelec*nelec*sizeof(double));
		double *hesssubm = (double *)malloc((nelec-1)*(nelec-1)*sizeof(double));
		double *eigenval = (double *)malloc((nelec-1)*sizeof(double));

		// calculate Hessian matrix
		VeeHessPhi(nelec, r, phi, hessfull);
		// omit first row and column in Hessian matrix since phi[0] is fixed
		for (j = 1; j < nelec; j++)
		{
			for (i = 1; i < nelec; i++)
			{
				hesssubm[(i-1)+(j-1)*(nelec-1)] = hessfull[i+j*nelec];
			}
		}
		// calculate eigenvalues
		dsyev('N', 'U', nelec-1, hesssubm, nelec-1, eigenval, &status);
		// check for errors
		if (status != 0) {
			fprintf(stderr, "Warning: 'dsyev' failed, info: %i\n", status);
			err.codes |= E_LAPACK;
		}
		if (eigenval[0] < 0) {
			fprintf(stderr, "Warning: Hessian has negative eigenvalue %g.\n", eigenval[0]);
			err.codes |= E_NEG_HESS;
			err.hess_eig = eigenval[0];
		}

		free(eigenval);
		free(hesssubm);
		free(hessfull);
	}

	return err;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Perform a global optimization of the polar angles with respect to Coulomb repulsion, given the radii
///
/// \param nelec number of electrons
/// \param r radii, array of length 'nelec'
/// \param phi serves as input (starting values) and output (optimized angles); array of length 'nelec';
///        phi[0] is not changed due to invariance under global rotation
/// \param optVee optimized Coulomb repulsion (output)
/// \return detailed error description, if any
///
error_desc_t OptimizeAnglesGlobal(const int nelec, const double *r, double *phi, double *optVee)
{
	int i, k, l;

	error_desc_t err = OptimizeAngles(nelec, r, phi);

	// temporary storage for angles
	double *chi	= (double *)malloc(nelec*sizeof(double));

	// try permuting angles to find minimum value
	(*optVee) = Vee(nelec, r, phi);
	for (i = 0; i < 32; i++)		// iterate several times
	{
		// for all pairs of electrons...
		for (l = 1; l < nelec; l++)		// angle of first electron always zero
		{
			for (k = 1; k < l; k++)
			{
				double curV;

				memcpy(chi, phi, nelec*sizeof(double));
				// permute k <-> l
				chi[k] = phi[l];
				chi[l] = phi[k];
				error_desc_t cur_err = OptimizeAngles(nelec, r, chi);

				curV = Vee(nelec, r, chi);
				if (curV < (*optVee))
				{
					// only accumulate error if we really use the optimized angles
					AccumulateError(&err, &cur_err);

					// save values
					memcpy(phi, chi, nelec*sizeof(double));
					(*optVee) = curV;
				}
			}
		}
	}

	// clean up
	free(chi);

	return err;
}


/*
//_______________________________________________________________________________________________________________________
//
// GSL wrappers of Vee() for global minimization by simulated annealing


typedef struct
{
	int nelec;
	const double *r;
	double phi[10];
}
saConfig_t;


static double saVee(void *xp)
{
	const saConfig_t *x = (saConfig_t *)xp;

	return Vee(x->nelec, x->r, x->phi);
}


static void saStep(const gsl_rng *randgen, void *xp, double step_size)
{
	int i;
	saConfig_t *x = (saConfig_t *)xp;

	// first angle always fixed
	for (i = 1; i < x->nelec; i++)
	{
		double u = gsl_rng_uniform(randgen);
		x->phi[i] += step_size * (u - 0.5);		// divide by nelec ?

		// [-pi,pi] range
		if (x->phi[i] < -M_PI)
		{
			x->phi[i] += 2*M_PI;
		}
		else if (x->phi[i] > M_PI)
		{
			x->phi[i] -= 2*M_PI;
		}
	}
}


static double saDistance(void *xp, void *yp)
{
	int i;

	const saConfig_t *x = (saConfig_t *)xp;
	const saConfig_t *y = (saConfig_t *)yp;

	// distance between angles
	double dphi[10];
	for (i = 0; i < x->nelec; i++)
	{
		dphi[i] = x->phi[i] - y->phi[i];

		// [-pi,pi] range
		if (dphi[i] < -M_PI)
		{
			dphi[i] += 2*M_PI;
		}
		else if (dphi[i] > M_PI)
		{
			dphi[i] -= 2*M_PI;
		}
	}

	assert(x->nelec == y->nelec);

	return Norm(dphi + 1, x->nelec - 1);	// skip first angle (always held fixed)
}


static void saPrint(void *xp)
{
	int i;
	const saConfig_t *x = (saConfig_t *)xp;
	
	printf("   angles phi: ");
	for (i = 0; i < x->nelec; i++) {
		printf("%f ", x->phi[i]);
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Globally optimize polar angles with respect to Coulomb repulsion, given the radii
///
/// \param nelec number of electrons
/// \param r radii, array of length 'nelec'
/// \param phi serves as input (starting values) and output (optimized angles); array of length 'nelec';
///        phi[0] is not changed due to invariance under global rotation
///
/// \return detailed error description, if any
///
error_desc_t OptimizeAnglesGlobal(const int nelec, const double *r, double *phi)
{
	const gsl_rng_type *T;
	gsl_rng *randgen;

	gsl_siman_params_t params;

	error_desc_t err = { 0 };

	// initial configuration
	saConfig_t config;
	config.nelec = nelec;
	config.r = r;
	memcpy(config.phi, phi, nelec*sizeof(double));

	gsl_rng_env_setup();

	// random number generator
	T = gsl_rng_default;
	randgen = gsl_rng_alloc(T);

	// parameters
	params.n_tries			= 200;				// how many points do we try before stepping
	params.iters_fixed_T	= 1000;				// how many iterations for each T?
	params.step_size		= 0.01 * M_PI;		// max step size in random walk
	params.k				= 1.0;				// virtual Boltzmann constant
	params.t_initial		= 0.001;			// initial temperature
	params.mu_t				= 1.004;			// damping factor for temperature
	params.t_min			= 2.0e-7;			// minimum temperature

	gsl_siman_solve(randgen, &config, saVee, saStep, saDistance, saPrint, NULL, NULL, NULL, sizeof(saConfig_t), params);

	// first angle must not change
	assert(config.phi[0] == phi[0]);

	// store optimized angles in 'phi'
	memcpy(phi, config.phi, nelec*sizeof(double));

	// clean up
	gsl_rng_free(randgen);


	// consistency check: all eigenvalues of Hessian matrix should be non-negative
	{
		int i, j;
		int status;
		double *hessfull = (double *)malloc(nelec*nelec*sizeof(double));
		double *hesssubm = (double *)malloc((nelec-1)*(nelec-1)*sizeof(double));
		double *eigenval = (double *)malloc((nelec-1)*sizeof(double));

		// calculate Hessian matrix
		VeeHessPhi(nelec, r, phi, hessfull);
		// omit first row and column in Hessian matrix since phi[0] is fixed
		for (j = 1; j < nelec; j++)
		{
			for (i = 1; i < nelec; i++)
			{
				hesssubm[(i-1)+(j-1)*(nelec-1)] = hessfull[i+j*nelec];
			}
		}
		// calculate eigenvalues
		dsyev('N', 'U', nelec-1, hesssubm, nelec-1, eigenval, &status);
		// check for errors
		if (status != 0) {
			fprintf(stderr, "Warning: 'dsyev' failed, info: %i\n", status);
			err.codes |= E_LAPACK;
		}
		if (eigenval[0] < 0) {
			fprintf(stderr, "Warning: Hessian has negative eigenvalue %g.\n", eigenval[0]);
			err.codes |= E_NEG_HESS;
			err.hess_eig = eigenval[0];
		}

		free(eigenval);
		free(hesssubm);
		free(hessfull);
	}

	return err;
}
*/