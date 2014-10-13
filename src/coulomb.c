/// \file coulomb.c
/// \brief Calculate pairwise Coulomb repulsion between electrons,
/// and corresponding gradients with respect to polar coordinates.
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

#include "coulomb.h"
#include "util.h"
#include <math.h>
#include <assert.h>


//_______________________________________________________________________________________________________________________
///
/// \brief Distance in polar coordinates
///
/// \param r0 radius of first particle
/// \param r1 radius of second particle
/// \param dphi relative angle between particles
double DistancePolar(const double r0, const double r1, const double dphi)
{
	if (r0 > r1)
	{
		double rd = r1/r0;
		return r0*sqrt(1.0 + rd*(rd - 2*cos(dphi)));
	}
	else
	{
		double rd = r0/r1;
		return r1*sqrt(1.0 + rd*(rd - 2*cos(dphi)));
	}
}

/// \brief Derivative of 1/distance with respect to r0
double InvDistanceDr(const double r0, const double r1, const double dphi)
{
	// save one evaluation of 'cos'
	double cosdphi = cos(dphi);
	double d;

	if (r0 > r1)
	{
		double rd = r1/r0;
		d = r0*sqrt(1.0 + rd*(rd - 2*cosdphi));
	}
	else
	{
		double rd = r0/r1;
		d = r1*sqrt(1.0 + rd*(rd - 2*cosdphi));
	}
	assert(d == DistancePolar(r0, r1, dphi));

	return (r1*cosdphi - r0)/(d*d*d);
}


/// \brief Derivative of 1/distance with respect to delta phi
double InvDistanceD1dphi(const double r0, const double r1, const double dphi)
{
	double d = DistancePolar(r0, r1, dphi);

	return -r0*r1*sin(dphi)/(d*d*d);
}


/// \brief Second derivative of 1/distance with respect to delta phi
double InvDistanceD2dphi(const double r0, const double r1, const double dphi)
{
	double d, cosdphi;

	if (r0 == 0 || r1 == 0)
	{
		return 0;
	}

	d = DistancePolar(r0, r1, dphi);
	cosdphi = cos(dphi);

	return -square(r0*r1)*(cosdphi*(r0/r1 + r1/r0) + square(cosdphi) - 3)/(d*square(square(d)));
}


//_______________________________________________________________________________________________________________________
//


/// \brief Pairwise Coulomb repulsion of electrons with positions (r,phi) in polar coordinates
double Vee(const int nelec, const double *r, const double *phi)
{
	int i, j;
	double v = 0;

	for (i = 0; i < nelec; i++)
	{
		for (j = i+1; j < nelec; j++)
		{
			v += 1/DistancePolar(r[i], r[j], phi[i] - phi[j]);
		}
	}

	return v;
}


/// \brief Derivative with respect to r[0]
double VeeDr(const int nelec, const double *r, const double *phi)
{
	int i;
	double dv = 0;

	// assuming that phi[0] == 0
	assert(phi[0] == 0);

	for (i = 1; i < nelec; i++)
	{
		dv += InvDistanceDr(r[0], r[i], phi[i]);
	}

	return dv;
}


/// \brief Gradient with respect to angles 'phi'
void VeeGradPhi(const int nelec, const double *r, const double *phi, double *grad)
{
	int i, j;

	for (i = 0; i < nelec; i++)
	{
		grad[i] = 0;

		for (j = 0; j < i; j++)
		{
			grad[i] += InvDistanceD1dphi(r[i], r[j], phi[i] - phi[j]);
		}
		for (j = i+1; j < nelec; j++)
		{
			grad[i] += InvDistanceD1dphi(r[i], r[j], phi[i] - phi[j]);
		}
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Hessian matrix (second derivative) with respect to angles 'phi'
///
/// \param nelec number of electrons
/// \param r radii
/// \param phi angles
/// \param hess return value; must point to an array of size 'nelec * nelec'
void VeeHessPhi(const int nelec, const double *r, const double *phi, double *hess)
{
	int i, j;

	for (i = 0; i < nelec; i++)
	{
		// diagonal entry
		double h = 0;
		for (j = 0; j < i; j++)
		{
			h += InvDistanceD2dphi(r[i], r[j], phi[i] - phi[j]);
		}
		for (j = i+1; j < nelec; j++)
		{
			h += InvDistanceD2dphi(r[i], r[j], phi[i] - phi[j]);
		}
		hess[i*nelec+i] = h;

		// off-diagonal entries
		for (j = i+1; j < nelec; j++)
		{
			hess[i*nelec+j] = hess[j*nelec+i] = -InvDistanceD2dphi(r[i], r[j], phi[i] - phi[j]);
		}
	}
}
