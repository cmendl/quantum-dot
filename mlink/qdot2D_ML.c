/// \file qdot2D_ML.c
///	\brief Main file of the MathLink Mathematica interface.
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

#include "comotion.h"
#include "angles.h"
#include "scf_step.h"
#include "error_desc.h"
#include <mathlink.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>



// similar to 'MLPutReal64Array', workaround to treat 1.#INF correctly
static void MLPutReal64Matrix(MLINK mlp, double *data, int dims[2])
{
	int i;

	MLPutFunction(mlp, "List", dims[0]);
	for (i = 0; i < dims[0]; i++)
	{
		MLPutReal64List(mlp, &data[i*dims[1]], dims[1]);
	}
}


//_______________________________________________________________________________________________________________________
//


void MLQuantumDotStep(int nelec, double omega, double *rho, long rho_len, double dr, int nshpt, double kBT, double *phiStart, long phiStart_len, int spinpol)
{
	int i;

	// check input arguments
	if (nelec <= 0 || omega <= 0 || rho_len <= 0 || dr <= 0 || nshpt <= 0 || kBT < 0 || phiStart_len != nelec)
	{
		MLEvaluate(stdlink, "Message[MLQuantumDotStep::InvalidArg]");
		// discard 'ReturnPacket'
		MLNextPacket(stdlink);
		MLNewPacket(stdlink);	// discard
		// final output
		MLPutSymbol(stdlink, "$Failed");
		return;
	}

	// perform computation
	potential_t pot;
	en_occlist_t occlist;
	double *rho_next = (double *)malloc(rho_len*sizeof(double));
	error_desc_t err = SelfConsistentFieldStep(nelec, omega, rho, dr, rho_len, nshpt, spinpol, kBT, phiStart, &pot, &occlist, rho_next);

	// report errors
	if ((err.codes & E_DENSITY_INT) != 0)
	{
		MLEvaluate(stdlink, "Message[MLQuantumDotStep::DensityIntError]");
		// discard 'ReturnPacket'
		MLNextPacket(stdlink);
		MLNewPacket(stdlink);	// discard
	}
	if ((err.codes & E_LAPACK) != 0)
	{
		MLEvaluate(stdlink, "Message[MLQuantumDotStep::LapackError]");
		// discard 'ReturnPacket'
		MLNextPacket(stdlink);
		MLNewPacket(stdlink);	// discard
	}
	if ((err.codes & E_ITERATION) != 0)
	{
		char buf[1024];
		sprintf(buf, "Message[MLQuantumDotStep::Iteration, \"%s\"]", err.msg);
		MLEvaluate(stdlink, buf);
		// discard 'ReturnPacket'
		MLNextPacket(stdlink);
		MLNewPacket(stdlink);	// discard
	}
	if ((err.codes & E_CONVERGENCE) != 0)
	{
		char buf[1024];
		sprintf(buf, "Message[MLQuantumDotStep::Convergence, %.16f]", err.conv_delta);
		MLEvaluate(stdlink, buf);
		// discard 'ReturnPacket'
		MLNextPacket(stdlink);
		MLNewPacket(stdlink);	// discard
	}
	if ((err.codes & E_NEG_HESS) != 0)
	{
		char buf[1024];
		sprintf(buf, "Message[MLQuantumDotStep::NegHessian, %.16f]", err.hess_eig);
		MLEvaluate(stdlink, buf);
		// discard 'ReturnPacket'
		MLNextPacket(stdlink);
		MLNewPacket(stdlink);	// discard
	}

	// return results to Mathematica
	MLPutFunction(stdlink, "List", 8);

	// new density
	MLPutReal64List(stdlink, rho_next, rho_len);

	// energy level occupation list
	MLPutFunction(stdlink, "List", occlist.length);
	for (i = 0; i < occlist.length; i++)
	{
		MLPutFunction(stdlink, "List", 3);

		MLPutReal64(stdlink, occlist.occupancy[i].en);

		MLPutFunction(stdlink, "List", 2);
		MLPutInteger(stdlink, occlist.occupancy[i].n);
		MLPutInteger(stdlink, occlist.occupancy[i].m);

		MLPutReal64(stdlink, occlist.occupancy[i].occ);
	}

	// SCE potential and derivative; omit last point at infinity
	MLPutReal64List(stdlink, pot.v,  nelec*nshpt);
	MLPutReal64List(stdlink, pot.dv, nelec*nshpt);

	// shell radii a_i
	double *a = malloc((nelec+1)*sizeof(double));
	for (i = 0; i <= nelec; i++) {
		a[i] = pot.rgrid[i*nshpt];
	}
	MLPutReal64List(stdlink, a, nelec+1);
	free(a);

	// radial co-motion functions and corresponding optimized angles
	int dims[2];
	dims[0] = nelec*nshpt;
	dims[1] = nelec;
	MLPutReal64Matrix(stdlink, pot.f,   dims);
	MLPutReal64Matrix(stdlink, pot.phi, dims);

	// g function
	MLPutReal64List(stdlink, pot.g, nelec*nshpt);

	MLEndPacket(stdlink);
	MLFlush(stdlink);

	// clean up
	FreeOccupancy(&occlist);
	FreePotential(&pot);
	free(rho_next);
}


//_______________________________________________________________________________________________________________________
//


void MLQuantumDotAngleOpt(int nelec, double omega, double *rho, long rho_len, double dr, int nshpt, double *phiStart, long phiStart_len)
{
	// check input arguments
	if (nelec <= 0 || omega <= 0 || rho_len <= 0 || dr <= 0 || nshpt <= 0 || phiStart_len != nelec || phiStart[0] != 0)
	{
		MLEvaluate(stdlink, "Message[MLQuantumDotAngleOpt::InvalidArg]");
		// discard 'ReturnPacket'
		MLNextPacket(stdlink);
		MLNewPacket(stdlink);	// discard
		// final output
		MLPutSymbol(stdlink, "$Failed");
		return;
	}

	error_desc_t err = { 0 };

	// integrate density
	dataNe_t data;
	err.codes |= IntegrateDensity(nelec, dr, rho_len, rho, nshpt, &data);

	// current 'f' and 'phi'
	double *f   = malloc(nelec*sizeof(double));	// radial part of co-motion functions
	double *phi = malloc(nelec*sizeof(double));	// polar angles

	// obtain radial parts
	EvaluateCoMotion(&data, 0, nshpt/2, f);
	memcpy(phi, phiStart, nelec*sizeof(double));

	// perform global angle optimization
	double optVee;
	error_desc_t cur_err = OptimizeAnglesGlobal(nelec, f, phi, &optVee);
	AccumulateError(&err, &cur_err);

	// report errors
	if ((err.codes & E_DENSITY_INT) != 0)
	{
		MLEvaluate(stdlink, "Message[MLQuantumDotAngleOpt::DensityIntError]");
		// discard 'ReturnPacket'
		MLNextPacket(stdlink);
		MLNewPacket(stdlink);	// discard
	}
	if ((err.codes & E_LAPACK) != 0)
	{
		MLEvaluate(stdlink, "Message[MLQuantumDotAngleOpt::LapackError]");
		// discard 'ReturnPacket'
		MLNextPacket(stdlink);
		MLNewPacket(stdlink);	// discard
	}
	if ((err.codes & E_ITERATION) != 0)
	{
		char buf[1024];
		sprintf(buf, "Message[MLQuantumDotAngleOpt::Iteration, \"%s\"]", err.msg);
		MLEvaluate(stdlink, buf);
		// discard 'ReturnPacket'
		MLNextPacket(stdlink);
		MLNewPacket(stdlink);	// discard
	}
	if ((err.codes & E_CONVERGENCE) != 0)
	{
		char buf[1024];
		sprintf(buf, "Message[MLQuantumDotAngleOpt::Convergence, %.16f]", err.conv_delta);
		MLEvaluate(stdlink, buf);
		// discard 'ReturnPacket'
		MLNextPacket(stdlink);
		MLNewPacket(stdlink);	// discard
	}
	if ((err.codes & E_NEG_HESS) != 0)
	{
		char buf[1024];
		sprintf(buf, "Message[MLQuantumDotAngleOpt::NegHessian, %.16f]", err.hess_eig);
		MLEvaluate(stdlink, buf);
		// discard 'ReturnPacket'
		MLNextPacket(stdlink);
		MLNewPacket(stdlink);	// discard
	}

	// return results to Mathematica
	MLPutFunction(stdlink, "List", 2);
	// optimized angles
	MLPutReal64List(stdlink, phi, nelec);
	MLPutReal64(stdlink, optVee);

	MLEndPacket(stdlink);
	MLFlush(stdlink);

	// clean up
	free(phi);
	free(f);
	FreeDensityInt(&data);
}


//_______________________________________________________________________________________________________________________
//


#if MACINTOSH_MATHLINK

int main(int argc, char* argv[])
{
	/* Due to a bug in some standard C libraries that have shipped with
	 * MPW, zero is passed to MLMain below.  (If you build this program
	 * as an MPW tool, you can change the zero to argc.)
	 */
	argc = argc; /* suppress warning */
	return MLMain(0, argv);
}

#elif WINDOWS_MATHLINK

#if __BORLANDC__
#pragma argsused
#endif

int PASCAL WinMain(HINSTANCE hinstCurrent, HINSTANCE hinstPrevious, LPSTR lpszCmdLine, int nCmdShow)
{
	char  buff[512];
	char FAR * buff_start = buff;
	char FAR * argv[32];
	char FAR * FAR * argv_end = argv + 32;

	hinstPrevious = hinstPrevious; /* suppress warning */

	if (!MLInitializeIcon(hinstCurrent, nCmdShow)) return 1;
	MLScanString(argv, &argv_end, &lpszCmdLine, &buff_start);
	return MLMain((int)(argv_end - argv), argv);
}

#else

int main(int argc, char* argv[])
{
	return MLMain(argc, argv);
}

#endif
