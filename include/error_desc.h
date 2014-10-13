//	Internal error codes and description.
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

#ifndef ERROR_DESC_H
#define ERROR_DESC_H

#include "util.h"


/// \brief error codes
enum error_codes_t
{
	E_NO_ERROR		= 0,
	E_DENSITY_INT	= 1 << 0,
	E_LAPACK		= 1 << 1,
	E_ITERATION		= 1 << 2,
	E_CONVERGENCE	= 1 << 3,
	E_NEG_HESS		= 1 << 4,
};


//_______________________________________________________________________________________________________________________
//


/// \brief detailed description of accumulated errors
typedef struct
{
	int codes;				//!< error codes, accumulated by bitwise OR
	const char *msg;		//!< error message
	double conv_delta;		//!< last delta in Newton-iteration, not converged if delta >= tol
	double hess_eig;		//!< negative eigenvalue of Hessian matrix
}
error_desc_t;


static inline void AccumulateError(error_desc_t *desc, const error_desc_t *desc_add)
{
	// collect error codes
	desc->codes |= desc_add->codes;

	// not optimal, cannot handle multiple (different) messages
	if (desc_add->msg != 0) {
		desc->msg = desc_add->msg;
	}

	// take larger delta (should be zero if there is no convergence error)
	desc->conv_delta = maxf(desc->conv_delta, desc_add->conv_delta);

	// take smaller (more negative) eigenvalue (should be zero if there is no negative eigenvalue error)
	desc->hess_eig = minf(desc->hess_eig, desc_add->hess_eig);
}



#endif
