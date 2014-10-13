//	One self-consistent field (SCF) iteration step.
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

#ifndef SCF_STEP_H
#define SCF_STEP_H

#include "potential.h"
#include "energy.h"
#include "util.h"


error_desc_t SelfConsistentFieldStep(const int nelec, const double omega, const double *rho, const double dr, const int ngrid, const int nshpt,
	const bool spinpol, const double kBT, const double *phiStart, potential_t *pot, en_occlist_t *occlist, double *rho_next);



#endif
