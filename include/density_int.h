//	Calculate integral Ne(r) of the density rho(r), as well as the inverse function of Ne.
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

#ifndef DENSITY_INT_H
#define DENSITY_INT_H


//_______________________________________________________________________________________________________________________
///
/// \brief Data structure for the density integral function 'Ne(r)' and its inverse
///
typedef struct
{
	int nelec;			//!< number of electrons
	int nshpt;			//!< number of radial grid points for each shell
	double *rgrid;		//!< radial grid points Ne^{-1}(i + j/nshpt), j = 0, 1, ..., nshpt-1 for each shell
}
dataNe_t;


//_______________________________________________________________________________________________________________________
//


int IntegrateDensity(const int nelec, const double dr, const int ngrid, const double *rho, const int nshpt, dataNe_t *data);

void FreeDensityInt(dataNe_t *data);



#endif
