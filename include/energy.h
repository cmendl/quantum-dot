//	Sort and fill energy levels to obtain occupancy of each level.
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

#ifndef ENERGY_H
#define ENERGY_H

#include "util.h"


/// \brief energy levels, for each 'm' quantum number
typedef struct
{
	double *en;		//!< array of size (m_max+1) x num
	int m_max;		//!< maximum 'm' quantum number (angular momentum)
	int num;		//!< number of energy levels for fixed 'm'
}
en_levels_t;


void AllocateEnergyLevels(const int m_max, const int num, en_levels_t *levels);

void FreeEnergyLevels(en_levels_t *levels);


//_______________________________________________________________________________________________________________________
//


/// \brief Occupancy of an individual energy level
typedef struct
{
	double en;		//!< energy value
	int m;			//!< 'm' quantum number
	int n;			//!< energy quantum number
	double occ;		//!< occupancy of this level (including degeneracy); allow fractional occupancies
}
en_occupancy_t;


/// \brief Occupancy list of energy levels
typedef struct
{
	int length;					//!< list length
	en_occupancy_t *occupancy;	//!< pointer to array of occupancies
}
en_occlist_t;


void CalculateOccupancy(const en_levels_t *levels, const int nelec, const bool spinpol, const double kBT, en_occlist_t *occlist);

void FreeOccupancy(en_occlist_t *occlist);



#endif
