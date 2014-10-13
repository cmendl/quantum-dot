//	Evaluate the radial co-motion functions f_i.
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

#ifndef COMOTION_H
#define COMOTION_H

#include "density_int.h"


// retrieve the grid index of electron 'k', given the grid index 'ir' of the first electron
int CoMotionIndexMap(const int nelec, const int nshpt, const int ir, const int k);

// assuming that f points to an array of 'nelec' doubles
void EvaluateCoMotion(const dataNe_t *data, const int shellindex, const int pointindex, double *f);

// return the shell of electron 'k' for first electron in shell 'shellindex'
int GetCoMotionShell(const int nelec, const int shellindex, const int k);

// return the index of the electron at infinity, or -1 if all electron radii are finite
int GetInfiniteCoMotionShell(const int nelec, const int shellindex);

// return the index of the electron at the origin, or -1 if all electron radii are larger than zero
int GetZeroCoMotionShell(const int shellindex);



#endif
