//	Optimize polar angles with respect to Coulomb repulsion, given the radii.
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

#ifndef ANGLES_H
#define ANGLES_H

#include "error_desc.h"


error_desc_t OptimizeAngles(const int nelec, const double *r, double *phi);

error_desc_t OptimizeAnglesGlobal(const int nelec, const double *r, double *phi, double *optVee);



#endif
