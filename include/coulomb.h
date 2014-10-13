//	Calculate pairwise Coulomb repulsion between electrons,
//	and corresponding gradients with respect to polar coordinates.
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

#ifndef COULOMB_H
#define COULOMB_H


//_______________________________________________________________________________________________________________________
//
// distance in polar coordinates
//
double DistancePolar(const double r0, const double r1, const double dphi);

// derivative of 1/distance with respect to r0
double InvDistanceDr(const double r0, const double r1, const double dphi);

// derivative of 1/distance with respect to delta phi
double InvDistanceD1dphi(const double r0, const double r1, const double dphi);

// second derivative of 1/distance with respect to delta phi
double InvDistanceD2dphi(const double r0, const double r1, const double dphi);


//_______________________________________________________________________________________________________________________
//
// pairwise Coulomb repulsion of electrons with positions (r,phi) in polar coordinates
//
double Vee(const int nelec, const double *r, const double *phi);

// derivative with respect to r[0]
double VeeDr(const int nelec, const double *r, const double *phi);

// gradient with respect to angles 'phi'
void VeeGradPhi(const int nelec, const double *r, const double *phi, double *grad);

// Hessian matrix (second derivative) with respect to angles 'phi'
// 'hess' must point to an array of size 'nelec * nelec'
void VeeHessPhi(const int nelec, const double *r, const double *phi, double *hess);



#endif
