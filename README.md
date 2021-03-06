Quantum dots in two dimensions and the "strictly correlated electrons" (SCE) framework
======================================================================================

C source code (with minor improvements) and Mathematica demonstration file `quantum_dot_demo.nb` for the calculations in Ref. 1: Kohn-Sham density functional theory using the "strictly correlated electrons" (SCE) functional applied to quantum dots in two dimensions. The program uses the Mathematica MathLink interface and requires a BLAS/LAPACK library, for example AMD's Core Math Library (ACML), as well as the GNU Scientific Library (GSL).

How to compile the source code:
* Windows: Visual Studio 2013 project files are provided in the *vcproj_mlink* subfolder
* Linux: a makefile is available in the *mlink* subfolder. You probably have to adapt paths and the BLAS/LAPACK function calls to your local installation


License
-------
Copyright (c) 2013-2014, Christian B. Mendl  
All rights reserved.  
http://christian.mendl.net

This program is free software; you can redistribute it and/or
modify it under the terms of the Simplified BSD License
http://www.opensource.org/licenses/bsd-license.php


References
----------
1. Christian B. Mendl, Francesc Malet, Paola Gori-Giorgi  
   Wigner localization in quantum dots from Kohn-Sham density functional theory without symmetry breaking  
   Phys. Rev. B 89, 125106 (2014), [arXiv:1311.6011](http://arxiv.org/abs/1311.6011)
2. Michael Seidl, Paola Gori-Giorgi, Andreas Savin  
   Strictly correlated electrons in density-functional theory: A general formulation with applications to spherical densities  
   Phys. Rev. A 75, 042511 (2007), [arXiv:cond-mat/0701025](http://arxiv.org/abs/cond-mat/0701025)
