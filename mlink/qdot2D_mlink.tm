
:Begin:
:Function:		MLQuantumDotStep
:Pattern:		MLQuantumDotStep[nelec_Integer, omega_Real, rho_List, dr_Real, nshpt_Integer, kBT_Real, phiStart_List, spinpol_Integer]
:Arguments:		{ nelec, omega, rho, dr, nshpt, kBT, phiStart, spinpol }
:ArgumentTypes:	{ Integer, Real, RealList, Real, Integer, Real, RealList, Integer }
:ReturnType:	Manual
:End:

:Evaluate: MLQuantumDotStep::usage = "MLQuantumDotStep[nelec_Integer, omega_Real, rho_List, dr_Real, nshpt_Integer, kBT_Real, phiStart_List, spinpol_Integer] performs a self-consistent field iteration step for a quantum dot model (2D harmonic oscillator) in the SCE formalism."
:Evaluate: MLQuantumDotStep::InvalidArg = "invalid argument"
:Evaluate: MLQuantumDotStep::OutOfMemory = "malloc failed, probably out of memory"
:Evaluate: MLQuantumDotStep::DensityIntError = "integral of density 'Ne' not strictly increasing"
:Evaluate: MLQuantumDotStep::LapackError = "an internal LAPACK function failed"
:Evaluate: MLQuantumDotStep::Iteration = "angle optimization iteration canceled, error message: `1`"
:Evaluate: MLQuantumDotStep::Convergence = "angle optimization not converged within maximum number of iterations, final delta = `1`"
:Evaluate: MLQuantumDotStep::NegHessian = "angle optimization: Hessian matrix has the negative eigenvalue `1`"



:Begin:
:Function:		MLQuantumDotAngleOpt
:Pattern:		MLQuantumDotAngleOpt[nelec_Integer, omega_Real, rho_List, dr_Real, nshpt_Integer, phiStart_List]
:Arguments:		{ nelec, omega, rho, dr, nshpt, phiStart }
:ArgumentTypes:	{ Integer, Real, RealList, Real, Integer, RealList }
:ReturnType:	Manual
:End:

:Evaluate: MLQuantumDotAngleOpt::usage = "MLQuantumDotAngleOpt[nelec_Integer, omega_Real, rho_List, dr_Real, nshpt_Integer, phiStart_List] optimize anglar part of co-motion functions for a quantum dot model (2D harmonic oscillator) in the SCE formalism."
:Evaluate: MLQuantumDotAngleOpt::InvalidArg = "invalid argument"
:Evaluate: MLQuantumDotAngleOpt::OutOfMemory = "malloc failed, probably out of memory"
:Evaluate: MLQuantumDotAngleOpt::DensityIntError = "integral of density 'Ne' not strictly increasing"
:Evaluate: MLQuantumDotAngleOpt::LapackError = "an internal LAPACK function failed"
:Evaluate: MLQuantumDotAngleOpt::Iteration = "angle optimization iteration canceled, error message: `1`"
:Evaluate: MLQuantumDotAngleOpt::Convergence = "angle optimization not converged within maximum number of iterations, final delta = `1`"
:Evaluate: MLQuantumDotAngleOpt::NegHessian = "angle optimization: Hessian matrix has the negative eigenvalue `1`"
