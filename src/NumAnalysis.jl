module NumAnalysis

using Printf
using FastGaussQuadrature

include("calculus.jl")
export simpson,
       trapezoid,
       adaptivequad,
       gaussianquad

include("direct_lin_algebra.jl")
export gaussbacksub,
       gausspivotpartial,
       gausspivotscaled,
       matrixfactorLU,
       matrixfactorPtLU,
       matrixfactorLDLt,
       matrixfactorLLt,
       systemsolveLU

include("eigen_estimate.jl")
export eigenpower,
       eigeninvpower,
       householder,
       hessenberg,
       qreigen

include("init_val_ode.jl")
export ivpsolve,
       euler,
       rungekuttaO4,
       midpoint

include("interpolate.jl")
export neville,
       newtondivdiff,
       hermite,
       natcubicspline,
       clampcubicspline

include("iter_lin_algebra.jl")
export jacobi,
       gaussseidel

include("lin_algebra.jl")
export gaussjordanelim,
       isposdef,
       issymmetric,
       istridiag,
       buildtridiag,
       IMatrix,
       typecheckfloat,
       swaprow,
       swapcol,
       addrow,
       reducecol,
       augmentmatrix,
       forwardsub,
       backsub,
       innerprod,
       maxlist,
       norm,
       norminfvector,
       domcomponent,
       norminfmatrix

include("nonlin_solve.jl")
export passto,
       fixedpointNL,
       newtonsNL

include("root_finder.jl")
export bisecmethod,
       fixedpoint,
       newtons,
       secantmethod,
       steffensen,
       horners

end