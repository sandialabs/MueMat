% Demonstrate how to use the MATLAB interface of ML (MLMEX)
%
% ML is the multilevel preconditioning package in Trilinos.
%
% MLMEX is ML's interface to the MATLAB environment. It allows
% access to ML's aggregation and solver routines from MATLAB.
%
% For more documentation, see also:
% ML 5.0 Smoothed Aggregation User's Guide
% http://trilinos.sandia.gov/packages/ml/mlguide5.pdf
%
% --
%
% HOW TO INSTALL MLMEX:
%
% Download Trilinos: http://trilinos.sandia.gov/download
%
% ML supports MATLAB 7.2 (R2006a) and up. For Matlab 64-bit
% (glnxa64), ML and all required Trilinos libraries must be
% compiled with the -fPIC option. BLAS and LAPACK libraries must
% also be compiled with this flag. Using static linking for BLAS
% and LAPACK prevents MATLAB's default libraries to take precedence.
%
% Script to build ML with MATLAB support :
% (You just need to change path to BLAS, LAPACK and MATLAB.)
%
%   #!/bin/sh
%
%   # cmake script for Trilinos.
%
%   export TRILINOS_HOME=/home/jhgaida/Trilinos/dev/Trilinos
%   EXTRA_ARGS=$@
%
%   rm CMakeCache.txt
%
%   cmake \
%     -D Trilinos_ENABLE_ML:BOOL=ON \
%   \
%     -D ML_ENABLE_Teuchos=ON       \
%     -D ML_ENABLE_Epetra=ON        \
%     -D ML_ENABLE_Galeri=ON        \
%     -D ML_ENABLE_Amesos=ON        \
%     -D ML_ENABLE_Ifpack=ON        \
%     -D ML_ENABLE_AztecOO=ON       \
%     -D ML_ENABLE_EpetraExt=ON     \
%   \
%     -D TPL_BLAS_LIBRARIES:STRING="/path-to/lapack-3.2.1/blas_LINUX.a"     \
%     -D TPL_LAPACK_LIBRARIES:STRING="/path-to/lapack-3.2.1/lapack_LINUX.a" \
%   \
%     -D TPL_ENABLE_MATLAB:BOOL=ON \
%     -D MATLAB_ROOT:STRING="/path-to/matlab" \
%     -D MATLAB_ARCH:STRING="glnxa64" \
%   \
%     -D CMAKE_CXX_FLAGS:STRING="-fPIC" \
%     -D CMAKE_C_FLAGS:STRING="-fPIC" \
%     -D CMAKE_Fortran_FLAGS:STRING="-fPIC" \
%   \
%     $EXTRA_ARGS \
%     ${TRILINOS_HOME}
%
% ML is interfaced with MATLAB via the MATLAB script ml.m. So add its
% folder to the search path of MATLAB:
% > addpath('/home/.../Trilinos-build/packages/ml/matlab')
%
% More installation documentation is available in ML 5.0 Smoothed
% Aggregation User's Guide.
%
% See also: ElasticityTest

if ~exist('ml','file'), fprintf('This example requires MLMEX. Type ''doc MLInterface'' for more information.\n'); return; end

clear all;

%% Define a Poisson problem
n  = 30;
A  = BuildLaplace2D(n);
b  = rand(A.GetRowMap().NDOFs(),1);
ns = BuildNullSpace(A);
ns = full(ns);

%% ML's aggregation mode.
% Visualize the result of ML's aggregation.
agg=ml('aggregate',A.GetMatrixData());

NODES(:,1) = reshape(meshgrid(1:n,1:n),   1, []);    % = [ 1111 2222 3333 4444 ]
NODES(:,2) = reshape(meshgrid(1:n,1:n).', 1, []).';  % = [ 1234 1234 1234 1234 ]
plot(NODES(:,1),NODES(:,2),'x')
dx=0.3;
h=plot_aggregate_boxes(NODES,agg+1,dx);

%% ML's Setup mode.
% Forms the problem setup for ML.This call returns a problem handle
% used to reference the problem in the future.
%
% See also: ElasticityTest for more advanced usage

% ML Smoothed Aggregation,
% Options: symmetric Gauss-Seidel, direct solver on coarsest level.
[h,OC]=ml('setup',A.GetMatrixData(),'PDE equations',1,...
         'null space: type','pre-computed','null space: dimension',1,'null space: vectors', ns,...
         'smoother: type','symmetric Gauss-Seidel', 'smoother: sweeps',2,'coarse: max size',100,'max levels',20,...
         'coarse: type','Amesos-KLU', 'ML output',10);

fprintf('Operator Complexity: %3.2f\n',OC);

%% ML's Solve mode.
% Given a problem handle and a right-hand side, ML solves the
% problem specified. On this example, ML is used as a preconditionner for CG.
% For more advanced usage, see also experiment/ElasticityTest
maxIts = 99;
tol    = 1e-12;

[x,flag,relRes,nIts,resVec] = pcg(A.GetMatrixData(),b,tol,maxIts,...
                              @(b)ml(h, ...
                                    A.GetMatrixData(), b, ...
                                    'krylov: type','fixed point',...
                                    'krylov: max iterations',1,'krylov: tolerance',1e-100,...
                                    'krylov: output level',0));
if flag == 0
  fprintf('pcg converged to the desired tolerance %g within %d iterations.\n', tol, nIts);
  fprintf('nIts=%d - relative residual norm=%g\n', nIts, relRes);
else
  error('Laplace.m:Solve', 'pcg flag');
end

%% Clean up
ml('cleanup');
