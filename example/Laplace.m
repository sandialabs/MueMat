% A Laplace example solves via a 2-level AMG hierarchy.
%
% Use a symmetric GaussSeidel smoother and
% AMG as a preconditioner to conjugate gradients.
%
% See also: BuildLaplace1D, BuildLaplace2D, BuildLaplace2DStretchedMesh
%
function Laplace()

  clear all;
  srand;

  % Solve a Laplace 1D problem
  n      = 100;
  A      = BuildLaplace1D(n);
  b      = rand(A.GetRowMap().NDOFs(),1);
  x      = Solve(A,b);

  % Solve a Laplace 2D problem
  n      = 30;
  A      = BuildLaplace2D(n);
  b      = rand(A.GetRowMap().NDOFs(),1);
  x      = Solve(A,b);

  % Solve a Laplace 2D (stretched mesh) problem
  n      = 30;
  epsi   = 1.e-3;
  A      = BuildLaplace2DStretchedMesh(n,epsi);
  b      = rand(A.GetRowMap().NDOFs(),1);
  x      = Solve(A,b);

end % function Laplace

function [x, nIts, relRes] = Solve(A, b)
      %SOLVE Solve Ax = b using CG preconditioned by AMG.
      %
      %   SYNTAX   [x, nIts, relRes] = Solve(A, b);
      %
      %     A         - matrix of the problem          (Operator)
      %     b         - right hand side                (MATLAB vector or MultiVector)
      %     x         - solution                       (same as b)
      %     nIts      - number of iterations performed (integer)
      %     relRes    - relative residual norm         (double)

  mue_include

  nDOFS = A.GetRowMap.NDOFs();

  %% Set AMG options
  numDesiredLevels = 2; % maximum number of AMG level

  % Setup prolongator factory
  AmalgamateDropFact = CoalesceDropFactory();
  AggFact            = AggregationFactory();
  %AggFact           = AggFact.SetAlgorithm('graph');
  Ptentfact          = TentativePFactory(AmalgamateDropFact,AggFact);
  Pfact              = SaPFactory(Ptentfact);

  % Setup restrictor factory
  Rfact              = TransPFactory();

  % combine prolongator and restrictor
  PRfact             = GenericPRFactory(Pfact,Rfact);
  PRfact.SetMaxCoarseSize(30);

  % Setup how to build coarse discretization
  Acfact             = RAPFactory(); % A_coarse = R * A * P

  % Setup smoother factory
  SmooFactory        = SmootherFactory(Smoother('GaussSeidel', 2, 1));

  %% Create the multigrid hierarchy
  % Construct and populate finest level with user information
  Finest = Level();
  Finest.Set('A', A);
  Finest.Set('NullSpace', BuildNullSpace(A));

  % Create a multigrid hierarchy
  MgHierarchy = Hierarchy();
  MgHierarchy.SetOutputLevel(1);
  MgHierarchy.SetLevel(Finest,1);
  MgHierarchy.FillHierarchy(PRfact, [], Acfact, 1, numDesiredLevels);
  MgHierarchy.SetSmoothers(SmooFactory);

  %% Iterative solve

  % CG parameters
  maxIts = 9;     % Maximum number of iterations
  tol    = 1e-12; % Tolerance parameter

  % AMG as a preconditioner to CG
  x    = zeros(nDOFS,1); SolStatus = ALLZEROS;
  zero = zeros(nDOFS,1);
  [x, flag, relRes, nIts, resVec] = pcg(A.GetMatrixData(),b,tol,maxIts,@(b)MgHierarchy.Iterate(b,1, zero,ALLZEROS)); SolStatus = NOTALLZEROS;

  % Test output of pcg
  if flag == 0
    fprintf('pcg converged to the desired tolerance %g within %d iterations.\n', tol, nIts);
    fprintf('nIts=%d - relative residual norm=%g\n', nIts, relRes);
    % relRes == norm(b-A.GetMatrixData()*x)/norm(b);
  elseif flag == 1
    fprintf('pcg iterated %d times but did not converge.\n', maxIts);
    % maxIts == nIts;
  elseif flag == 2
    fprintf('Preconditioner M was ill-conditioned.\n');
  elseif flag == 3
    fprintf('pcg stagnated. (Two consecutive iterates were the same.)\n');
  elseif flag == 4
    fprintf('One of the scalar quantities calculated during pcg became too small or too large to continue computing.\n');
  else
    error('Laplace.m:Solve', 'unknow pcg flag');
  end

end % function Solve
