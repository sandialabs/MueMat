% A Laplace example that tests a fairly special smoother on the
% finest level. In particular, the fine level pre smoother
% corresponds to
%              Smoother('GaussSeidel', 1, 1.) 
%           TipSmoother('GaussSeidel', 1, 1.)
%
% while the fine level post smoother is 
%
%           TipSmoother('GaussSeidel', 1, 1.)
%              Smoother('GaussSeidel', 1, 1.) 
%
% TipSmoother is actually a block Gauss-Seidel with dofs
% making up blocks determined by a user supplied Collection
% which is put into 'Tips' in the Level structure.
%
function MergedSmootherWithSpecialBlkGS()

  clear all;
  srand;

  % Solve a Laplace 2D problem
  n      = 60; 
  A      = BuildLaplace2D(n);
  b      = rand(A.GetRowMap().NDOFs(),1);

  mue_include

  nDOFS = A.GetRowMap.NDOFs();

  %% Set AMG options
  numDesiredLevels = 3; % maximum number of AMG level

  % Setup prolongator factory
  %AmalgamateDropFact = CoalesceDropFactory();
  %AggFact            = AggregationFactory();
  Pfact              = SaPFactory();

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
  MgHierarchy.FillHierarchy(PRfact,[], Acfact, 1, numDesiredLevels);
  if (MgHierarchy.GetNumLevel() > 1)
    MgHierarchy.SetSmoothers(SmooFactory, 2, MgHierarchy.GetNumLevel()-2);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Special Smoother does a combination of regular Gauss-Seidel         %
  % and block Gauss-Seidel where blocks correspond to tip dofs          %
  % (which are just faked as the following dofs: 1, 3, 7). The          %
  % main trick is that the post smoother must reverse the order         %
  % in which the regular and block Gauss-Seidel are applied.            %
  Collection.NSubsets = 1;                                              %
  Collection.Subsets(1)=CreateDOFSubset(A.GetRowMap(),'Scattered',-1,-1,...
                                          [1 3 7]);  % made up tip dofs %
  Finest.Set('Tips', Collection);  % TipSmoother looks for this.        %
  tipsmo = TipSmoother('GaussSeidel', 1, 1.);                           %
  smothers = {Smoother('GaussSeidel', 1, 1.) tipsmo};                   %
  MergedPre = MergedSmoother( smothers );                               %
  MergedPost= MergedPre.Copy();                                         %
  MergedPost.ReverseOrder();                                            %
  NewSmooFactory = SmootherFactory(MergedPre,MergedPost);               %
  MgHierarchy.SetSmoothers(NewSmooFactory, 1, 1);  % set on finest level%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

