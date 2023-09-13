% Does not work?

% A 1D Laplace example solves via a 2-level AMG hierarchy.
% Use a symmetric GaussSeidel smoother.
%
% The resolution is done twice. The first time with a standalone
% AMG method and the second time using AMG as a preconditioner to
% conjugate gradients.
%
% See also: Laplace2D BuildLaplace1D

srand;
clear all;
mue_include

SolStatus = NOTALLZEROS;

% Build problem matrix
%n      = 9;
%Amat   = BuildLaplace1D(n);
%
%n      = 10;
%Amat   = BuildLaplace2D(n);
%
%Amat = BuildLaplace1DBlk(4,0,20);
%n=size(Amat,1)
n=50;
[A,nullspace]=BuildElasticity2D(n,1e9,.25);
Amat = Operator(A.GetMatrixData(),2,2);


RowMap = Amat.GetRowMap();
ndofs = RowMap.NDOFs();

%nullspace = BuildNullSpace(Amat);
%nullspace = [nullspace (1:ndofs)']; % for tests
%nullspace = [nullspace nullspace];  % for tests
%nullspace = [(1:ndofs)' (ndofs+1:2*ndofs)']; % for tests
%nullspace = [repmat([1;0],ndofs/2,1),repmat([0;1],ndofs/2,1)];
%nullspace = [nullspace (1:ndofs)',[1:ndofs].^2']; % for tests

% Set options
numDesiredLevels = 2;         % number of AMG levels
options.NoQR=1;
%options.PtentRootModifications='4c';

% Options for Ray's PtentExperiment
%options.PtentSA=1;
options.NCoarseDofPerNode=2;
options.CoarseNullspaceWeightingMultiplier=[sqrt(10000),1,1];

eminSteps = 2;

AmalgamateDropFact = CoalesceDropFactory();
AggFact            = AggregationFactory();
CNSFact            = CoarseNSFactory();
Pinitfact          = TentativePFactoryEx(AmalgamateDropFact, AggFact, CNSFact, options);
PatFact            = AP_PatternFactory();
PatFact.SetPatternType('AP');
FCSplittingFact = FCSplittingFactory(AggFact);
ConstraintFact     = ConstraintFactory(PatFact, FCSplittingFact, options);
Pfact              = EminPFactory(PatFact, ConstraintFact, CGEminSolver(eminSteps), Pinitfact, options);

Rfact              = TransPFactory();
Acfact             = RAPFactory();
GSFactory          = SmootherFactory(Smoother('GaussSeidel', 2, 1));


%
%  Construct and populate finest level with user information
%
Finest = Level();
Finest.Set('A', Amat);
Finest.Set('NullSpace', nullspace);

MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(1);
MgHierarchy.SetLevel(Finest,1);

MgHierarchy.FillHierarchy(Pfact, Rfact, Acfact, 1, numDesiredLevels);
MgHierarchy.SetSmoothers(GSFactory);

rhs = rand(RowMap.NDOFs(),1);

% Maximum number of iterations
maxit=9;

% Tolerance parameter (used only in preconditioned conjugate gradients method)
% tol = [];  % default = 1e-6
tol = eps*2; % near machine precision

% AMG without CG
fprintf('AMG without CG\n');
sol = zeros(RowMap.NDOFs(),1); SolStatus = ALLZEROS;
sol = MgHierarchy.Iterate(rhs, maxit, sol, SolStatus); SolStatus = NOTALLZEROS;
norm(rhs-Amat.GetMatrixData()*sol)/norm(rhs)

% AMG as a preconditioner to CG
fprintf('CG preconditioned by AMG\n');
sol   = zeros(RowMap.NDOFs(),1); SolStatus = ALLZEROS;
zeroGuess = zeros(RowMap.NDOFs(),1);
[sol,flag,relres,iter,resvec] = pcg(Amat.GetMatrixData(),rhs,tol,maxit,...
             @(rhs)MgHierarchy.Iterate(rhs,1, zeroGuess,ALLZEROS,VCYCLE));
SolStatus = NOTALLZEROS;
%for ii=2:length(resvec),
%  fprintf('%3d-: ||r||=%e\n',ii-1,resvec(ii));
%end
norm(rhs-Amat.GetMatrixData()*sol)/norm(rhs)
