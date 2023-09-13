% clear all;

% Solve a Laplace 1D problem
n      = 600; 
A      = BuildLaplace1D(n);
b      = ones(n,1);
numDesiredLevels = 4;

mue_include

%% Set AMG options

AmalgamateDropFact= CoalesceDropFactory();
AmalgamateDropFact.SetPostDropSpecifications([],[]);
AmalgamateDropFact.SetPreDropSpecifications(@ExperimentalDrop,.7);
AggFact        = AggregationFactory();
AggFact.SetRemoveAggsThisSizeOrLess(1);

options.DropConstraintsMethod = 'fpoints';
Pinitfact          = TentativePFactoryEx(AmalgamateDropFact,AggFact); %TentativePFactory();
%Pfact              = EminPFactory2(DecoalescedGraph_PatternFactory(), ConstraintFactory([],[],options), CGEminSolver(), Pinitfact);
Pfact              = EminPFactory(DecoalescedGraph_PatternFactory(), ConstraintFactory([],[],options), GMRESEminSolver(1), Pinitfact);
Rfact              = TransPFactory();
PRfact             = GenericPRFactory(Pfact,Rfact);
PRfact.SetMaxCoarseSize(30);
Acfact             = RAPFactory();
SmooFactory        = SmootherFactory(Smoother('GaussSeidel', 2, 1));

% Create the multigrid hierarchy
% Construct and populate finest level with user information
Finest = Level(); Finest.KeepAll(false); 
Finest.Set('A', A);
Finest.Set('NullSpace', BuildNullSpace(A));

CoarseLevel = Level(); CoarseLevel.KeepAll(false);
% Note: Instead of creating a first coarselevel, users can set the Keep()
% options directly on the finest level.


MgHierarchy = Hierarchy(); MgHierarchy.SetOutputLevel(1);
MgHierarchy.SetLevel(Finest,1);
MgHierarchy.SetLevel(CoarseLevel,2);

status = MgHierarchy.FillHierarchy(PRfact, Acfact, 1, numDesiredLevels);
MgHierarchy.SetSmoothers(SmooFactory);

if (MgHierarchy.GetNumLevel() ~= numDesiredLevels)
error('MgHierarchy.GetNumLevel() ~= numDesiredLevels: change the problem size for this test');
end

% AMG as a preconditioner to CG
maxIts = 9;     % Maximum number of iterations
tol    = 1e-12; % Tolerance parameter

x    = zeros(n,1); SolStatus = ALLZEROS;
zero = zeros(n,1);
[x, flag, relRes, nIts, resVec] = pcg(A.GetMatrixData(),b,tol,maxIts,@(b)MgHierarchy.Iterate(b,1, zero,ALLZEROS)); SolStatus = NOTALLZEROS;

