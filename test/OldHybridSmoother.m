% a stupid hybrid smoother test

srand;
clear all;
global VERBOSE_CPY
VERBOSE_CPY=1;

mue_include

% Set options
numDesiredLevels = 2;         % number of AMG levels

AmalgamateDropFact = CoalesceDropFactory();
AggFact            = AggregationFactory();
%AggFact           = AggFact.SetAlgorithm('graph');
Ptentfact          = TentativePFactory(AmalgamateDropFact,AggFact);
Pfact              = SaPFactory(Ptentfact);
Rfact              = TransPFactory();
PRfact             = GenericPRFactory(Pfact,Rfact);
PRfact.SetMaxCoarseSize(30);
Acfact             = RAPFactory();
%SmooFactory          = SmootherFactory(Smoother('GaussSeidel', 2, 1));
%SmooFactory          = SmootherFactory(ChebySmoother(2,1/30));

NumR = 50;
nMainIts = 5;
MiddleOne = 5;
StartTwo = 5;
MiddleTwo = 5;
EndTwo = 5;

SmootherOne = Smoother('GaussSeidel', 2, 1);
SmootherTwo = Smoother('GaussSeidel', 2, 1);

Smoo = OldHybrid2x2Smoother(nMainIts, MiddleOne, ...
                      StartTwo, MiddleTwo, EndTwo, ...
                      SmootherOne, SmootherTwo);

SmooFactory = OldHybrid2x2SmootherFactory(Smoo);

%
%  Construct and populate finest level with user information
%
n      = 30; 
Amat   = BuildLaplace2D(n);
Finest = Level();
Finest.Set('A', Amat);
Finest.Set('NullSpace', BuildNullSpace(Amat));

Bsize = 1;
AmatData=Amat.GetMatrixData();
AOne = Operator(AmatData(1:NumR,1:NumR),Bsize,Bsize);
Finest.Set('Arr', AOne);

MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(1);
MgHierarchy.SetLevel(Finest,1);
status = MgHierarchy.FillHierarchy(PRfact, Acfact,1,numDesiredLevels);
MgHierarchy.SetSmoothers(SmooFactory);

RowMap = Amat.GetRowMap();
rhs = rand(RowMap.NDOFs(),1);
sol = zeros(RowMap.NDOFs(),1); SolStatus = ALLZEROS;
sol = MgHierarchy.Iterate(rhs, 9, sol, SolStatus,1,1); SolStatus = NOTALLZEROS;

tol = eps(2); % near machine precision
maxit=99;

% Multigrid as a preconditioner to CG
sol   = zeros(RowMap.NDOFs(),1); SolStatus = ALLZEROS;
zeros = zeros(RowMap.NDOFs(),1);
[sol,flag,relres,iter,resvec]= pcg(Amat.GetMatrixData(),rhs,tol,maxit,...
                   @(rhs)MgHierarchy.Iterate(rhs,1, zeros,ALLZEROS));
SolStatus = NOTALLZEROS;
