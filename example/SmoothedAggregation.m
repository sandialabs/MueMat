% A 2D Linear Elasticity problem solves via a 2-level AMG hierarchy.
%
% We use the knowledge of the near nullspace modes of the initial
% matrix to constuct the tentative prolongator of the Smoothed
% Aggregation AMG.
%
% See also: BuildElasticity2D

srand;
clear all;
mue_include
SolStatus = NOTALLZEROS;

n = 20;

E = 1e5;
nu = 0.3;

[Amat, NullSpace] = BuildElasticity2D(n, E, nu);
Amat =  Operator(Amat.GetMatrixData(),2,2);

% Set options
numDesiredLevels = 4;         % number of AMG levels

AmalgamateDropFact= CoalesceDropFactory();
AggFact        = AggregationFactory();
Ptentfact      = TentativePFactory(AmalgamateDropFact,AggFact);
Pfact          = SaPFactory(Ptentfact);
Rfact             = TransPFactory();
Acfact            = RAPFactory();
GSFactory         = SmootherFactory(Smoother('GaussSeidel', 2, 1));

%
%  Construct and populate finest level with user information
%
Finest = Level();
Finest.Set('A', Amat);
Finest.Set('NullSpace', NullSpace);

MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(1);
MgHierarchy.SetLevel(Finest,1);
MgHierarchy.FillHierarchy(Pfact, Rfact, Acfact, 1, numDesiredLevels);
MgHierarchy.SetSmoothers(GSFactory);

rhs = rand(Amat.GetRowMap().NDOFs(),1);
sol = zeros(Amat.GetRowMap().NDOFs(),1);               SolStatus = ALLZEROS;
sol = MgHierarchy.Iterate(rhs, 15, sol,SolStatus); SolStatus = NOTALLZEROS;

norm(rhs-Amat.GetMatrixData()*sol)/norm(rhs)
