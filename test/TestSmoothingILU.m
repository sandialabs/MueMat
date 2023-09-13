%% Demonstration of different ways to construct smoothers.
%
% This example builds a 2D Laplacian and solves via a 2-level multigrid hierarchy.
% The example shows how ILU smoothers can be used.

% *Problem setup*
srand;
mue_include
InitGuessStatus = NOTALLZEROS;

Amat = BuildLaplace2D(20);

rowmap = Amat.GetRowMap();
rhs  = rand(rowmap.NDOFs(),1);
sol  = zeros(rowmap.NDOFs(),1); InitGuessStatus = ALLZEROS;

%%
% *Factories to define multigrid*
%
Pfact             = SaPFactory();
Rfact             = TransPFactory();
PRfact            = GenericPRFactory(Pfact,Rfact);
Acfact            = RAPFactory();
%%
% *Construct and populate finest level with user information*
%
Finest = Level(); Finest.KeepAll(false); 
Finest.Set('A', Amat);
%Finest.Set('NullSpace', ones(rowmap.NDOFs(),1));
Finest.Set('NullSpace', BuildNullSpace(Amat));

MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(1);
MgHierarchy.SetLevel(Finest,1);
MgHierarchy.FillHierarchy(PRfact, [], Acfact, 1, 2);

%%
% *Different ways to set and exercise smoothers*
%
sweeps = 1;               % number of smoothing sweeps (pre or post)
startLevel = 1;           % start indexing at level 1
numDesiredLevels = 2;     % max number of levels in hierarchy
iterations = 9;           % multigrid iterations
smoothCoarsest = false;   % if true, run smoother on coarsest problem
                          % (must be true if a one-level method and you
                          %  don't want a direct solve)

%%
% *Symmetric Gauss-Seidel*
fprintf('Running symmetric Gauss-Seidel V(1,1) with w = .99\n');
SFact   = SmootherFactory(Smoother('GaussSeidel', sweeps, .99));
NewMg = MgHierarchy.Copy();
NewMg.SetSmoothers(SFact);
if (smoothCoarsest), NewMg.SetCoarsestSolver(SFact); end;
newsol   = NewMg.Iterate(rhs, iterations, sol, InitGuessStatus);

%%
%%
% *ILU*
fprintf('Running ILU with default parameters (currently ILU(0))\n');
SFact  = SmootherFactory(ILUSmoother());
NewMg = MgHierarchy.Copy();
NewMg.SetSmoothers(SFact);
if (smoothCoarsest), NewMg.SetCoarsestSolver(SFact); end;
newsol   = NewMg.Iterate(rhs, iterations, sol, InitGuessStatus);

%%
% *ILU*
fprintf('Running ILUC with droptol= 0.1\n');
ILU.type = 'crout';
ILU.droptol = 0.1;
SFact  = SmootherFactory(ILUSmoother(ILU));
NewMg = MgHierarchy.Copy();
NewMg.SetSmoothers(SFact);
if (smoothCoarsest), NewMg.SetCoarsestSolver(SFact); end;
newsol   = NewMg.Iterate(rhs, iterations, sol, InitGuessStatus);

%%
% *ILU*
fprintf('Running ILUC with droptol= 0.01\n');
ILU.type = 'crout';
ILU.droptol = 0.01;
SFact  = SmootherFactory(ILUSmoother(ILU));
NewMg = MgHierarchy.Copy();
NewMg.SetSmoothers(SFact);
if (smoothCoarsest), NewMg.SetCoarsestSolver(SFact); end;
newsol   = NewMg.Iterate(rhs, iterations, sol, InitGuessStatus);

%%
% *DirectSolve Smoother*
fprintf('Running DirectSolveSmoother\n');
SFact = SmootherFactory(DirectSolveSmoother());
NewMg = MgHierarchy.Copy();
NewMg.SetSmoothers(SFact);
if (smoothCoarsest), NewMg.SetCoarsestSolver(SFact); end;
newsol   = NewMg.Iterate(rhs, iterations, sol, InitGuessStatus);
