% Test Cheby with block diag and point diag

clear classes;

%%
% *Problem setup*
srand;
mue_include
InitGuessStatus = NOTALLZEROS;

%% Set test case parameters
% Parameters for 304 Stainless Steel
n=40;
ELASTIC_MODULUS=193e9; 
POISSONS_RATIO=.305;
%[Amat,nullspace]=BuildElasticity2D(n, ELASTIC_MODULUS, POISSONS_RATIO);
%Amat =  Operator(Amat.GetMatrixData(),2,2); %TODO
 
Amat = BuildLaplace2D(n);
Amat = Operator(Amat.GetMatrixData(),2,2);

rowmap = Amat.GetRowMap();
rhs  = rand(rowmap.NDOFs(),1);
sol  = zeros(rowmap.NDOFs(),1); InitGuessStatus = ALLZEROS;

%%
% *Factories to define AMG*
%
Pfact             = SaPFactory();
Rfact             = TransPFactory();
Acfact            = RAPFactory();
%%
% *Construct and populate finest level with user information*
%
Finest = Level();
Finest.KeepAll(false);
Finest.Set('A', Amat);
%Finest.Set('NullSpace', ones(rowmap.NDOFs(),1));
Finest.Set('NullSpace', BuildNullSpace(Amat));
%Finest.Set('NullSpace', nullspace);
  
MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(1);
MgHierarchy.SetLevel(Finest,1);
MgHierarchy.FillHierarchy(Pfact, Rfact, Acfact, 1, 2);

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
fprintf('Running Cheby with default diag\n');
SFact   = SmootherFactory(ChebySmoother());
NewMg   = MgHierarchy.Copy();
NewMg.SetSmoothers(SFact);
if (smoothCoarsest), NewMg.SetCoarsestSolver(SFact); end;
%newsol  = NewMg.Iterate(rhs, iterations, sol, InitGuessStatus);

%% Solve Ax=b (AMG as a preconditioner of an iterative method)
maxIts = 30;
tol  = 1e-13;

[x,flag,relres,iter,resvec] = ...                                % pcg() is the Matlab conjugate gradients method.
    pcg(                           ... %
        Amat.GetMatrixData(), rhs,                   ... % parameters: * Matrix and right-hand side
        tol, maxIts,               ... %             * Conjugate Gradient parameters
        @(rhs)NewMg.Iterate(rhs, 1) ... %             * AMG Preconditioner
       );


%%
fprintf('Running Cheby with point diag\n');
mySmoother=ChebySmoother();
mySmoother.SetDiagonalView('point');
SFact   = SmootherFactory(mySmoother);
NewMg   = MgHierarchy.Copy();
NewMg.SetSmoothers(SFact);
if (smoothCoarsest), NewMg.SetCoarsestSolver(SFact); end;
%newsol  = NewMg.Iterate(rhs, iterations, sol, InitGuessStatus);

%% Solve Ax=b (AMG as a preconditioner of an iterative method)
maxIts = 30;
tol  = 1e-13;

[x,flag,relres,iter,resvec] = ...                                % pcg() is the Matlab conjugate gradients method.
    pcg(                           ... %
        Amat.GetMatrixData(), rhs,                   ... % parameters: * Matrix and right-hand side
        tol, maxIts,               ... %             * Conjugate Gradient parameters
        @(rhs)NewMg.Iterate(rhs, 1) ... %             * AMG Preconditioner
       );
