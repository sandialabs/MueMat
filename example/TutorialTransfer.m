%% A tutorial that demonstrates the usage of different transfer operators
% This tutorial presents the basic usage of different transfer operator
% strategies (e.g. SA-AMG, PG-AMG...).

%%
% First we reinitialize the MATLAB workspace and include some basic MueMat
% variables.
srand;
clear all;
mue_include;
ITERS = zeros(3,1);
OC    = zeros(3,1);

%% define a linear problem
% We're using a nonsymmetric linear
% system from a finite element discretization of a scalar convection
% diffusion equation
% load problem data
SetHomeDir
mytests = { [MUEMAT_ROOT_DIR '/data/TutorialTransfer.mat']};
load(mytests{1});
n     = Amat.GetRowMap().NDOFs();
guess = zeros(n,1);

%%
% Now we're solving the linear system Amat x = rhs for x using the MueMat
% AMG solver with different transfer operators.

%% Basic AMG setup
% # allocate a |Level| object for the finest level and associate the fine
%   level matrix with that data bucket.
% # setup the AMG smoothers. In this example we only use three sweeps with
%   a damped Jacobi iteration on each multigrid level
% # prepare multgrid hierarchy
Finest=Level();
Finest.Set('A', Amat);

Sfact = SmootherFactory(Smoother('Jacobi',3,0.5));

MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(10);
MgHierarchy.SetLevel(Finest,1);

%% Setup of multigrid transfer operators
% The |MgHierarchy.FillHierarchy()| function needs a |PRFactory| derived
% object that handles the transfer operators between the multigrid levels.
% In MueMat you can easily combine different prolongation and
% restriction operators (that derive form |PFactory| and |RFactory|
% respectively) and put them together using the |GenericPRFactory| object.


%% Example 1: Plain-Aggregation (PA-AMG)
%%
% generate 3 multigrid levels using plain aggregation and set smoothers
status = MgHierarchy.FillHierarchy(TentativePFactory(), TransPFactory(), RAPFactory(), 1, 3);
MgHierarchy.SetSmoothers(Sfact);

%%
% Invoke the multilevel scheme either as solver or preconditioner.
[sol,flag,relres,iter,resvec] = gmres(Amat.GetMatrixData(),rhs,[],1e-7,n,@(rhs)MgHierarchy.Iterate(rhs,1,zeros(n,1),ALLZEROS),[],guess);
ITERS(1) = iter(2);
OC(1)    = status.OperatorComplexity;

%% Example 2: Smoothed aggregation (SA-AMG)
% The same as before. Now we're applying the smoothed aggregation
% prolongator and use the transposed of the prolongator for the restriction
PRfact = GenericPRFactory(SaPFactory(), TransPFactory());

%%
% again fill multigrid hierarchy using the new transfer operators
status = MgHierarchy.FillHierarchy(PRfact, RAPFactory(), 1, 3);
MgHierarchy.SetSmoothers(Sfact);
[sol,flag,relres,iter,resvec] = gmres(Amat.GetMatrixData(),rhs,[],1e-7,n,@(rhs)MgHierarchy.Iterate(rhs,1,zeros(n,1),ALLZEROS),[],guess);
ITERS(2) = iter(2);
OC(2)    = status.OperatorComplexity;

%% Example 3: PG-AMG
% For nonsymmetric problems it's not the best idea to use the transposed of
% the prolongation operator for the restriction. There are more advanced
% strategies for transfer operators, that are designed for nonsymmetric
% problems and itself provide separate prolongation and restriction
% operators.
% |PgPRFactory| implements both the methods for the prolongation and the
% restriction operator.
PRfact = PgPRFactory();
status = MgHierarchy.FillHierarchy(PRfact, RAPFactory(), 1, 3);
MgHierarchy.SetSmoothers(Sfact);
[sol,flag,relres,iter,resvec] = gmres(Amat.GetMatrixData(),rhs,[],1e-7,n,@(rhs)MgHierarchy.Iterate(rhs,1,zeros(n,1),ALLZEROS),[],guess);
ITERS(3) = iter(2);
OC(3)    = status.OperatorComplexity;

%% compare effect of different transfer operator strategies
fprintf('\n');
fprintf('        |  PA-AMG   |  SA-AMG    |   PG-AMG    \n');
fprintf(' SZ LVL | ITS   OC  | ITS   OC   |   ITS   OC  \n');
fprintf('-----------------------------------------------------\n');
fprintf('%3d %2d  | %2d  %4.2f  | %2d  %4.2f   |   %2d  %4.2f  \n',[20*20,3,ITERS(1),OC(1),ITERS(2),OC(2),ITERS(3),OC(3)]')
