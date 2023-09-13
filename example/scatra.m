
% small example script that uses the PG-AMG transfer operators for solving
% a convection dominated nonsymmetric system of equations

srand;
clear all;
mue_include;

%% prepare example, load data from file (data generated with scatra2d code)
load('dgl_48x48_1x100.mat'); % load example matrix and rhs
Amat = Operator(A,Map(size(A,1),1),Map(size(A,2),1),@MatlabApply, 'Scatra2D doubleglazing');
rhs = b;
clear A;clear K;clear C;clear D;clear F;clear b;clear xy;

%% prepare finest Level
Finest      = Level();                          % Allocate a smoothed aggregation data bucket.
Finest.Set('A', Amat);                                     % Associate the fine level matrix with the data bucket.

%% prepare PG-AMG transfer operators
CoalesceDropFact    = CoalesceDropFactory();
AggFact             = AggregationFactory();
PRfact              = PgPRFactory(CoalesceDropFact,AggFact);

%% setup multigrid hierarchy
MgHierarchy = Hierarchy();                             % Allocate the AMG hierarchy
MgHierarchy.SetOutputLevel(10);                         % Verbose
MgHierarchy.SetLevel(Finest,1);                       % Associate the data bucket with the finest level.
MgHierarchy.FillHierarchy(PRFact,RAPFactory(),1,5);      % note that restrictor operator is set by PFact, too!!! -> design issue?? a separate restrictor factory doesn't make sense...
MgHierarchy.SetSmoothers();                            % Set the default smoothers.

n     = Amat.GetRowMap().NDOFs();
guess = zeros(n,1);
CGtol = 100.*eps(2);
CGits = 99;

%sol = MgHierarchy.Iterate(rhs, 50, guess, ALLZEROS);                         % Use AMG as a solver.

%[sol,flag,relres,iter,resvec]= pcg(Amat.GetMatrixData(),rhs,CGtol,CGits,...  % Use AMG as a preconditioner.
%                   @(rhs)MgHierarchy.Iterate(rhs,1, zeros(n,1), ALLZEROS));

%% solve with GMRES (preconditioned with AMG)
[sol,flag,relres,iter,resvec] = gmres(Amat.GetMatrixData(),rhs,size(Amat.GetMatrixData(),1),1e-7,size(Amat.GetMatrixData(),1),@(rhs)MgHierarchy.Iterate(rhs,1,zeros(n,1),ALLZEROS),[],guess);
