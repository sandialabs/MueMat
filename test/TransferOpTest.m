
srand;
clear all;
mue_include;

%% define problem
% sc = Scatra2D(47,1/800); % initialize doubleglazing problem with 47 x 47 elements (=48x48 nodes)
% [Amat,rhs,problem] = sc.Build(); % generate problem (matrix, rhs, information about problem definition)
% n     = Amat.GetRowMap().NDOFs();
% guess = zeros(n,1);


SetHomeDir
mytests = { [MUEMAT_ROOT_DIR '/data/TransferOpTest.mat']};
load(mytests{1});
Amat = Operator(Amat);
fprintf('\nTransferOpTests:%30s \n=====================\n','data/TransferOpTest.mat');


%% container
ITERS = zeros(4,1);
OC = zeros(4,1);

nlevels = 3;  % number of multigrid levels


%% MINDESCENT with simple thresholding

% setup aggregates
AggFact = AggregationFactory();
FCFact  = FCSplittingFactory(AggFact);
CoalesceFact = CoalesceDropFactory();

Pinitfact = TentativePFactory(CoalesceFact,AggFact);

% setup pattern with filter
Thsld = Thresholding(0.005,FCFact);  % 0.0034 = 151 0.003 = 147 0.0029 = 143 0.0028 = 141 0.005 = 128
PatFact = AffInvAfc_PatternFactory(Thsld,FCFact);
PatFact.SetAffInverseIterations(25);

% setup transfer operator factory 
Pfact = EminPFactory(PatFact, ConstraintFactory(), MinDescentSolver(8,0.1,'global'), Pinitfact);
Rfact = GenericRFactory(Pfact);
% PRFact = GenericPRFactory(Pfact,GenericRFactory(Pfact));
% PRFact.SetMaxCoarseSize(1);

% setup smoothers
SFact = SmootherFactory(Smoother('Jacobi',3,0.5));

Finest=Level();                          % Allocate a smoothed aggregation data bucket.
Finest.KeepAll(false);
Finest.Set('A', Amat);                         % Associate the fine level matrix with the data bucket.

% setup multigrid hierarchy
MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(10);
MgHierarchy.SetLevel(Finest,1);
MgHierarchy.SetMaxCoarseSize(1);
status = MgHierarchy.FillHierarchy(Pfact, Rfact, RAPFactory(), 1, nlevels);

% set smoothers
MgHierarchy.SetSmoothers(SFact);                            % Set the default smoothers.

% use AMG as a precondtioner within GMRES
[sol,flag,relres,iter,resvec] = gmres(Amat.GetMatrixData(),rhs,[],1e-7,n,@(rhs)MgHierarchy.Iterate(rhs,1,zeros(n,1),ALLZEROS),[],guess);

ITERS(1) = iter(2);
OC(1) = status.OperatorComplexity;

%% MINDESCENT with AP pattern
Finest=Level();                          % Allocate a smoothed aggregation data bucket.
Finest.KeepAll(false);
Finest.Set('A', Amat);                         % Associate the fine level matrix with the data bucket.

MgHierarchy2 = Hierarchy();
MgHierarchy2.SetOutputLevel(10);
MgHierarchy2.SetLevel(Finest,1);
MgHierarchy2.SetMaxCoarseSize(1);
Pinitfact = TentativePFactory(CoalesceFact,AggFact);
Thsld2 = Thresholding(0.0, FCFact);
PatFact = AP_PatternFactory(Thsld2);
PatFact.SetDegree(1);
Pfact = EminPFactory(PatFact, ConstraintFactory(), MinDescentSolver(8,0.1,'global'), Pinitfact);
% PRfact = GenericPRFactory(Pfact,GenericRFactory(Pfact));
% PRfact.SetMaxCoarseSize(1);
status = MgHierarchy2.FillHierarchy(Pfact, GenericRFactory(Pfact), RAPFactory(), 1, nlevels);
MgHierarchy2.SetSmoothers(SFact);                            % Set the default smoothers.
[sol,flag,relres,iter,resvec] = gmres(Amat.GetMatrixData(),rhs,[],1e-7,n,@(rhs)MgHierarchy2.Iterate(rhs,1,zeros(n,1),ALLZEROS),[],guess);
ITERS(2) = iter(2);
OC(2) = status.OperatorComplexity;

%% PG-AMG
Finest=Level();                          % Allocate a smoothed aggregation data bucket.
Finest.KeepAll(false);
Finest.Set('A', Amat);                         % Associate the fine level matrix with the data bucket.
MgHierarchy3 = Hierarchy();
MgHierarchy3.SetOutputLevel(10);
MgHierarchy3.SetLevel(Finest,1);

% setup transfer operator factory
Pfact = PgPFactory(Pinitfact);
Rfact = GenericRFactory(Pfact);
PRfact = GenericPRFactory(Pfact,Rfact);
status = MgHierarchy3.FillHierarchy(PRfact, [], RAPFactory(), 1, nlevels);
MgHierarchy3.SetSmoothers(SFact);                            % Set the default smoothers.
[sol,flag,relres,iter,resvec] = gmres(Amat.GetMatrixData(),rhs,[],1e-7,n,@(rhs)MgHierarchy3.Iterate(rhs,1,zeros(n,1),ALLZEROS),[],guess);

ITERS(3) = iter(2);
OC(3) = status.OperatorComplexity;

%% SA-AMG

% setup transfer operator factory
PRFact = GenericPRFactory(SaPFactory(),TransPFactory()); %PgPRFactory();
status = MgHierarchy.FillHierarchy(PRFact, [], RAPFactory(), 1, nlevels);
MgHierarchy.SetSmoothers(SFact);                            % Set the default smoothers.
[sol,flag,relres,iter,resvec] = gmres(Amat.GetMatrixData(),rhs,[],1e-7,n,@(rhs)MgHierarchy.Iterate(rhs,1,zeros(n,1),ALLZEROS),[],guess);

ITERS(4) = iter(2);
OC(4) = status.OperatorComplexity;

%% plot solution
sc.plot_vector(sol);


fprintf('\n');
fprintf('        |  MinDesc  | MinDesc (AP) |   PG-AMG    |   SA-AMG   \n');
fprintf(' SZ LVL | ITS   OC  | ITS   OC     |   ITS   OC  |  ITS    OC  \n');
fprintf('---------------------------------------------------------------\n');
fprintf('%3d %2d | %2d  %4.2f  | %2d  %4.2f     |   %2d  %4.2f  |  %2d  %4.2f \n',[48*48,nlevels,ITERS(1),OC(1),ITERS(2),OC(2),ITERS(3),OC(3),ITERS(4),OC(4)]')

if ITERS(1)~=66 || ITERS(2)~=82 || ITERS(3)~=93 || ITERS(4)~=129
   error('number of iterations are wrong. something seems to be broken'); 
end

% note 6/2/2011: PG-AMG needs 93 iterations
%                with fallback strategy for 3 wiggly prolongation operator
%                basis functions we only need 85 iterations -> experimental
