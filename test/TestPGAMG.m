
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
load(mytests{1}); Amat = Operator(Amat);
fprintf('\nTransferOpTests:%30s \n=====================\n','data/TransferOpTest.mat');


%% container
ITERS = zeros(3,1);
OC = zeros(3,1);

nlevels = 3;  % number of multigrid levels

%% prepare AMG

Finest=Level();                          % Allocate a smoothed aggregation data bucket.
Finest.Set('A', Amat);                         % Associate the fine level matrix with the data bucket.
Finest.KeepAll(false);   % note: PgPRFactory needs KeepAll(true), otherwise the nullspace is recalculated on the second level!!!


%% PG-AMG (old implementation, tentative prolongator as initial guess)

% setup smoothers
SFact = SmootherFactory(Smoother('Jacobi',6,0.5));

% setup multigrid hierarchy
MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(10);
MgHierarchy.SetLevel(Finest,1);

% setup transfer operator factory
Pfact = PgPFactory();
Rfact = GenericRFactory(Pfact);
PRfact = GenericPRFactory(Pfact,Rfact);

guess = zeros(n,1);
status = MgHierarchy.FillHierarchy(PRfact, [], RAPFactory(), 1, nlevels);
MgHierarchy.SetSmoothers(SFact);                            % Set the default smoothers.
[sol,flag,relres,iter,resvec] = gmres(Amat.GetMatrixData(),rhs,[],1e-7,n,@(rhs)MgHierarchy.Iterate(rhs,1,zeros(n,1),ALLZEROS),[],guess);

ITERS(1) = iter(2);
OC(1) = status.OperatorComplexity;

iter
resvec
status

%% new PG-AMG implementation with tentative prolongator as initial guess

% setup multigrid hierarchy
Finest=Level();                          % Allocate a smoothed aggregation data bucket.
Finest.Set('A', Amat);                         % Associate the fine level matrix with the data bucket.
Finest.KeepAll(false);

MgHierarchy2 = Hierarchy();
MgHierarchy2.SetOutputLevel(10);
MgHierarchy2.SetLevel(Finest,1);

Pinitfact = TentativePFactory();
Pfact = PgPFactory(Pinitfact);
%Pfact = PgPFactory();
Rfact = GenericRFactory(Pfact);
% PRfact2 = GenericPRFactory(Pfact,Rfact);

guess = zeros(n,1);
status = MgHierarchy2.FillHierarchy(Pfact, Rfact, RAPFactory(), 1, nlevels);
MgHierarchy2.SetSmoothers(SFact);                            % Set the default smoothers.
[sol,flag,relres,iter,resvec2] = gmres(Amat.GetMatrixData(),rhs,[],1e-7,n,@(rhs)MgHierarchy2.Iterate(rhs,1,zeros(n,1),ALLZEROS),[],guess);

ITERS(2) = iter(2);
OC(2) = status.OperatorComplexity;


iter
resvec2
status

% %% new PG-AMG implementation with SA-AMG prolongator as initial guess
% 
% % % setup multigrid hierarchy
% % Finest=Level();                          % Allocate a smoothed aggregation data bucket.
% % Finest.Set('A', Amat);                         % Associate the fine level matrix with the data bucket.
% % Finest.KeepAll(false);
% % 
% % MgHierarchy = Hierarchy();
% % MgHierarchy.SetOutputLevel(10);
% % MgHierarchy.SetLevel(Finest,1);
% % 
% % Pinitfact = SaPFactory();
% % Pfact = PgPFactory(Pinitfact);
% % Rfact = GenericRFactory(Pfact);
% % PRfact3 = GenericPRFactory(Pfact,Rfact);
% % 
% % status = MgHierarchy.FillHierarchy(PRfact3, RAPFactory(), 1, nlevels);
% % MgHierarchy.SetSmoothers(SFact);                            % Set the default smoothers.
% % [sol,flag,relres,iter,resvec3] = gmres(Amat.GetMatrixData(),rhs,[],1e-7,n,@(rhs)MgHierarchy.Iterate(rhs,1,zeros(n,1),ALLZEROS),[],guess);
% % 
% % iter 
% % resvec3
% % status
% % 
% % ITERS(3) = iter(2);
% % OC(3) = status.OperatorComplexity;
% 


%% print out results
fprintf('\n');
fprintf(' | PgPFact(Tent)      | PgPFact2(Tent)     |   PgPFact(SaPFact)  |\n');
fprintf(' | ITS       OC       | ITS         OC     |   ITS         OC    |\n');
fprintf(' +--------------------+--------------------+---------------------+\n');
fprintf(' | %2d      %4.4f     | %2d        %4.4f   |    %2d       %4.2f    |\n',[ITERS(1),OC(1),ITERS(2),OC(2),ITERS(3),OC(3)]')


if ITERS(1) ~= 54 || ITERS(2) ~= 54 || ... %ITERS(3) ~= 51 || ...
   OC(1) ~= 1.1848 || OC(2) ~= 1.1848  %|| (abs(OC(3)-1.360333333333333) >1e-4)
   error('something seems to be broken'); 
end

% note: 6/2/2011: the operator complexity of PG-AMG is 1.1848
%                 using some fallback stragegy for wiggly PG-AMG
%                 prolongator basis functions the operator complexity is
%                 slightly lower
