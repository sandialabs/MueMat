clear all;

%% define problem
SetHomeDir();

if exist('result.mat') == 2
    load('results.mat');
    clear all
    save('results.mat');
else
    save('results.mat');
end


if ~exist(sprintf('vg_32_32.mat'),'file')
    sc = Scatra2D_vg1(31,31); 
    [A,b,problem] = sc.Build();
    n     = A.GetRowMap().NDOFs();
    guess = zeros(n,1);
    numDesiredLevels = 3;
    mue_include 
    save (sprintf('vg_32_32.mat'));
else
    load (sprintf('vg_32_32.mat'));
end

%% PA-AMG
clear all
load (sprintf('vg_32_32.mat'));

% setup multigrid
AggFact   = AggregationFactory();
Pfact = TentativePFactory(CoalesceDropFactory(), AggFact)
Rfact     = TransPFactory();
Acfact    = RAPFactory();

Smofact   = SmootherFactory(Smoother('Jacobi',1,0.7));

% prepare multigrid level
Finest = Level();
Finest.Set('A',A);
Finest.Set('NullSpace', BuildNullSpace(A));

MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(0);
MgHierarchy.SetLevel(Finest,1);
MgHierarchy.SetMaxCoarseSize(1);
MgHierarchy.FillHierarchy(Pfact, Rfact, Acfact, 1, numDesiredLevels);
MgHierarchy.SetSmoothers(Smofact);

if (MgHierarchy.GetNumLevel() ~= numDesiredLevels)
error('MgHierarchy.GetNumLevel() ~= numDesiredLevels: change the problem size for this test');
end

% AMG as a preconditioner to CG
maxIts = 100;     % Maximum number of iterations
tol    = 1e-8; % Tolerance parameter

params(1) = tol; params(2) = maxIts;
%mysol = MgHierarchy.Iterate(b,10,zeros(length(b),1),ALLZEROS,1);
[mysol,flag,iter,resvec]= mygmres(A, @(v)MgHierarchy.Iterate(v,1,zeros(length(b),1),ALLZEROS,1), b,params,zeros(length(b),1));
fprintf('Nits = %d,||r||=%e\n',length(resvec)-1,norm(b-A.GetMatrixData()*mysol)/norm(b));

% plot results
subplot(2,2,1);
sc.plot_vector(mysol);

subplot(2,2,2);
PP = MgHierarchy.Levels_{2}.Get('P').GetMatrixData();
surf(reshape(-PP(:,57),32,32));
subplot(2,2,4);
RR = MgHierarchy.Levels_{2}.Get('R').GetMatrixData()';
surf(reshape(-RR(:,57),32,32));

subplot(2,2,3);
AggInfo = MgHierarchy.Levels_{1}.Get('Aggregates');
fpoints = ones(PP,1); fpoints(AggInfo.Roots) = 0; fpoints = find(fpoints);
fdofs = Node2DOF(fpoints,A.GetRowMap());
rdofs = Node2DOF(AggInfo.Roots,A.GetRowMap());
AA = A.GetMatrixData();
Acc = AA(rdofs,rdofs);
Aff = AA(fdofs,fdofs);
Afc = AA(fdofs,rdofs);
Acf = AA(rdofs,fdofs);
SS  = Acc - Acf*inv(Aff)*Afc;
A2  = MgHierarchy.Levels_{2}.Get('A').GetMatrixData();
surf(SS-A2);
pause

%% setup multigrid
AggFact   = AggregationFactory();
Pinitfact = TentativePFactory(CoalesceDropFactory(), AggFact)
Pfact     = OptimalPFactory(AggFact); %PgPFactory(Pinitfact);
Rfact     = GenericRFactory(Pfact);
PRfact    = GenericPRFactory(Pfact,Rfact);
PRfact.SetMaxCoarseSize(1);
Acfact    = RAPFactory();

Smofact   = SmootherFactory(Smoother('Jacobi',1,0.7));

% prepare multigrid level
Finest = Level();
Finest.Set('A',A);
Finest.Set('NullSpace', BuildNullSpace(A));

MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(0);
MgHierarchy.SetLevel(Finest,1);

MgHierarchy.FillHierarchy(PRfact,[], Acfact, 1, numDesiredLevels);
MgHierarchy.SetSmoothers(Smofact);

if (MgHierarchy.GetNumLevel() ~= numDesiredLevels)
error('MgHierarchy.GetNumLevel() ~= numDesiredLevels: change the problem size for this test');
end

% AMG as a preconditioner to CG
maxIts = 100;     % Maximum number of iterations
tol    = 1e-8; % Tolerance parameter

params(1) = tol; params(2) = maxIts;
%mysol = MgHierarchy.Iterate(b,10,zeros(length(b),1),ALLZEROS,1);
[mysol,flag,iter,resvec]= mygmres(A, @(v)MgHierarchy.Iterate(v,1,zeros(length(b),1),ALLZEROS,1), b,params,zeros(length(b),1));
fprintf('Nits = %d,||r||=%e\n',length(resvec)-1,norm(b-A.GetMatrixData()*mysol)/norm(b));

% plot results
subplot(2,2,1);
sc.plot_vector(mysol);

subplot(2,2,2);
PP = MgHierarchy.Levels_{2}.Get('P').GetMatrixData();
surf(reshape(PP(:,57),32,32));
subplot(2,2,4);
RR = MgHierarchy.Levels_{2}.Get('R').GetMatrixData()';
surf(reshape(RR(:,57),32,32));

subplot(2,2,3);
AggInfo = MgHierarchy.Levels_{1}.Get('Aggregates');
fpoints = ones(PP,1); fpoints(AggInfo.Roots) = 0; fpoints = find(fpoints);
fdofs = Node2DOF(fpoints,A.GetRowMap());
rdofs = Node2DOF(AggInfo.Roots,A.GetRowMap());
AA = A.GetMatrixData();
Acc = AA(rdofs,rdofs);
Aff = AA(fdofs,fdofs);
Afc = AA(fdofs,rdofs);
Acf = AA(rdofs,fdofs);
SS  = Acc - Acf*inv(Aff)*Afc;
A2  = MgHierarchy.Levels_{2}.Get('A').GetMatrixData();
surf(SS-A2);

pause

%% PG-AMG
clear all
load (sprintf('vg_32_32.mat'));

% setup multigrid
AggFact   = AggregationFactory();
Pinitfact = TentativePFactory(CoalesceDropFactory(), AggFact)
Pfact     = PgPFactory(Pinitfact);
Rfact     = GenericRFactory(Pfact);
Acfact    = RAPFactory();

Smofact   = SmootherFactory(Smoother('Jacobi',1,0.7));

% prepare multigrid level
Finest = Level();
Finest.Set('A',A);
Finest.Set('NullSpace', BuildNullSpace(A));

MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(0);
MgHierarchy.SetLevel(Finest,1);
MgHierarchy.SetMaxCoarseSize(1);
MgHierarchy.FillHierarchy(Pfact,Rfact, Acfact, 1, numDesiredLevels);
MgHierarchy.SetSmoothers(Smofact);

if (MgHierarchy.GetNumLevel() ~= numDesiredLevels)
error('MgHierarchy.GetNumLevel() ~= numDesiredLevels: change the problem size for this test');
end

% AMG as a preconditioner to CG
maxIts = 100;     % Maximum number of iterations
tol    = 1e-8; % Tolerance parameter

params(1) = tol; params(2) = maxIts;
%mysol = MgHierarchy.Iterate(b,10,zeros(length(b),1),ALLZEROS,1);
[mysol,flag,iter,resvec]= mygmres(A, @(v)MgHierarchy.Iterate(v,1,zeros(length(b),1),ALLZEROS,1), b,params,zeros(length(b),1));
fprintf('Nits = %d,||r||=%e\n',length(resvec)-1,norm(b-A.GetMatrixData()*mysol)/norm(b));

% plot results
subplot(2,2,1);
sc.plot_vector(mysol);

subplot(2,2,2);
PP = MgHierarchy.Levels_{2}.Get('P').GetMatrixData();
surf(reshape(-PP(:,57),32,32));
subplot(2,2,4);
RR = MgHierarchy.Levels_{2}.Get('R').GetMatrixData()';
surf(reshape(-RR(:,57),32,32));

subplot(2,2,3);
AggInfo = MgHierarchy.Levels_{1}.Get('Aggregates');
fpoints = ones(PP,1); fpoints(AggInfo.Roots) = 0; fpoints = find(fpoints);
fdofs = Node2DOF(fpoints,A.GetRowMap());
rdofs = Node2DOF(AggInfo.Roots,A.GetRowMap());
AA = A.GetMatrixData();
Acc = AA(rdofs,rdofs);
Aff = AA(fdofs,fdofs);
Afc = AA(fdofs,rdofs);
Acf = AA(rdofs,fdofs);
SS  = Acc - Acf*inv(Aff)*Afc;
A2  = MgHierarchy.Levels_{2}.Get('A').GetMatrixData();
surf(SS-A2);
pause

%% MinDescent
clear all
load (sprintf('vg_32_32.mat'));

% setup multigrid
AggFact   = AggregationFactory();
Pinitfact = TentativePFactory(CoalesceDropFactory(), AggFact)
Pfact     = EminPFactory(AP_PatternFactory, ConstraintFactory(), MinDescentSolver(8,1.0,'global'), Pinitfact);
Rfact     = GenericRFactory(Pfact);
PRfact    = GenericPRFactory(Pfact,Rfact);
PRfact.SetMaxCoarseSize(1);
Acfact    = RAPFactory();

Smofact   = SmootherFactory(Smoother('Jacobi',1,0.7));

% prepare multigrid level
Finest = Level();
Finest.Set('A',A);
Finest.Set('NullSpace', BuildNullSpace(A));

MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(0);
MgHierarchy.SetLevel(Finest,1);

MgHierarchy.FillHierarchy(PRfact,[], Acfact, 1, numDesiredLevels);
MgHierarchy.SetSmoothers(Smofact);

if (MgHierarchy.GetNumLevel() ~= numDesiredLevels)
error('MgHierarchy.GetNumLevel() ~= numDesiredLevels: change the problem size for this test');
end

% AMG as a preconditioner to CG
maxIts = 100;     % Maximum number of iterations
tol    = 1e-8; % Tolerance parameter

params(1) = tol; params(2) = maxIts;
%mysol = MgHierarchy.Iterate(b,10,zeros(length(b),1),ALLZEROS,1);
[mysol,flag,iter,resvec]= mygmres(A, @(v)MgHierarchy.Iterate(v,1,zeros(length(b),1),ALLZEROS,1), b,params,zeros(length(b),1));
fprintf('Nits = %d,||r||=%e\n',length(resvec)-1,norm(b-A.GetMatrixData()*mysol)/norm(b));

% plot results
subplot(2,2,1);
sc.plot_vector(mysol);

subplot(2,2,2);
PP = MgHierarchy.Levels_{2}.Get('P').GetMatrixData();
surf(reshape(-PP(:,57),32,32));
subplot(2,2,4);
RR = MgHierarchy.Levels_{2}.Get('R').GetMatrixData()';
surf(reshape(-RR(:,57),32,32));

subplot(2,2,3);
AggInfo = MgHierarchy.Levels_{1}.Get('Aggregates');
fpoints = ones(PP,1); fpoints(AggInfo.Roots) = 0; fpoints = find(fpoints);
fdofs = Node2DOF(fpoints,A.GetRowMap());
rdofs = Node2DOF(AggInfo.Roots,A.GetRowMap());
AA = A.GetMatrixData();
Acc = AA(rdofs,rdofs);
Aff = AA(fdofs,fdofs);
Afc = AA(fdofs,rdofs);
Acf = AA(rdofs,fdofs);
SS  = Acc - Acf*inv(Aff)*Afc;
A2  = MgHierarchy.Levels_{2}.Get('A').GetMatrixData();
surf(SS-A2);
pause

%% MinDescent (AffInv)
clear all
load (sprintf('vg_32_32.mat'));

% setup multigrid
AggFact   = AggregationFactory();
Pinitfact = TentativePFactory(CoalesceDropFactory(), AggFact)
Filter    = InlineNNZperRowFilter('0.9*LevelID*LevelID');
PatFact   = AffInvAfc_PatternFactory(Filter,Pinitfact,[]);
Pfact     = EminPFactory(PatFact, ConstraintFactory(), MinDescentSolver(8,1.0,'global'), Pinitfact);
Rfact     = GenericRFactory(Pfact);
PRfact    = GenericPRFactory(Pfact,Rfact);
PRfact.SetMaxCoarseSize(1);
Acfact    = RAPFactory();

Smofact   = SmootherFactory(Smoother('Jacobi',1,0.7));

% prepare multigrid level
Finest = Level();
Finest.Set('A',A);
Finest.Set('NullSpace', BuildNullSpace(A));

MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(0);
MgHierarchy.SetLevel(Finest,1);

MgHierarchy.FillHierarchy(PRfact,[], Acfact, 1, numDesiredLevels);
MgHierarchy.SetSmoothers(Smofact);

if (MgHierarchy.GetNumLevel() ~= numDesiredLevels)
error('MgHierarchy.GetNumLevel() ~= numDesiredLevels: change the problem size for this test');
end

% AMG as a preconditioner to CG
maxIts = 100;     % Maximum number of iterations
tol    = 1e-8; % Tolerance parameter

params(1) = tol; params(2) = maxIts;
%mysol = MgHierarchy.Iterate(b,10,zeros(length(b),1),ALLZEROS,1);
[mysol,flag,iter,resvec]= mygmres(A, @(v)MgHierarchy.Iterate(v,1,zeros(length(b),1),ALLZEROS,1), b,params,zeros(length(b),1));
fprintf('Nits = %d,||r||=%e\n',length(resvec)-1,norm(b-A.GetMatrixData()*mysol)/norm(b));

% plot results
subplot(2,2,1);
sc.plot_vector(mysol);

subplot(2,2,2);
PP = MgHierarchy.Levels_{2}.Get('P').GetMatrixData();
surf(reshape(-PP(:,57),32,32));
subplot(2,2,4);
RR = MgHierarchy.Levels_{2}.Get('R').GetMatrixData()';
surf(reshape(-RR(:,57),32,32));

subplot(2,2,3);
AggInfo = MgHierarchy.Levels_{1}.Get('Aggregates');
fpoints = ones(PP,1); fpoints(AggInfo.Roots) = 0; fpoints = find(fpoints);
fdofs = Node2DOF(fpoints,A.GetRowMap());
rdofs = Node2DOF(AggInfo.Roots,A.GetRowMap());
AA = A.GetMatrixData();
Acc = AA(rdofs,rdofs);
Aff = AA(fdofs,fdofs);
Afc = AA(fdofs,rdofs);
Acf = AA(rdofs,fdofs);
SS  = Acc - Acf*inv(Aff)*Afc;
A2  = MgHierarchy.Levels_{2}.Get('A').GetMatrixData();
surf(SS-A2);
pause

%% MinDescent (AffInv)
clear all
load (sprintf('vg_32_32.mat'));

% setup multigrid
AggFact   = AggregationFactory();
Pinitfact = TentativePFactory(CoalesceDropFactory(), AggFact)
Filter    = InlineNNZperRowFilter('5*LevelID*LevelID');
PatFact   = AffInvAfc_PatternFactory(Filter,Pinitfact,[]);
Pfact     = EminPFactory(PatFact, ConstraintFactory(), MinDescentSolver(4,0.1,'adaptive'), Pinitfact);
Rfact     = GenericRFactory(Pfact);
PRfact    = GenericPRFactory(Pfact,Rfact);
PRfact.SetMaxCoarseSize(1);
Acfact    = RAPFactory();

Smofact   = SmootherFactory(Smoother('Jacobi',1,0.7));

% prepare multigrid level
Finest = Level();
Finest.Set('A',A);
Finest.Set('NullSpace', BuildNullSpace(A));

MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(0);
MgHierarchy.SetLevel(Finest,1);

MgHierarchy.FillHierarchy(PRfact,[], Acfact, 1, numDesiredLevels);
MgHierarchy.SetSmoothers(Smofact);

if (MgHierarchy.GetNumLevel() ~= numDesiredLevels)
error('MgHierarchy.GetNumLevel() ~= numDesiredLevels: change the problem size for this test');
end

% AMG as a preconditioner to CG
maxIts = 100;     % Maximum number of iterations
tol    = 1e-8; % Tolerance parameter

params(1) = tol; params(2) = maxIts;
%mysol = MgHierarchy.Iterate(b,10,zeros(length(b),1),ALLZEROS,1);
[mysol,flag,iter,resvec]= mygmres(A, @(v)MgHierarchy.Iterate(v,1,zeros(length(b),1),ALLZEROS,1), b,params,zeros(length(b),1));
fprintf('Nits = %d,||r||=%e\n',length(resvec)-1,norm(b-A.GetMatrixData()*mysol)/norm(b));

% plot results
subplot(2,2,1);
sc.plot_vector(mysol);

subplot(2,2,2);
PP = MgHierarchy.Levels_{2}.Get('P').GetMatrixData();
surf(reshape(-PP(:,57),32,32));
subplot(2,2,4);
RR = MgHierarchy.Levels_{2}.Get('R').GetMatrixData()';
surf(reshape(-RR(:,57),32,32));

subplot(2,2,3);
AggInfo = MgHierarchy.Levels_{1}.Get('Aggregates');
fpoints = ones(PP,1); fpoints(AggInfo.Roots) = 0; fpoints = find(fpoints);
fdofs = Node2DOF(fpoints,A.GetRowMap());
rdofs = Node2DOF(AggInfo.Roots,A.GetRowMap());
AA = A.GetMatrixData();
Acc = AA(rdofs,rdofs);
Aff = AA(fdofs,fdofs);
Afc = AA(fdofs,rdofs);
Acf = AA(rdofs,fdofs);
SS  = Acc - Acf*inv(Aff)*Afc;
A2  = MgHierarchy.Levels_{2}.Get('A').GetMatrixData();
surf(SS-A2);

%% Emin (GMRES)
clear all
load (sprintf('vg_32_32.mat'));

% setup multigrid
AggFact   = AggregationFactory();
Pinitfact = TentativePFactory(CoalesceDropFactory(), AggFact)
Pfact     = EminPFactory(AP_PatternFactory, ConstraintFactory(), GMRESEminSolver(100), Pinitfact);
Rfact     = GenericRFactory(Pfact);
PRfact    = GenericPRFactory(Pfact,Rfact);
PRfact.SetMaxCoarseSize(1);
Acfact    = RAPFactory();

Smofact   = SmootherFactory(Smoother('Jacobi',1,0.7));

% prepare multigrid level
Finest = Level();
Finest.Set('A',A);
Finest.Set('NullSpace', BuildNullSpace(A));

MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(0);
MgHierarchy.SetLevel(Finest,1);

MgHierarchy.FillHierarchy(PRfact,[], Acfact, 1, numDesiredLevels);
MgHierarchy.SetSmoothers(Smofact);

if (MgHierarchy.GetNumLevel() ~= numDesiredLevels)
error('MgHierarchy.GetNumLevel() ~= numDesiredLevels: change the problem size for this test');
end

% AMG as a preconditioner to CG
maxIts = 100;     % Maximum number of iterations
tol    = 1e-8; % Tolerance parameter

params(1) = tol; params(2) = maxIts;
%mysol = MgHierarchy.Iterate(b,10,zeros(length(b),1),ALLZEROS,1);
[mysol,flag,iter,resvec]= mygmres(A, @(v)MgHierarchy.Iterate(v,1,zeros(length(b),1),ALLZEROS,1), b,params,zeros(length(b),1));
fprintf('Nits = %d,||r||=%e\n',length(resvec)-1,norm(b-A.GetMatrixData()*mysol)/norm(b));

% plot results
subplot(2,2,1);
sc.plot_vector(mysol);

subplot(2,2,2);
PP = MgHierarchy.Levels_{2}.Get('P').GetMatrixData();
surf(reshape(-PP(:,57),32,32));
subplot(2,2,4);
RR = MgHierarchy.Levels_{2}.Get('R').GetMatrixData()';
surf(reshape(-RR(:,57),32,32));

subplot(2,2,3);
AggInfo = MgHierarchy.Levels_{1}.Get('Aggregates');
fpoints = ones(PP,1); fpoints(AggInfo.Roots) = 0; fpoints = find(fpoints);
fdofs = Node2DOF(fpoints,A.GetRowMap());
rdofs = Node2DOF(AggInfo.Roots,A.GetRowMap());
AA = A.GetMatrixData();
Acc = AA(rdofs,rdofs);
Aff = AA(fdofs,fdofs);
Afc = AA(fdofs,rdofs);
Acf = AA(rdofs,fdofs);
SS  = Acc - Acf*inv(Aff)*Afc;
A2  = MgHierarchy.Levels_{2}.Get('A').GetMatrixData();
surf(SS-A2);