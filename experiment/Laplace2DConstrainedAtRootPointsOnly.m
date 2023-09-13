% The goal of this experiment is to demonstrate that a 2D Laplace problem constrained only at root point do not interpolate correctly constant vector.
%
% For this experiment, we use the energy minimization algorithm with:
% - rectangle aggregates
% - options.DropConstraintsMethod = 'fpoints' and Pfact.SetConstraintWeight(0)

srand;
clear all;
mue_include;

rootPointsOnly = true; % default = true, (switch to 'false' to compare with default constraints)
if (rootPointsOnly), fprintf('Constraints at root point only\n'); else fprintf('Constraints: default\n'); end

n = 36;
Amat = BuildLaplace2D(n);
nullspace = ones(n*n,1);

%% Aggregation
AggFact            = AggregationFactory();
AggFact.SetAlgorithm('rectangle');
AggPtsPerDim(1)     = 6; AggPtsPerDim(2) = 6;
AggFact.SetAggPtsPerDim(AggPtsPerDim);
junk = [n n];
AggFact.SetDomainPtsPerDim(junk);

%% P Pattern
PatFact = AP_PatternFactory(); PatFact.SetPatternType('AP');
PatFact.SetDegree(1);

%% PFactory
options.NoQR=1;
options.NCoarseDofPerNode=1;
options.CoarseNullspaceWeightingMultiplier=[1];

Pinitfact = TentativePFactory(CoalesceDropFactory(),AggFact);

if (rootPointsOnly), options.DropConstraintsMethod = 'fpoints'; end % Drop constraints at Fpoints

Constraintfact = ConstraintFactory([],[],options);
if (rootPointsOnly), Constraintfact.SetConstraintWeight(0); end     % Drop constraints at Fpoints

Pfact = EminPFactory(PatFact,Constraintfact,CGEminSolver(100), Pinitfact, options);

%%
Rfact              = TransPFactory();
PRfact             = GenericPRFactory(Pfact,Rfact);
Acfact             = RAPFactory();
GSFactory          = SmootherFactory(Smoother('GaussSeidel', 2, 1));

%% Construct and populate finest level with user information
numDesiredLevels = 2; % number of AMG levels

Finest = Level();
Finest.Set('A', Amat);
Finest.Set('NullSpace', nullspace);

MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(1);
MgHierarchy.SetLevel(Finest,1);
MgHierarchy.FillHierarchy(PRfact, Acfact, 1, numDesiredLevels);
MgHierarchy.SetSmoothers(GSFactory);

%% Solve (AMG as a preconditioner to CG)
ndofs = Amat.GetRowMap().NDOFs();

rhs = rand(ndofs,1);

maxit=9;     % Maximum number of iterations
tol = eps*2; % Tolerance parameter (near machine precision)

sol   = zeros(ndofs,1); SolStatus = ALLZEROS;
zeroGuess = zeros(ndofs,1);
[sol,flag,relres,iter,resvec] = pcg(Amat.GetMatrixData(),rhs,tol,maxit,...
             @(rhs)MgHierarchy.Iterate(rhs,1, zeroGuess,ALLZEROS,VCYCLE));
SolStatus = NOTALLZEROS;
for ii=2:length(resvec),
  fprintf('%3d: ||r||=%e\n',ii-1,resvec(ii));
end
norm(rhs-Amat.GetMatrixData()*sol)/norm(rhs)

%% Display P(:,15) in 3D
P = MgHierarchy.Levels_{2}.Get('P').GetMatrixData();
mesh(reshape(P(:,15),n,n));

pause(5);

%% Display P*(constant vector) in 3D
coarseOne = P*ones(n,1);
mesh(reshape(coarseOne/6,n,n)); % hard coded scaling

%% Misc
% pause;
% ptent = MgHierarchy.Levels_{2}.Get('P',Pinitfact).GetMatrixData();
% mesh(reshape(ptent(:,15),n,n));
%
% pause(5);
% coarseOne = ptent * ones(n,1);
% mesh(reshape(coarseOne/6,n,n));

%cnull = MgHierarchy.Levels_{2}.Get('NullSpace');
%plot(cnull)
