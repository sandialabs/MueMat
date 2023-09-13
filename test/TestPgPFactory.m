% clear all;

% Solve a Laplace 1D problem
n      = 600; 
A      = BuildLaplace1D(n);
b      = rand(n,1);
numDesiredLevels = 4;

mue_include

%% Set AMG options

Pinitfact          = TentativePFactory();
Pfact              = PgPFactory(Pinitfact);
%Rfact              = TransPFactory();
Rfact              = GenericRFactory(Pfact);
Acfact             = RAPFactory();
SmooFactory        = SmootherFactory(Smoother('GaussSeidel', 2, 1));
% FPSmoother = FPointSmoother(Smoother('Jacobi',5,0.1));
% SmooFactory = SmootherFactory([],FPSmoother);

%% RUN 1
% Create the multigrid hierarchy
% Construct and populate finest level with user information
Finest = Level(); Finest.KeepAll(false); 
Finest.Keep('Aggregates');
Finest.Set('A', A);
Finest.Set('NullSpace', BuildNullSpace(A));


MgHierarchy = Hierarchy(); MgHierarchy.SetOutputLevel(1);
MgHierarchy.SetLevel(Finest,1);
MgHierarchy.SetMaxCoarseSize(30);

MgHierarchy.FillHierarchy(Pfact,Rfact, Acfact, 1, numDesiredLevels);
MgHierarchy.SetSmoothers(SmooFactory);

if (MgHierarchy.GetNumLevel() ~= numDesiredLevels)
error('MgHierarchy.GetNumLevel() ~= numDesiredLevels: change the problem size for this test');
end

% AMG as a preconditioner to CG
maxIts = 100;     % Maximum number of iterations
tol    = 1e-8; % Tolerance parameter

x    = zeros(n,1); SolStatus = NOTALLZEROS; %SolStatus = ALLZEROS;
zero = zeros(n,1);
[x, flag, relRes, nIts, resVec] = pcg(A.GetMatrixData(),b,tol,maxIts,@(b)MgHierarchy.Iterate(b,1, zero,ALLZEROS)); 

