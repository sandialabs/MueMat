% Example for developing an interface between level data and factories.
%
% Example 1 from document "Level Interface: Getting at Level Data in MueMat/MueLu".
% Smoothed aggregation with coalesced A for aggregation
% == Standard SA (requiring only defaults)

clear all
srand
mue_include

%% Set problem
n=50;
[A,nullspace]=BuildElasticity2D(n,1e9,.25); A = Operator(A.GetMatrixData(),2,2);
nDofs =  A.GetRowMap().NDOFs();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set options

% By default, aggregation uses a coalesced graph issued by CoalesceDropFactory.

Pfact = SaPFactory();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Construct and populate finest level with user information
Finest = Level();
Finest.Set('A', A);
Finest.Set('NullSpace', nullspace);

%% Construction of the MG hierarchy
numDesiredLevels = 3; % number of AMG levels
MgHierarchy = Hierarchy(Finest); 
MgHierarchy.SetOutputLevel(1);
MgHierarchy.FillHierarchy(GenericPRFactory(Pfact), RAPFactory(), 1, numDesiredLevels);
MgHierarchy.SetSmoothers(SmootherFactory(Smoother('GaussSeidel')));

%% AMG as a preconditioner to CG
rhs = rand(nDofs,1);

% CG parameters
maxit = 9;
cgTol = 1.0e-10; % tolerance parameter

fprintf('CG preconditioned by AMG\n');
sol   = zeros(nDofs,1);
zeroGuess = zeros(nDofs,1);
[sol,flag,relres,iter,resvec] = pcg(A.GetMatrixData(),rhs,cgTol,maxit,...
             @(rhs)MgHierarchy.Iterate(rhs,1, zeroGuess,ALLZEROS,VCYCLE));
for ii=2:length(resvec),
  fprintf('%3d: ||r||=%e\n',ii-1,resvec(ii));
end
norm(rhs-A.GetMatrixData()*sol)/norm(rhs)
