% Example for developing an interface between level data and factories.
%
% Example 2 from document "Level Interface: Getting at Level Data in MueMat/MueLu".
% Smoothed aggregation with 
% - Filtered + Coalesced A for aggregation
% - Filtered A for prolongation smoother

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

%% Filtered + Coalesced A for aggregation
%   CoalesceDropFactory is responsible to build the coalesced graph used by the aggregation.
%   An option of CoalesceDropFactory allows to filter A before generating the graph (PreDropping option).

CoalesceDropfact = CoalesceDropFactory();
CoalesceDropfact.SetPreDropSpecifications(@SAdrop, 0.01);

%% Filtered A for prolongation smoother
%   By default, SaPFactory use the matrix A for smoothing.
%   An option of SaPFactory allows to use another matrix (SetAForSmoothing).
%
%  Additional notes:
%   - By default, the name of the matrix filtered by CoalesceDropFactory is 'Afiltred'
%     (it can be changed by using CoalesceDropfact.SetAFiltredName(''))
%   - If this matrix is used by another factory, it will be saved in
%     the Level (automatically, by using the reference count mechanism of Level)

Pfact = SaPFactory(TentativePFactory(CoalesceDropfact,AggregationFactory()));
Pfact.SetAForSmoothing('Afiltered');

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
