% Example for developing an interface between level data and factories.
%
% Example 4 from document "Level Interface: Getting at Level Data in MueMat/MueLu".
% Smoothed aggregation with 
% - Auxiliary graph for aggregation
% - Filtered A where pattern matches decoalesced version of auxiliary graph for prolongator smoother

clear all
for test=1:2

srand
mue_include

Finest = Level();

%% Set problem
n=50;
[A,nullspace, coords]=BuildElasticity2D(n,1e9,.25); A = Operator(A.GetMatrixData(),2,2);
Finest.Set('xcoords', coords(:,1));
Finest.Set('ycoords', coords(:,2));
nDofs =  A.GetRowMap().NDOFs();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set options

%% Auxiliary graph for aggregation

% Create an Auxiliary graph using distance laplacian.
% Note: This graph is already coalesced
tmpCoalesceDropfact = CoalesceDropFactory();
AmalgMat = tmpCoalesceDropfact.Build_(A);
AuxMat = BuildDistanceLaplacian(AmalgMat,coords);

% Save it
Finest.Set('AuxMatrix', AuxMat);

% How to project it: it is regenerate (rather than coarsen).
Finest.Set('AuxMatrixFunc', @BuildDistanceLaplacian); % Option used by RAPFactory.m    

% Use it
% I see two way to do it
if (test==1)
  % 1) Do not use any CoalesceDropFactory. Define directly which Graph must be used in Aggregation:
  CoalesceDropfact = [];
  Aggfact = AggregationFactory([]);
  Aggfact.SetGraphName('AuxMatrix');
else
  % 2) Use a Coalesce factory (allows to use a AuxMatrix with additional dropping etc.)
  CoalesceDropfact = CoalesceDropfactory();
  CoalesceDropfact.SetAName('AuxMatrix');
  Aggfact = AggregationFactory(CoalesceDropfact);
end

%% Filtered A where pattern matches decoalesced version of auxiliary graph for prolongator smoother
%   By default, SaPFactory use the matrix A for smoothing.
%   An option of SaPFactory allows to use another matrix (SetAForSmoothing).
%   See also: Ex2_SA_FiltredA.m
%
%   The difference with Ex2_SA_FiltredA is that we need to build the Afiltred ourself
%

% Build Afiltred for level 1
filteredA = BuildFiltredMatrix(A, AuxMat);
Finest.Set('myAfiltred', filteredA);

% Define how to project Afiltred
% *** TODO: it's hard coded in RAP factory right now...  Maybe we need
% *** a more global mechanism for projecting stuffs.

Pfact = SaPFactory(TentativePFactory(CoalesceDropfact,Aggfact));
Pfact.SetAForSmoothing('myAfiltred');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Construct and populate finest level with user information
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
disp(norm(rhs-A.GetMatrixData()*sol)/norm(rhs));


end % for test=1:2
