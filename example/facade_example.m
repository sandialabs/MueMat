%% Facade classes
% This minimalistic example solves a linear system $Ax=b$ using smoothed
% aggregation multigrid
% The usage of MueMat facade classes for an easy setup of the AMG preconditioner/solver
% is demonstrated
%

%% Header
clear all;
mue_include;

%% Setup the multigrid Solver
Amat = gallery('poisson',30); % The test matrix
Amat = Operator(Amat,1,1,@MatlabApply);

MgHierarchy = PA_AMG.Setup(Amat);

%% Generate a random right-hand side ('b')
% 'srand' initializes the random number generator to ensure reproducible results
srand;
n = size(Amat,1);
b = rand(n,1);

%% Solve Ax=b (AMG used directly as a solver)
% A fixed number of iterations is performed
nIts = 9;
fprintf('Solving Ax=b with multigrid...\n\n');
x    = MgHierarchy.Iterate(b, nIts);

fprintf('\n||r_0|| / ||r_final|| = %g\n\n',norm(b-Amat*x)/norm(b));

%% Solve Ax=b (AMG as a preconditioner of an iterative method)
maxIts = nIts;
tol  = 1e-8;

MgHierarchy.SetOutputLevel(0);
fprintf('Solving Ax=b again with conjugate gradients preconditioned by multigrid...\n\n');
[x,flag,relres,iter,resvec] = ...                                % pcg() is the Matlab conjugate gradients method.
    pcg(                           ... %
        Amat.GetMatrixData(), b,                   ... % parameters: * Matrix and right-hand side
        tol, maxIts,               ... %             * Conjugate Gradient parameters
        @(v)MgHierarchy.Iterate(v, 1) ... %             * AMG Preconditioner
       );
for ii=1:length(resvec),
  fprintf('  %d: ||r||=%g\n',ii,resvec(ii));
end

fprintf('\n||r_0|| / ||r_final|| = %g\n\n',norm(b-Amat*x)/norm(b));

fprintf('######################################################################\n');

MgHierarchy = Poisson_AMG.Setup(Amat); % AMG for Poisson problems
MgHierarchy.SetOutputLevel(10);

%% Solve Ax=b (AMG used directly as a solver)
% A fixed number of iterations is performed
nIts = 9;
fprintf('Solving Ax=b with multigrid...\n\n');
x    = MgHierarchy.Iterate(b, nIts);

fprintf('\n||r_0|| / ||r_final|| = %g\n\n',norm(b-Amat*x)/norm(b));

%% Solve Ax=b (AMG as a preconditioner of an iterative method)
maxIts = nIts;
tol  = 1e-8;

MgHierarchy.SetOutputLevel(0);
fprintf('Solving Ax=b again with conjugate gradients preconditioned by multigrid...\n\n');
[x,flag,relres,iter,resvec] = ...                                % pcg() is the Matlab conjugate gradients method.
    pcg(                           ... %
        Amat.GetMatrixData(), b,                   ... % parameters: * Matrix and right-hand side
        tol, maxIts,               ... %             * Conjugate Gradient parameters
        @(v)MgHierarchy.Iterate(v, 1) ... %             * AMG Preconditioner
       );
for ii=1:length(resvec),
  fprintf('  %d: ||r||=%g\n',ii,resvec(ii));
end

fprintf('\n||r_0|| / ||r_final|| = %g\n\n',norm(b-Amat*x)/norm(b));