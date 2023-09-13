%% A minimalist example of the MueMat interface
% This example solves a linear system $Ax=b$ using the default
% parameters of MueMat.
%
% See also: Tutorial, DefaultParameters, Laplace1D

%% Header
clear all;
mue_include;

%% Setup the multigrid Solver
Amat = gallery('poisson',30); % The test matrix

% Create a multigrid solver with default parameters
mySolver = Hierarchy(Amat);
mySolver.SetOutputLevel(1);
status = mySolver.FillHierarchy();
mySolver.SetSmoothers();


%% Generate a random right-hand side ('b')
% 'srand' initializes the random number generator to ensure reproducible results
srand;
n = size(Amat,1);
b = rand(n,1);

%% Solve Ax=b (AMG used directly as a solver)
% A fixed number of iterations is performed
nIts = 9;
disp(norm(b))
fprintf('Solving Ax=b with multigrid...\n\n');
x    = mySolver.Iterate(b, nIts);

fprintf('\n||r_0|| / ||r_final|| = %g\n\n',norm(b-Amat*x)/norm(b));

%% Solve Ax=b (AMG as a preconditioner of an iterative method)
maxIts = nIts;
tol  = 1e-8;

mySolver.SetOutputLevel(0);
fprintf('Solving Ax=b again with conjugate gradients preconditioned by multigrid...\n\n');
[x,flag,relres,iter,resvec] = ...                                % pcg() is the Matlab conjugate gradients method.
    pcg(                           ... %
        Amat, b,                   ... % parameters: * Matrix and right-hand side
        tol, maxIts,               ... %             * Conjugate Gradient parameters
        @(v)mySolver.Iterate(v, 1) ... %             * AMG Preconditioner
       );
for ii=1:length(resvec),
  fprintf('  %d: ||r||=%g\n',ii,resvec(ii));
end

fprintf('\n||r_0|| / ||r_final|| = %g\n\n',norm(b-Amat*x)/norm(b));
