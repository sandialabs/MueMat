%% Example 1
% This minimalistic example solves a linear system $Ax=b$ using smoothed
% aggregation multigrid
%

%% Header
clear all;
mue_include;

%% Setup the multigrid Solver
Amat = gallery('poisson',50); % The test matrix
Amat = Operator(Amat,1,1,@MatlabApply); % transform matrix Amat to MueMat operator

% create SaPFactory for smoothed aggregation transfer operators
Pfact = SaPFactory(TentativePFactory());
Rfact = TransPFactory();    % use the transposed of P for restriction

% define level smoothers
Sfact = SmootherFactory(Smoother('Jacobi',3,1.0));

% define finest multigrid level
Finest = Level();
Finest.KeepAll(true);   %% free variables as early as possible
Finest.Set('A',Amat);   % set fine level matrix

MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(10); % be verbose
MgHierarchy.SetLevel(Finest,1);
MgHierarchy.SetMaxCoarseSize(1);
status = MgHierarchy.FillHierarchy(Pfact, Rfact, RAPFactory(), 1, 3);
MgHierarchy.SetSmoothers(Sfact);


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
