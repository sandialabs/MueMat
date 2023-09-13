% Preconditioning indefinite Helmholtz systems with AMG
%
% Use Shifted Laplacian-AMG as a preconditioner to GMRES.
%
% Notes: Prolongation/Restriction operators are built with just the stiffness matrix.
%        Complex shifts are different at every level (uses RAPShiftFactory)
%
% See BuildHelmholtz2D

clear all;
srand;
mue_include

i=sqrt(-1);
alphas = [1 1 1 1 1 1 1];
%betas  = [0.5 0.4 0.3 0.2 0.1];
%betas  = alphas/2;
betas   = [0.5 0.7 0.9 1.1 1.3 1.5 1.7];
shifts = alphas+i*betas;

% Build a 2D Helmholtz problem
waves=40; npts=10; medium=0; N=waves*npts;
lambda=1/waves;
omega=2*pi/lambda;
[Kmat, Mmat, Nullspace, Nx, Ny] = BuildHelmholtz2D(waves,npts,medium);

% Form sparse matrices and operators
Amat = Kmat-(omega^2)*Mmat;
Smat = Kmat-shifts(1)*(omega^2)*Mmat;
A=Operator(Amat,1,1); S=Operator(Smat,1,1);
K=Operator(Kmat,1,1); M=Operator(Mmat,1,1);

[nrows,ncols]=size(Amat);
nDOFS=nrows

%----------------------------------------------
% Smoothed Aggregation
%----------------------------------------------

% maximum number of AMG level
numDesiredLevels = length(shifts);

% Setup factories
AmalgamateDropFact = CoalesceDropFactory();
AggFact            = AggregationFactory();
Ptentfact          = TentativePFactory(AmalgamateDropFact,AggFact); Ptentfact.TentativeWithQR(0);
Pfact              = SaPFactory(Ptentfact); Pfact.SetAForSmoothing('Afiltered');
Rfact              = TransPFactory();
Acfact             = RAPShiftFactory(); % Special RAP factory for Helmholtz
SmooFactory        = SmootherFactory(ILUSmoother());
Acfact.SetOmega(omega);
Acfact.SetShifts(shifts);

% Fill fine leve information
Finest = Level();
Finest.Keep('A');
Finest.Set('A', S);
Finest.Keep('Afiltered');
Finest.Set('Afiltered', K);
Finest.Keep('M');
Finest.Set('M', M);

% Create a multigrid hierarchy
MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(3);
MgHierarchy.SetLevel(Finest,1);
MgHierarchy.FillHierarchy(Pfact, Rfact, Acfact, 1, numDesiredLevels);
MgHierarchy.SetSmoothers(SmooFactory);

% Iterative solve
% GMRES parameters
maxIts  = 100;   
tol     = 1e-6; 
restart = 100;
params(1) = tol; params(2) = 1000; params(3)=1;

% AMG as a preconditioner to GMRES or FGMRES (if GMRES is smoother)
x    = zeros(nDOFS,1);
zero = zeros(nDOFS,1);
b    = zeros(nDOFS,1);
h    = lambda/npts;
idx  = reshape(1:nDOFS,N,N);
b(idx(end/2,end/2)) = 1/(h*h);
%[x, flag, relRes, nIts, resVec] = gmres(A.GetMatrixData(),b,restart,tol,maxIts,@(v)MgHierarchy.Iterate(v,1, zero,ALLZEROS));
[x, flag, relRes, nIts, resVec]  = bicgstab(A.GetMatrixData(),b,tol,maxIts,@(v)MgHierarchy.Iterate(v,1, zero,ALLZEROS));
fprintf('converged to the desired tolerance %g within %d iterations.\n', tol, nIts);
imagesc(real(reshape(x,size(idx)))); colorbar; axis equal; axis tight;
