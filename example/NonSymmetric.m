function NonSymmetric(param)
if ~varexist('param'), param = 0, end

srand;
%clear all;
mue_include
SolStatus = NOTALLZEROS;

n       = 100;

Amat = BuildLaplace1D(n); %TODO: change this problem for an unsymmetric problem.

% Set options
numDesiredLevels = 4;         % number of AMG levels

AmalgamateDropFact= CoalesceDropFactory();
AggFact        = AggregationFactory();
PtentFact      = TentativePFactory(AmalgamateDropFact,AggFact);
Pfact          = SaPFactory(PtentFact);
if param % DEBUG
  Rfact             = TransPFactory();
else

  Rfact            = GenericRFactory(Pfact);
end

Acfact            = RAPFactory();
GSFactory         = SmootherFactory(Smoother('GaussSeidel', 2, 1));

%
%  Construct and populate finest level with user information
%
Finest = Level();
Finest.Set('A', Amat);
Finest.Set('NullSpace', BuildNullSpace(Amat));

MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(1);
MgHierarchy.SetLevel(Finest,1);
MgHierarchy.FillHierarchy(Pfact, Rfact, Acfact, 1, numDesiredLevels);
MgHierarchy.SetSmoothers(GSFactory);

rhs = rand(Amat.GetRowMap().NDOFs(),1);
sol = zeros(Amat.GetRowMap().NDOFs(),1);           SolStatus = ALLZEROS;

%% Solve Ax=b (AMG as a preconditioner of an iterative method)
maxIts = 10;
tol  = 1e-12;

x = ...                                    % pcg() is the Matlab conjugate gradients method.
    pcg(                               ... %
        Amat.GetMatrixData(), rhs,     ... % parameters: * Matrix and right-hand side
        tol, maxIts,                   ... %             * Conjugate Gradient parameters
        @(b)MgHierarchy.Iterate(b, 1)  ... %             * AMG Preconditioner
       );
