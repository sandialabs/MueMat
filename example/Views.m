% Example to illustrate how to use different views of an Operator.
%
mue_include

Amat = BuildLaplace1DBlk(1,-1,100);
rowmap = Amat.GetRowMap();
srand; rhs  = rand(rowmap.NDOFs(),1);
sol  = zeros(rowmap.NDOFs(),1); InitGuessStatus = ALLZEROS;

%
% Factories to define AMG
%
Pfact              = SaPFactory();
Pfact.SetDiagonalView('point'); %important for this test
Rfact              = TransPFactory();
PRfact             = GenericPRFactory(Pfact,Rfact);
Acfact             = RAPFactory();
%
%  Construct and populate finest level with user information
%
Finest = Level();
Finest.KeepAll(false);
Finest.Set('A', Amat);
Finest.Set('NullSpace', ones(rowmap.NDOFs(),1));

MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(1);
MgHierarchy.SetLevel(Finest,1);
MgHierarchy.FillHierarchy(PRfact,[], Acfact, 1, 2);

sweeps = 1;               % number of smoothing sweeps (pre or post)
startLevel = 1;           % start indexing at level 1
numDesiredLevels = 2;     % max number of levels in hierarchy
iterations = 9;           % multigrid iterations

% First run point relaxation.
fprintf('Running symmetric Gauss-Seidel V(1,1) with w = .99\n');
SFact   = SmootherFactory(Smoother('GaussSeidel', sweeps, .99));
MgHierarchy    = MgHierarchy.Copy();
MgHierarchy.SetSmoothers(SFact);
newsol   = MgHierarchy.Iterate(rhs, iterations, sol, InitGuessStatus);

% Now run block relaxation, where the diagonal blocks are 2x2.
fprintf('Running symmetric Block Gauss-Seidel V(1,1) with w = .99\n');
fprintf('diagonal blocks are 2x2\n');
FineA    = MgHierarchy.GetLevel(1).Get('A');
% TODO we should be able to query Operator FineA to get #rows, etc.
map = Map(size(FineA.GetMatrixData(),1)/2,2);
FineA.CreateView('constblock', map, map, FineA.GetApply());
previousView = FineA.SwitchToView('constblock');
% create a new block diagonal
Diag = FineA.GetDiagonal();
SFact   = SmootherFactory(Smoother('GaussSeidel', sweeps, .99));
MgHierarchy.SetSmoothers(SFact);
newsol   = MgHierarchy.Iterate(rhs, iterations, sol, InitGuessStatus);

% Now switch back to the point view.
% Iterates should be the same as using point Gauss-Seidel.
fprintf('Running symmetric Gauss-Seidel V(1,1) with w = .99\n');
fprintf('diagonal blocks are 1x1\n');
FineA = MgHierarchy.GetLevel(1).Get('A');
FineA.SwitchToView(previousView);
SFact   = SmootherFactory(Smoother('GaussSeidel', sweeps, .99));
MgHierarchy.SetSmoothers(SFact);
newsol   = MgHierarchy.Iterate(rhs, iterations, sol, InitGuessStatus);


% TODO: add examples without view switching (using Smoother(...'constblock'))
