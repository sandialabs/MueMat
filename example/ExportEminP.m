% Build an energy-minimizing prolongator that can be used outside of
% MueMat for experiments. Factories are called directly instead of
% running the Hierarchy setup phase.

% Test driver
function EminP = ExportEminP(FineNullSpace, CoarseNullSpace, Ppattern_in, PinitFact_in, EnergyMatrix_in)

% TODO Nullspace

% Data storage
fineLevel   = Level();
coarseLevel = Level();
fineLevel.KeepAll();
coarseLevel.KeepAll();

% NullSpace
fineLevel.Set('NullSpace', FineNullSpace);
coarseLevel.Set('NullSpace', CoarseNullSpace);

% Sparsity pattern of P
%  The sparsity pattern is already available, so we don't need a
%  factory to build it. The class NoFactory is used as a
%  placeholder.
PatternFact = NoFactory('Pattern Factory');
coarseLevel.Set('Ppattern', Ppattern_in, PatternFact);

% Initial guess for the energy minimization process
InitPFact = NoFactory('Initial P guess');
coarseLevel.Set('P', PinitFact_in, InitPFact);

% Emin solver
eminSteps  = 2; % number of energy mininization iteration
EminSolver = CGEminSolver(eminSteps);

% Options - I think this argument is unused
options = [];

% Matrix defining the energy norm
EnergyMatrix = 'EnergyMatrix'; % == string
fineLevel.Set(EnergyMatrix, EnergyMatrix_in);

% Factory building the energy-minimizing prolongator
% Args: EminPFactory(PatternFact, ConstraintFact, EminSolver, InitPFact, options, EnergyMatrix, diagonalView)
Pfact = EminPFactory(PatternFact, ConstraintFactory(), EminSolver, InitPFact, options, EnergyMatrix);

% Build the energy-minimizing prolongator
Pfact.Build(fineLevel, coarseLevel);
EminP = coarseLevel.Get('P', Pfact);

end
