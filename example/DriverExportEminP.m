clear all;
mue_include;

dim = 1;
Amat = gallery('poisson',30);
A = Operator(Amat, dim, dim);
FineNullSpace = BuildNullSpace(A);

% Data storage
fineLevel   = Level();
coarseLevel = Level();
fineLevel.KeepAll();
coarseLevel.KeepAll();

% Populate fine level
fineLevel.Set('A', A);
fineLevel.Set('NullSpace', FineNullSpace);

% Compute InitPFact
InitPFact = TentativePFactory();
InitPFact.Build(fineLevel, coarseLevel);
InitP = coarseLevel.Get('P', InitPFact);
CoarseNullSpace = coarseLevel.Get('NullSpace'); % should be CoarseNullSpace = coarseLevel.Get('NullSpace', InitPFact);

% Compute Ppattern
PpatternFact = AP_PatternFactory([], InitPFact);
PpatternFact.Build(fineLevel, coarseLevel);
Ppattern = coarseLevel.Get('Ppattern', PpatternFact);

% Get EminP
EminP = ExportEminP(FineNullSpace, CoarseNullSpace, Ppattern, InitP, A);