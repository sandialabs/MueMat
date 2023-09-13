% Load Amat and rhs:
load('data/SmallHetero'); Amat = Operator(Amat);

mue_include;

numDesiredLevels = 3         % number of multigrid levels

AmalgamateDropFact= CoalesceDropFactory();
AmalgamateDropFact.SetPostDropSpecifications([],[]);
AmalgamateDropFact.SetPreDropSpecifications(@ExperimentalDrop,.185);
AggFact        = AggregationFactory();
AggFact.SetRemoveAggsThisSizeOrLess(1);
FCFact         = FCSplittingFactory(AggFact);

Patfact             = AP_PatternFactory();
Patfact.SetPatternType('AP');
% Patfact.SetUseAfiltered(true);
NSFactory           = CoarseNSFactory();
options.NCoarseDofPerNode = 1;
options.PtentSA = 1;
options.CoarseNullspaceWeightingMultiplier = [1];
options.DropConstraintsMethod = 'nullspace';
Pinitfact = TentativePFactoryEx(CoalesceDropFactory(), AggFact, NSFactory);
Pfact    = EminPFactory(Patfact, ConstraintFactory([],FCFact,options), GMRESEminSolver(5), Pinitfact); 
AltPfact = EminPFactory(Patfact, ConstraintFactory([],FCFact,options), GMRESEminSolver(5), Pinitfact);

Rfact             = GenericRFactory(AltPfact);
% PRfact             = GenericPRFactory(Pfact,Rfact);
Acfact            = RAPFactory();
GSFactory         = SmootherFactory(Smoother('GaussSeidel', 1, 1.25));

%
%  Construct and populate finest level with user information
%
Finest = Level();
Finest.KeepAll(false);
Finest.Set('A', Amat);

MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(1);
MgHierarchy.SetLevel(Finest,1);
Needs.SaveAggregates = 'all';
status = MgHierarchy.FillHierarchy(Pfact, Rfact, Acfact, 1, numDesiredLevels);
MgHierarchy.SetSmoothers(GSFactory, 1, numDesiredLevels);

[sol,flag,relres,iter,resvec]=gmres(Amat.GetMatrixData, rhs,[],1e-8,100, @(v)MgHierarchy.Iterate(v,1,zeros(length(rhs),1),ALLZEROS,2),[],zeros(length(rhs),1));

fprintf('Nits = %d, ||r|| = %e \n', length(resvec)-1,...
norm(rhs-Amat.GetMatrixData()*sol)/norm(rhs));
