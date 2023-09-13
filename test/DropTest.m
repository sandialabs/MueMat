% Load Amat and rhs:
load('data/SmallHetero'); Amat = Operator(Amat);

mue_include;

numDesiredLevels = 4;         % number of multigrid levels

AmalgamateDropFact= CoalesceDropFactory();
AmalgamateDropFact.SetPostDropSpecifications([],[]);
AmalgamateDropFact.SetPreDropSpecifications(@ExperimentalDrop,.185);
AggFact        = AggregationFactory();
AggFact.SetRemoveAggsThisSizeOrLess(1);
if (1 == 0) 
   % Use Smoothed Aggregation
   Ptentfact      = TentativePFactory(AmalgamateDropFact,AggFact);
   Pfact          = SaPFactory(Ptentfact);
   Pfact.SetUseAfiltered(true);
else
  % Use Emin
   Patfact             = AP_PatternFactory();
   Patfact.SetPatternType('AP');
   Patfact.SetUseAfiltered(true);
   CSFactory           = ConstraintFactory();
   options.NCoarseDofPerNode = 1;
   options.CoarseNullspaceWeightingMultiplier = [1];
   options.SmoothNullspace = 1;
   Pinit = TentativePFactoryEx(AmalgamateDropFact, AggFact, SaCoarseNSFactory, options);
   Pfact = EminPFactory(Patfact, CSFactory, CGEminSolver(5), Pinit, options);
   %Pfact.SetSolverIterations(5);
end
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
status = MgHierarchy.FillHierarchy(GenericPRFactory(Pfact),[], Acfact, 1, numDesiredLevels);
MgHierarchy.SetSmoothers(GSFactory);

params(1) = 1e-7; params(2) = 100;

[mysol,flag,relres,iter,resvec]=pcg(Amat.GetMatrixData(), rhs, 1e-7, 100, @(v)MgHierarchy.Iterate(v,1,zeros(length(rhs),1),ALLZEROS,2));

fprintf('Nits = %d,||r||=%e\n',length(resvec)-1,norm(rhs-Amat.GetMatrixData()*mysol)/norm(rhs));
