% this script is currently not in the non-regression test
% todo: refactorize reuse scenarios
function TestReUseOptions(fast)
  if ~varexist('fast'), fast = true; end;

   if fast
    
    for numDesiredLevel=3 % it's important to test with numLevels>2
      TestSaPFactoryReUse(false, false, false, false, numDesiredLevel); % no reuse
      TestSaPFactoryReUse(true,  false, false, false, numDesiredLevel); % reuse P

      % FIXME: no reuse of Ptent possible...
      %TestSaPFactoryReUse(false, true,  false, false, numDesiredLevel); % reuse Ptent
      TestSaPFactoryReUse(false, false, true,  false, numDesiredLevel); % reuse Aggregates
      TestSaPFactoryReUse(false, false, false, true,  numDesiredLevel); % reuse Graph
    end
    
  else
    
    for reUseP=[false, true]
      for reUsePtent=[false, true]
        for reUseAggregates=[false, true]
          for reUseGraph=[false, true]
            for numDesiredLevel=2:3
              TestSaPFactoryReUse(reUseP, reUsePtent, reUseAggregates, reUseGraph, numDesiredLevel);
            end
          end
        end
      end
    end
    
  end
  
end  

function TestSaPFactoryReUse(ReUseP, ReUsePtent, ReUseAggregates, ReUseGraph, numDesiredLevels)
  fprintf('TestSaPFactoryReUse(%d, %d, %d, %d, %d)', ReUseP, ReUsePtent, ReUseAggregates, ReUseGraph, numDesiredLevels);
   
  srand;
  
%   ReUseP = 0;
%   ReUsePtent = 0;
%   ReUseAggregates = 1;
%   ReUseGraph = 0;
  
  % Solve a Laplace 1D problem
  n      = 100; 
  A      = BuildLaplace1D(n);
  b      = rand(n,1);
  numDesiredLevels = 3;

  mue_include

  %% Set AMG options

  AmalgamateDropFact = CoalesceDropFactory();
  AggFact            = AggregationFactory();
  PtentFact          = TentativePFactory(AmalgamateDropFact,AggFact);
  Pfact              = SaPFactory(PtentFact);
  Rfact              = TransPFactory();
  PRfact             = GenericPRFactory(Pfact,Rfact);
  PRfact.SetMaxCoarseSize(30);
  Acfact             = RAPFactory();
  SmooFactory        = SmootherFactory(Smoother('GaussSeidel', 2, 1));  

  %% RUN 1
  % Create the multigrid hierarchy
  % Construct and populate finest level with user information
  Finest = Level(); Finest.KeepAll(false); 
  Finest.Set('A', A);
  Finest.Set('NullSpace', BuildNullSpace(A));
  Finest.Keep('NullSpace'); % Keeped for the second run
  
  CoarseLevel = Level(); CoarseLevel.KeepAll(false);
  % Note: Instead of creating a first coarselevel, users can set the Keep()
  % options directly on the finest level.
  
  if ReUseP
    CoarseLevel.Keep('P',Pfact); % it's the default ...
    CoarseLevel.Keep('NullSpace');
  end
  if ReUsePtent 
    Finest.Keep('P',PtentFact);
    CoarseLevel.Keep('P',PtentFact);
    CoarseLevel.Keep('NullSpace');
  end
  if ReUseAggregates
    Finest.Keep('Aggregates',AggFact);
    CoarseLevel.Keep('Aggregates',AggFact);
  end
  if ReUseGraph
    Finest.Keep('Graph',AmalgamateDropFact);
    CoarseLevel.Keep('Graph',AmalgamateDropFact);
  end
  
  MgHierarchy = Hierarchy(); MgHierarchy.SetOutputLevel(1);
  MgHierarchy.SetLevel(Finest,1);
  MgHierarchy.SetLevel(CoarseLevel,2);
  MgHierarchy.FillHierarchy(PRfact,[], Acfact, 1, numDesiredLevels);
  MgHierarchy.SetSmoothers(SmooFactory);
  
  if (MgHierarchy.GetNumLevel() ~= numDesiredLevels)
    error('MgHierarchy.GetNumLevel() ~= numDesiredLevels: change the problem size for this test');
  end
  
  % AMG as a preconditioner to CG
  maxIts = 9;     % Maximum number of iterations
  tol    = 1e-12; % Tolerance parameter

  x    = zeros(n,1); SolStatus = ALLZEROS;
  zero = zeros(n,1);
  [x, flag, relRes, nIts, resVec] = pcg(A.GetMatrixData(),b,tol,maxIts,@(b)MgHierarchy.Iterate(b,1, zero,ALLZEROS)); SolStatus = NOTALLZEROS;

  %% RUN 2
  
  AmalgamateDropFact2 = CoalesceDropFactory();
  AggFact2            = AggregationFactory();
  PtentFact2          = TentativePFactory(AmalgamateDropFact2,AggFact2);
  Pfact2              = SaPFactory(PtentFact2);
  Rfact2              = TransPFactory();%GenericRFactory(Pfact2);
  Acfact2             = RAPFactory();
  

  if ReUsePtent, PtentFact2 = PtentFact; Pfact2 = PgPFactory(PtentFact2); Rfact2 = GenericRFactory(Pfact2); end;
  if ReUseAggregates, AggFact2 = AggFact; PtentFact2 = TentativePFactory(AmalgamateDropFact2,AggFact2); end;
  if ReUseGraph, AmalgamateDropFact2 = AmalgamateDropFact; end;
  PtentFact2          = TentativePFactory(AmalgamateDropFact2,AggFact2);
  Pfact2              = SaPFactory(PtentFact2);
  if ReUseP,     Pfact2 = Pfact; Rfact2 = Rfact;  end;
  
  clear PRfact;
  PRfact2      = GenericPRFactory(Pfact2,Rfact2);
  PRfact2.SetMaxCoarseSize(30);
  
  % Some cleanup to avoid warnings: Coarse A, P and R are always kept in a Level but
  % will be overwritten in the second run... So I manually removed them...
  % No idea how to do that cleaner
  for i=2:MgHierarchy.GetNumLevel()
    MgHierarchy.Levels_{i}.Delete('A'); MgHierarchy.Levels_{i}.Keep('A');
    if ~ReUseP, MgHierarchy.Levels_{i}.Delete('P'); MgHierarchy.Levels_{i}.Keep('P'); end;
    %if ~ReUsePtent, MgHierarchy.Levels_{i}.Delete('P',PtentFact2); MgHierarchy.Levels_{i}.Keep('P',PtentFact2); end;
    MgHierarchy.Levels_{i}.Delete('R'); MgHierarchy.Levels_{i}.Keep('R');
  end
  
  MgHierarchy.FillHierarchy(PRfact2,[], Acfact2, 1, numDesiredLevels);

  % AMG as a preconditioner to CG
  x    = zeros(n,1); SolStatus = ALLZEROS;
  zero = zeros(n,1);
  [x, flag, relRes, nIts, resVec] = pcg(A.GetMatrixData(),b,tol,maxIts,@(b)MgHierarchy.Iterate(b,1, zero,ALLZEROS)); SolStatus = NOTALLZEROS;
end