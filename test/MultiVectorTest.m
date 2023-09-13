function MultiVectorTest()
  % Test MultiVector in a variety of settings.  This is just an amalgation of
  % SmoothingTest.m and ElasticityTest.m.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Smoothing Test
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  mue_include
  SetHomeDir
  InitGuessStatus = NOTALLZEROS;

  newstream = RandStream.create('mrg32k3a','NumStreams',1);
  oldStream = RandStream.setGlobalStream(newstream);

  Amat = BuildLaplace1DBlk(1,-1,100);

  rowmap = Amat.GetRowMap();
  %rhs  = rand(RandStream.create('mrg32k3a','NumStreams',1),rowmap.NDOFs(),1);
  %sol  = zeros(rowmap.NDOFs(),1); InitGuessStatus = ALLZEROS;
  rhs = MultiVector(rowmap.NDOFs(),1);
  rhs.Random();
  sol  = MultiVector(rowmap.NDOFs(),1);
  sol.PutScalar(0); InitGuessStatus = ALLZEROS;

  %%
  % *Factories to define multigrid*
  %
  AmalgamateDropFact= CoalesceDropFactory();
  AggFact        = AggregationFactory();
  Ptentfact         = TentativePFactory(AmalgamateDropFact,AggFact);
  Pfact             = SaPFactory(Ptentfact);
  Rfact             = TransPFactory();
  PRfact            = GenericPRFactory(Pfact,Rfact);
  Acfact            = RAPFactory();
  %%
  % *Construct and populate finest level with user information*
  %
  Finest = Level();
  Finest.Set('A', Amat);
  %Finest.Set('NullSpace', ones(rowmap.NDOFs(),1));
  Finest.Set('NullSpace', BuildNullSpace(Amat));

  MgHierarchy = Hierarchy();
  MgHierarchy.SetOutputLevel(1);
  MgHierarchy.SetLevel(Finest,1);
  MgHierarchy.FillHierarchy(PRfact, [], Acfact, 1, 2);

  %%
  % *Different ways to set and exercise smoothers*
  %
  sweeps = 1;               % number of smoothing sweeps (pre or post)
  startLevel = 1;           % start indexing at level 1
  numDesiredLevels = 2;     % max number of levels in hierarchy
  iterations = 9;           % multigrid iterations
  smoothCoarsest = false;   % if true, run smoother on coarsest problem
                          % (must be true if a one-level method and you
                          %  don't want a direct solve)

  %%
  % *Symmetric Gauss-Seidel*
  fprintf('Running symmetric Gauss-Seidel V(1,1) with w = .99\n');
  SFact   = SmootherFactory(Smoother('GaussSeidel', sweeps, .99, 'point'));
  NewMg = MgHierarchy.Copy();
  NewMg.SetSmoothers(SFact);
  if (smoothCoarsest), NewMg.SetCoarsestSolver(SFact); end;
  newsol   = NewMg.Iterate(rhs, iterations, sol, InitGuessStatus);

  %%
  % *Symmetric Block Gauss-Seidel, pre and post smoothing*
  fprintf('Running symmetric Block Gauss-Seidel V(1,1) with w = .99\n');
  SFact   = SmootherFactory(Smoother('GaussSeidel', sweeps, .99, 'default'));
  NewMg = MgHierarchy.Copy();
  NewMg.SetSmoothers(SFact);
  if (smoothCoarsest), NewMg.SetCoarsestSolver(SFact); end;
  newsol   = NewMg.Iterate(rhs, iterations, sol, InitGuessStatus);

  %%
  % *Jacobi*
  fprintf('Running Jacobi V(1,1) with w = .7\n');
  SFact  = SmootherFactory(Smoother('Jacobi', sweeps, .7, 'point'));
  NewMg = MgHierarchy.Copy();
  NewMg.SetSmoothers(SFact);
  if (smoothCoarsest), NewMg.SetCoarsestSolver(SFact); end;
  newsol   = NewMg.Iterate(rhs, iterations, sol, InitGuessStatus);

  %%
  % *Block Jacobi*
  fprintf('Running Block Jacobi V(1,1) with w = .7\n');
  SFact  = SmootherFactory(Smoother('Jacobi', sweeps, .7, 'default'));
  NewMg = MgHierarchy.Copy();
  NewMg.SetSmoothers(SFact);
  if (smoothCoarsest), NewMg.SetCoarsestSolver(SFact); end;
  newsol   = NewMg.Iterate(rhs, iterations, sol, InitGuessStatus);

  %%
  % *Overlapping additive Schwarz, pre and post smoothing, random domains*
  fprintf('Running overlapping additive domain decomp. V(1,1) where domains are chosen randomly\n');
  SFact  = SmootherFactory(Smoother('Jacobi', sweeps, .7,'Random NonOverlapping'));
  NewMg = MgHierarchy.Copy();
  NewMg.SetSmoothers(SFact);
  if (smoothCoarsest), NewMg.SetCoarsestSolver(SFact); end;
  newsol   = NewMg.Iterate(rhs, iterations, sol, InitGuessStatus);

  %%
  % *Overlapping multiplicative Schwarz, pre and post smoothing, random domains*
  fprintf('Running symmetric overlapping multiplicative domain decomp. V(1,1) where domains are chosen randomly\n');
  SFact   = SmootherFactory(Smoother('GaussSeidel', sweeps, .99,'Random Overlapping'));
  NewMg = MgHierarchy.Copy();
  NewMg.SetSmoothers(SFact);
  if (smoothCoarsest), NewMg.SetCoarsestSolver(SFact); end;
  newsol   = NewMg.Iterate(rhs, iterations, sol, InitGuessStatus);

  %%
  % *Alternative setup method for overlapping multiplicative Schwarz, pre and post smoothing, random domains*
  fprintf('Running symmetric overlapping multiplicative domain decomp. V(1,1) where domains are chosen randomly but ...\n');
  fprintf('setting the blocks in a more direct but low-level fashion\n');
  NewMg = MgHierarchy.Copy();
  Finest      = NewMg.GetLevel(1);
  FineA       = Finest.Get('A');
  
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Elasticity Test
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  dim = 2;
  NLIST = 40;
  stretch = 1;

  %% Set test case parameters
  % Parameters for 304 Stainless Steel
  ELASTIC_MODULUS=193e9; 
  POISSONS_RATIO=.305;
  
  global numNS_EMIN;
  numNS       = 3;
  numNS_NOROT = 2;
  numNS_EMIN  = 2;
  stretchvec = [stretch 1];
  wgts = [1 1 1];

  %% Set solver parameters
  eminSteps = 10;
  AggDropTol=.01;

  %% Init output
  ITERS=zeros(length(NLIST),4);
  OC=ITERS;
  NLEVELS=ITERS;

  for I=1:length(NLIST),
  n=NLIST(I);
  nstr = num2str(n);
  sstr = num2str(stretch);
  dimstr = num2str(dim);

  if dim == 2, sizestr = [nstr 'x' nstr]; else sizestr = [nstr 'x' nstr 'x' nstr]; end

  fprintf('***Running %s %d/%d***\n',sizestr,I,length(NLIST));

  % Build matrix, block it and grab parameters
  filename = [MUEMAT_ROOT_DIR '/data/Elas' dimstr 'D-' sizestr];
  if(stretch==1)
    filename = [filename '.mat'];
  else % dim ==3
    filename = [filename '.s.' sstr '.mat'];    
  end
  
  fid = fopen(filename,'r');
  if fid < 0
    if dim == 2
      [Amat,nullspace,coords]=BuildElasticity2D(n,ELASTIC_MODULUS, POISSONS_RATIO,stretchvec);  
    else % dim == 3
      [Amat,nullspace,coords]=BuildElasticity3D(n,ELASTIC_MODULUS, POISSONS_RATIO,stretchvec);  
    end
    Amat = Amat.GetMatrixData();
    %save(filename, 'Amat', 'nullspace', 'coords');
  else
    fclose(fid);
    load(filename);
  end
  Amat =  Operator(Amat,dim,dim);
  RowMap = Amat.GetRowMap();
  ndofs = RowMap.NDOFs();

  % Create RHS, run parameters
  newstream = RandStream.create('mrg32k3a','NumStreams',1);
  oldStream = RandStream.setGlobalStream(newstream);
  rhs = MultiVector(RowMap.NDOFs(),1);
  rhs.Random();
  %rhs = rhs / norm(rhs);  %FIXME can't do this right now
  maxit=50;
  tol = 1e-10; 
  zeroGuess  = MultiVector(RowMap.NDOFs(),1);
  zeroGuess.PutScalar(0);

  % Run Smoothed aggregation
  J=2;
  LABEL{J}='sa';
  [MgHierarchy,OC(I,J)]=build_amg_hierarchy(dim,Amat.Copy(),nullspace(:,1:numNS),'sa',numNS,wgts,eminSteps,coords,AggDropTol);
  [sol,flag,relres,ITERS(I,J),RESID{I,J}] = pcg(Amat.GetMatrixData(),rhs.GetVectors(),tol,maxit,...
                                          @(b)MgHierarchy.Iterate(b,1, zeroGuess.GetVectors(),ALLZEROS,VCYCLE));
  
  % Run Ray's special sauce
  J=3;
  LABEL{J}=['svd-' num2str(numNS_EMIN)];
  [MgHierarchy,OC(I,J),NLEVELS(I,J)]=build_amg_hierarchy(dim,Amat.Copy(),nullspace(:,1:numNS),'svd',numNS,wgts,eminSteps,coords,AggDropTol);
  [sol,flag,relres,ITERS(I,J),RESID{I,J}] = pcg(Amat.GetMatrixData(),rhs.GetVectors(),tol,maxit,...
                                          @(b)MgHierarchy.Iterate(b,1, zeroGuess.GetVectors(),ALLZEROS,VCYCLE));

  % Output diagnostics
  fprintf('\nstatistics: ');
  for J=1:length(LABEL),fprintf('%s %2d(%3.2f)   ',LABEL{J},ITERS(I,J),OC(I,J));end
  fprintf('\n');

  end

  fprintf('DIM=%d\n',dim);
  fprintf('STRETCH=%f\n',stretch(1));

  % Useful for printing
  %fprintf('%2d %2d %4.2f\n',[NLIST',ITERS(:,3),OC(:,3)]')
  fprintf('\n');
  fprintf('       |  SA-NR  |    SA   |   EMIN  |    ML  \n');
  fprintf('SZ LVL | ITS  OC | ITS  OC | ITS  OC | ITS  OC\n');
  fprintf('----------------------------------------------\n');
  fprintf('%2d %2d  | %2d %4.2f | %2d %4.2f | %2d %4.2f | %2d %4.2f\n',[NLIST',NLEVELS(:,3),ITERS(:,1),OC(:,1),ITERS(:,2),OC(:,2),ITERS(:,3),OC(:,3),ITERS(:,4),OC(:,4)]')

  % Reset random number generator.
  RandStream.setGlobalStream(oldStream);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MgHierarchy,oc,NLevels]=build_amg_hierarchy(dim, Amat,nullspace,method_type,numNS,wgts,eminSteps,coords,AggDropTol)
  %  Construct and populate finest level
  Finest = Level();
  Finest.Set('A', Amat);
  Finest.Set('NullSpace', nullspace);
  Finest.Set('xcoords', coords(:,1));
  Finest.Set('ycoords', coords(:,2));
  if dim == 3
    Finest.Set('zcoords', coords(:,3));
  end
  
  % Set options
  numDesiredLevels = 10;
  %aggoptions.AggTargetSize = 3;
  %SAoptions.NoQR=1;
  %SAoptions.CNullMode='inject partialFN';
  SAoptions.CNullMode='inject CMS';
  %SAoptions.CNullMode='inject SVD';
  
  % Build agglomeration factory
  AmalgamateDropFact = CoalesceDropFactory();

  % Amalgamate Amat to build auxillary coordinate Laplacian for aggregation
  % Note: We do this *before* we do the dropping since we want the
  % full distance Laplacian.  Also set the AuxMatrixFunc to
  % regenerate (rather than coarsen) the distance Laplacian.
  AmalgMat = AmalgamateDropFact.Build_(Amat);
  AuxMat = BuildDistanceLaplacian(AmalgMat,coords);
  Finest.Set('AuxMatrix', AuxMat);
  Finest.Set('AuxMatrixFunc', @BuildDistanceLaplacian);      
  AmalgamateDropFact.SetAName('AuxMatrix');

  
  % Set dropping
  AmalgamateDropFact.SetPostDropSpecifications(@SAdrop,AggDropTol);

  % Aggregation
  % Note: Dropping should be done in ADF, not ML aggregation
  AggFact = AggregationFactory();
  %AggFact = AggFact.SetAlgorithm('ml'); % use the same aggregation as ML (need MLMEX)
  %AggFact = AggFact.SetMLOptions('%ML output',10);

  if(strcmp(method_type,'sa')),
    % Standard Smoothed Aggregation
    Ptentfact          = TentativePFactory(AmalgamateDropFact,AggFact);
    Pfact              = SaPFactory(Ptentfact);
  elseif(strcmp(method_type,'svd')),
     % Ray's PtentExperiment
    global numNS_EMIN
    SAoptions.NCoarseDofPerNode=numNS_EMIN;
    wgts = wgts(1:numNS);
    fprintf('weights = [%s]\n',num2str(wgts));
    SAoptions.CoarseNullspaceWeightingMultiplier = wgts;
    SAoptions.DropConstraintsMethod='nullspace';
    %SAoptions.PFatten = true; 
    SAoptions.FilteredAP=true;
    %SAoptions.PtentRootModifications = '4b';
    CNSFact = CoarseNSFactory();

    Pinitfact      = TentativePFactoryEx(AmalgamateDropFact, AggFact, CNSFact, SAoptions);
    PatFact        = AP_PatternFactory([],Pinitfact,SAoptions); % no filter, only options
    Pfact          = EminPFactory(PatFact,ConstraintFactory(),CGEminSolver(eminSteps), Pinitfact, SAoptions);
  else
    fprintf('ERROR: Bogus Method type\n');
    return;    
  end
  PRfact = GenericPRFactory(Pfact);
  PRfact.SetMaxCoarseSize(100);
  Rfact              = TransPFactory();
  Acfact             = RAPexFactory();
  GSFactory          = SmootherFactory(Smoother('GaussSeidel', 2, 1));
  
  MgHierarchy = Hierarchy();
  MgHierarchy.SetOutputLevel(1);
  MgHierarchy.SetLevel(Finest,1);
  Needs.SaveAggregates = '1';
  status = MgHierarchy.FillHierarchy(PRfact, [], Acfact, 1, numDesiredLevels);
  MgHierarchy.SetSmoothers(GSFactory);

  NLevels=length(MgHierarchy.Levels_);
  fprintf('Levels Used = %d\n',NLevels);
  
  oc=status.OperatorComplexity;

  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Useful utility function from Chris' stash.
function s=cdate
v=datevec(date);
s=sprintf('%2.2d%2.2d%4d',v(2),v(3),v(1));
