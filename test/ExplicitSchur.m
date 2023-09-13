%
%   The following files used by this driver have closely related
%   equivalents in XMC. If major changes are made to either
%   the XMC or MueMat/test version, it might be wise to make suitable
%   updates so that both versions are roughly in sync.
%
%      XMC file name       MueMat file name           Difference 
%      ========================================================================
%      NoCrossCrackDrop.m  TestNoCrossCrackDrop.m   s/NoCrossCr/TestNoCrossCr/g
%      XfemLevel.m         TestXfemLevel.m          s/XfemLevel/TestXfemLevel/g
%      Hybrid2x2PFactory.m TestHybrid2x2PFactory.m  s/Hybrid2x2/TestHybrid2x2/g
%      SchurDrop.m         TestSchurDrop.m          s/SchurDrop/TestSchurDrop/g
%
%      BestAMG.m           ImplicitSchur.m          Abbreviated setup/load for
%                                                   ImplicitSchur, references
%                                                   to X/X1 data structure are
%                                                   removed so things like 
%                                                   X1.node are just node. Pre-
%                                                   append references to 
%                                                   NoCrossCrackDrop,XfemLevel
%                                                   and Hybrid2x2PFactory
%                                                   with 'Test'.
%
%      TestAMG.m           Implicit/ExplictSchur.m  A bit more substantial.

function ExplicitSchur()

CGtol = 1e-7;
CGtol = 1e-9;
CGits = 20;
NumEminSolverIts = 1;

    
runAMG(CGtol,CGits,NumEminSolverIts);
          
end

function [resvec,MgHierarchy] = runAMG(CGtol,CGits,NumEminSolverIts)

   % retrieve XFEM info
   SetHomeDir
   mytests = { [MUEMAT_ROOT_DIR '/data/Xfem_1cr_al_prop3_090x91.mat']};
   load(mytests{1});
   fprintf('\nSCHUR:%30s \n=====================\n','data/Xfem_1cr_al_prop3_090x91.mat');

        A = S;
        f = f(rdof);
        K11 = A(rdof,rdof);
        n11 = length(rdof);
        
        % The Operator
        Amat = Operator(A,2,2); 
        Finest = Level;
        Finest.KeepAll(true);
        Finest.Set('A', Amat);
        Finest.Set('NullSpace', nsp);
        Finest.Set('Arr', Operator(K11,2,2)); 
        CoalesceDropFact    = CoalesceDropFactory();
        mydata.dof2node     = dof2node;
        mydata.COORD        = node;
        mydata.xCr          = xCr;
        mydata.Schur        = 0;
        mydata.K11          = K11;
        mydata.droptol      = 0.001;

        CoalesceDropFact.SetPreDropSpecifications(@TestSchurDrop,mydata.droptol);
        Smoo               = Smoother('GaussSeidel',1,1.0);
        Smoo.SetDiagonalView('block');
        SFact = SmootherFactory(Smoo);

options.NCoarseDofPerNode = 3;
options.PtentSA = 1;
options.CoarseNullspaceWeightingMultiplier = [1,1,1];      
        
% aggregates and f/c coarsening and coarse grid nullspace
AggFact             = AggregationFactory();
NSFact              = SaCoarseNSFactory();

% initial prolongator
PtentFact           = TentativePFactoryEx(CoalesceDropFact, AggFact, NSFact, options);

% sparsity pattern constraint
PatNewfact             = AP_PatternFactory([],PtentFact,options);
PatNewfact.SetPatternType('AP');
PatNewfact.SetUseAfiltered(true);

% constraints (nullspace preservation)
ConstraintFact      = ConstraintFactory();

PFact=EminPFactory(PatNewfact, ConstraintFact, CGEminSolver(NumEminSolverIts), PtentFact, options);
PFact.SetEnergyMatrix('Arr');

mue_include;
MgHierarchy = Hierarchy(); 
MgHierarchy.SetLevel(Finest,1);
status      = MgHierarchy.FillHierarchy(GenericPRFactory(PFact),[],RAPXfemFactory(),1,3); 
MgHierarchy.SetSmoothers(SFact);
MgHierarchy.SetOutputLevel(1);

[x,flag,relres,iter,resvec]= pcg(A,f,CGtol,CGits,...
        @(b)MgHierarchy.Iterate(b,1,zeros(length(f),1),ALLZEROS,1));  %#ok<*ASGLU>
fprintf('flag = %d; res = %6.3g after %d iterations \n',flag,resvec(end),length(resvec)-1);
end

function Arr = LvlStaticGetArr(obj)
  Arr = obj.Get('Arr');
end