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


function ImplicitSchur()




CGtol     = 1e-8;% tolerance for CG method
CGits     = 45;  % # of CG iterations to solve via AMG preconditioning
NumEminIts= 2;   % # of iterations of emin algorithm for generating prolongator
Nlevel    = 3;   % # of AMG levels


runAMG(CGtol,CGits,NumEminIts,Nlevel);
          
end

function [resvec,MgHierarchy] = runAMG(CGtol,CGits,NumEminIts,Nlevel)

        % retrieve XFEM info
        SetHomeDir
        mytests = { [MUEMAT_ROOT_DIR '/data/Xfem_1cr_al_prop3_090x91.mat']};
        load(mytests{1});

        fprintf('\n Explicit AMG :%30s \n=====================\n','data/Xfem_1cr_al_prop3_090x91.mat');

        % define how Afiltered is constructed. Afiltered is used for
        % the aggregation and the prolongator sparsity pattern associated
        % with the A_rr block.  To do this, the XFEM data structure
        % must be passed into the function TestNoCrossCrackDrop which is
        % used to remove edges in Arr that cross a crack.

        mydata.COORD        = node;
        mydata.xCr          = xCr;
        mydata.droptol      = 0.08;

        CoalesceDropFact    = CoalesceDropFactory();
        CoalesceDropFact.SetPreDropSpecifications(@TestNoCrossCrackDrop,mydata);

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Everything in this block is pretty standard AMG-Emin for elasticity.
        
        options.NCoarseDofPerNode = 3;
        options.PtentSA = 1;
        options.CoarseNullspaceWeightingMultiplier = [1,1,1];
        
        AggFact      = AggregationFactory();
        FCFact       = FCSplittingFactory(AggFact);

        NSFact       = SaCoarseNSFactory();
        PtentFact    = TentativePFactoryEx(CoalesceDropFact(), AggFact, NSFact, options);
        PatFact      = AP_PatternFactory([],PtentFact, options);
        PatFact.SetUseAfiltered(true);

        EFact=EminPFactory(PatFact, ConstraintFactory(PatFact,FCFact), CGEminSolver(NumEminIts), ...
                           PtentFact, options);

        %%%%%%%%%%%%%%%%%%  End of standard AMG-Emin block %%%%%%%%%%%%%%%%%%

        % Build a block diagonal prolongator = [P 0; 0 I] where I corresponds
        % to the special dofs. Emin is used to generate the (1,1) block
        % of the prolongator. Hiding inside this factory is something to use
        % Arr instead of the larger 2x2 system when applying emin.

        PFact = TestHybrid2x2PFactory(EFact,[]);

        GS     = Smoother('GaussSeidel',1,1.0);
        GS.SetDiagonalView('point');
        SFact = SmootherFactory(GS);

        mue_include;

        % Create a Level. In addition to the usual stuff like the matrix
        % and null space, also store Arr and the dof2node map. dof2node is 
        % used inside of the function 'TestNoCrossCrackDrop' and Arr is given to 
        % Emin to build prolongators. Both dof2node and Arr need to be 
        % projected on coarse levels. This is done by the RAPXfemFactory!

        Amat = Operator(A,1,1); % could do variable block?
        fr   = f;
        Finest = Level;
        Finest.KeepAll(false);
        Finest.Keep('Arr'); % needed for Hybrid2x2Factory and computation of subblock complexity (below)
        Finest.Set('A', Amat);
        Finest.Set('NullSpace', nsp);
        Finest.Set('Arr', Operator(A(rdof,rdof),2,2)); 
        Finest.Set('dof2node', dof2node);

        %%%%%%%%%%%%%%%%%%  standard AMG %%%%%%%%%%%%%%%%%%
        MgHierarchy = Hierarchy(); 
        MgHierarchy.SetLevel(Finest,1);
        status = MgHierarchy.FillHierarchy(GenericPRFactory(PFact),[],...
                                           RAPXfemFactory(FCFact),1,Nlevel); 
        MgHierarchy.SetSmoothers(SFact);
        MgHierarchy.SetOutputLevel(1);

        [x,flag,relres,iter,resvec]= pcg(A,fr,CGtol,CGits,...
           @(b)MgHierarchy.Iterate(b,1,zeros(length(fr),1),ALLZEROS,2)); 

        % compute sub-block complexity
        thesum = 0;
        for i=1:status.lastLevel
           nn   = MgHierarchy.Levels_{i}.Get('Arr').GetRowMap().NDOFs;
          temp = MgHierarchy.Levels_{i}.Get('A').GetMatrixData();
          thesum = thesum + nnz(temp) - nnz(temp(nn+1:end,nn+1:end));
          if i==1, A0 = thesum; end;
        end;
        fprintf('flag=%d; res = %7.4g after %d iterations.\nAMG Operator Complexity = %e, sub complexity = %e\n',flag,resvec(end),length(resvec)-1,status.OperatorComplexity,thesum/A0);

end

