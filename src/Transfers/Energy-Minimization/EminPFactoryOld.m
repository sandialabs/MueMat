%TODO: InitialGuess_

classdef EminPFactory < PFactory
  % This factory creates a prolongator that is the solution of a constrained
  % energy minimization problem.

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Properties                                                                %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  properties (Access=private)
    SparsityFact_            = []
    CoarseNSFact_            = []
    Solver_                  = []
    SolverIterations_        = []
    AggFact_                 = []
    InitialGuess_            = []
    diagonalView_            = 'current' % diagonal view label (default == current view)
    options_                 = []
    ConstraintWgt_           = []
    GetEnergyMatrix_         = [];
    GetEnergyFunc_   = false;
    reUseP_          = false;
  end %properties (private)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Public methods                                                            %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  methods

    function [this] = EminPFactory(SparsityFact, CoarseNSFact, Solver, ...
                                 CoalesceFact, AggregationFact, Options, InitialGuess, diagonalView)
      %EMINPFACTORY Constructor
      %
      %   SYNTAX   obj = EminPFactory(SparsityFact, CoarseNSFact, Solver, CoalesceFact, AggregationFact, Options, InitialGuess);
      %
      %     SparsityFact    -  sparsity pattern factory
      %     CoarseNSFact    -  coarse nullspace factory
      %     Solver          -  solver handle for solving the 2x2 system
      %     CoalesceFact    -  coalesce (amalgation) factory
      %     AggregationFact -  aggregation factory
      %     Options         -  options structure
      %     InitialGuess    -  initial prolongator guess (optional)
      %     DiagonalView    -  view label

      this.SolverIterations_   = 5;

      if varexist('SparsityFact'), this.SparsityFact_ = SparsityFact;
      else this.SparsityFact_ = AP_PatternFactory(); this.SparsityFact_.SetPatternType('AP'); end;

      if varexist('CoarseNSFact'), this.CoarseNSFact_ = CoarseNSFact;
      else this.CoarseNSFact_ = SaCoarseNSFactory(); end;

      if varexist('Solver'), this.Solver_ = Solver;
      else this.Solver_ = @CGemin; end;

      if varexist('AggregationFact'), this.AggFact_ = AggregationFact;
      else this.AggFact_ = AggregationFactory(); end;

      if varexist('CoalesceFact') && ~isempty(CoalesceFact), this.AggFact_.SetCoalesceFactory(CoalesceFact);end

      if varexist('InitialGuess'), this.InitialGuess_ = InitialGuess; end;

      if varexist('diagonalView'), this.SetDiagonalView(diagonalView); end
      if varexist('Options'), this.options_ = Options; end;
      this.ConstraintWgt_      = 1;

      % internal logic: check if pattern factory knows its Factory
      % for the initial prolongator
      % if not: set the same InitPFact that is used for this PFactory
      if isempty(this.SparsityFact_.GetInitialPFactory())
          this.SparsityFact_.SetInitialPFactory(this.InitialGuess_);
      end

    end


    function [ToF] = SetGetEnergyMatrix(this, GetEnergyMatrixFunc)
       this.GetEnergyFunc_   = true;
       this.GetEnergyMatrix_ = GetEnergyMatrixFunc;
    end

    function [ToF] = ReUseP(this, ToF)
        if varexist('ToF'), this.reUseP_ = ToF; end
        ToF = this.reUseP_;
    end

    function [ToF] = SupportsRestrictionMode(this)
        ToF = false;
    end

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

    function [z] = GetNeeds(this)
      %GETNEEDS
      %
      %   SYNTAX   z = obj.GetNeeds();
      %
      %     z - Needs object
       z = GetNeeds@PFactory(this); % NeedsObject?

       if ~isempty(this.AggFact_),
         z =CrossFactory.MergeNeeds(z, this.AggFact_.GetNeeds());
       end
       if ~isempty(this.SparsityFact_),
         z =CrossFactory.MergeNeeds(z, this.SparsityFact_.GetNeeds());
       end
    end %GetNeeds()

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    function SetNeeds(this, FineLevel, CoarseLevel)
      % Obtain any cross factory specifications

       if ~isempty(this.AggFact_),
         this.AggFact_.SetNeeds(FineLevel);
       end
       if ~isempty(this.SparsityFact_),
         this.SparsityFact_.SetNeeds(FineLevel,CoarseLevel);
       end

       FineLevel.Request('Aggregates');
       FineLevel.Request('NullSpace');  % move me to initialP!!!

    end

    function SetConstraintWeight(this, alpha)
      %SETCONSTRAINTWEIGHT Sets a scalar weighting factor that scales all constraint equations
      % corresponding rows of P that have too few nonzeros to exactly interpolate
      % the fine nullspace.
      %
      %   SYNTAX   obj = obj.SetConstraintWeight(alpha);
      %
      %     alpha -
      %     obj   -
      this.ConstraintWgt_ = alpha;
    end

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

     function SetDiagonalView(this, diagonalView)
     % indicates use of either point or block diagonal prolongator smoothing.
        %this.diagonalView_ = diagonalView;
      warning('SetDiagonalView() has no effect right now');
     end

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

    function SetSparsityFactory(this,SparsityFact)
      %SETSPARSITYFACTORY Set the factory that will generate the prolongator sparsity pattern.
      %
      %   SYNTAX   obj = obj.SetSparsityFactory(SparsityFact);
      %
      %     SparsityFact -
      %     obj          -
      this.SparsityFact_ = SparsityFact;
    end

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

    function SetCoarseNSFactory(this,CoarseNSFact)
      %SETCOARSENSFACTORY Set the factory that will generate the coarse grid nullspace.
      %
      %   SYNTAX   obj = obj.SetCoarseNSFactory(CoarseNSFact);
      %
      %     CoarseNSFact -
      %     obj          -
      this.CoarseNSFact_ = CoarseNSFact;
    end

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

    function SetSolver(this,Solver)
      %SETSOLVER Set the function handle of the solver that will solve the 2x2 optimization problem.
      %
      %   SYNTAX   obj = obj.SetSolver(Solver);
      %
      %     Solver -
      %     obj    -
      this.Solver_             = Solver;
    end

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

    function SetInitialGuess(this,InitialGuess)
      %SETINITIALGUESS
      %
      %   SYNTAX   obj = obj.SetInitialGuess(InitialGuess);
      %
      %     InitialGuess -
      %     obj          -
      this.InitialGuess_       = InitialGuess;
    end

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

    function SetSolverIterations(this,nits)
      %SETSOLVERITERATIONS Specify number of iterations the solver should use in solving the 2x2 system corresponding
      % to the optimization problem.
      % Defaults to 10.
      %
      %   SYNTAX   obj = obj.SetSolverIterations(nits);
      %
      %     nits -
      %     obj  -
      this.SolverIterations_ = nits;
    end

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

    function flag = BuildP(this,FineLevel,CoarseLevel)
      %BUILD Solve the optimization problem to produce an energy-minimized prolongator.
      %
      %   SYNTAX   [FineLevel, CoarseLevel, flag] = obj.Build(FineLevel, CoarseLevel);
      %
      %     FineLevel   -
      %     CoarseLevel -
      %     flag        -
      flag = true;
      %if ~isempty('this.maxCoarseSize_') && FineLevel.Get('A').GetRowMap().NNodes() <= this.MaxCoarseSize_
      %  flag = false;
      %  return;
      %end
      % rst: Checks to see if data is available. Seems a little strange here.
      % rst: It is not real clear that ReUseP() fully makes sense. It would
      % rst: be better to have
      % rst:       this.GetInitialGuess() which points to Ptent, P, or perhaps
      % rst:                              something from the Afc world.
      % rst:                              The main trick is that if this lives
      % rst:                              in the EminPFactory class, how can
      % rst:                              user's easily define their own data
      % rst:                              that this can point to ... which must
      % rst:                              work on each level in the hierarchy.
      % rst:                              Should this be some string which is
      % rst:                              searched in the big hash table.
      % rst: Also, there is no explicit need to have a coarse null space at
      % rst: this point. In general, the Fine/Coarse Null space serve only one
      % rst: purpose and that is to define constraints. They can be optionally
      % rst: used to generate initial guess.  Thus, I would prefer that there
      % rst: is no reference to NullSpace inside this class.
      % rst:
      % rst:       this.GetCoarseNull()   could point to an existing coarse
      % rst:                              null space otherwise there needs
      % rst:                              to be an algorithm to compute it.
      % rst: Suggestion: remove this
      if this.ReUseP()
          if ~CoarseLevel.IsAvailable('P') || ~CoarseLevel.IsAvailable('NullSpace'),
              warning ('NoReuse', ' P or NewNull not found');
          end
          P = CoarseLevel.Get('P');
          NewNull =  CoarseLevel.Get('NullSpace');
          CoarseLevel.Release('NullSpace');
          if ~isempty(P) && ~isempty(NewNull), return; end
          warning('NoReuse','P or NewNull not found');
      end

      if ~FineLevel.IsAvailable('A')
          fprintf('EminFactory.Build:: Amat not found\n');
          keyboard;
      end
      Amat    = FineLevel.Get('A');
      % rst: Suggestion: move this to SparsityPattern (if Ptent is needed) and to
      % rst: BuildConstraints if not already done and we have at least some
      % rst: constraints.
      if ~FineLevel.IsAvailable('NullSpace')
          NS = BuildNullSpace(Amat);
          FineLevel.Set('NullSpace', NS);
      end
      NS      = FineLevel.Get('NullSpace');
      FineLevel.Release('NullSpace');

      P = []; cnull = [];

      % rst: Aggregates only needed to generate sparsity patterns. There should
      % rst: be no reference to aggregation here.  It looks like we are adding
      % rst: an option to save the graph that is used to build the aggregates
      % rst: Not sure why this is needed?
      if isempty(P) || isempty(cnull),
        AggregationPtent = 1;
        if AggregationPtent,
          opts.SaveGraph='1'; % 1 ??
%DEPRECATED

           % Build aggregates
           this.AggFact_.Build(FineLevel);
           AggInfo = FineLevel.Get('Aggregates');
           FineLevel.Release('Aggregates');

          % assert
          if ~FineLevel.IsAvailable('Graph'), keyboard; end;
       % rst: Really not sure why this fine level null space stuff is here.
       % rst: It would be better if this gets buried somewhere in the thing that
       % rst: builds the constraints.
       % rst: Suggestion: Same as line 239 suggestion.
          % Jacob's nullspace massage to improve energy minimization
          % near the boundaries
          if isfield(this.options_,'SmoothNullspace'),
            fprintf('EminFactory: Smoothing nullspace before energy minimization\n');
            Adata=Amat.GetMatrixData();
            D=spdiags(diag(Adata),0,size(Adata,1),size(Adata,2));
            DA=D\Adata;
            NSold=NS;
            lambda_max = Amat.GetDinvALambda(this.diagonalView_);

            for iii=1:this.options_.SmoothNullspace
               NS=NS - 4/(3*lambda_max)*(DA*NS);
            end;
            NS=normalize(NS);
          end

          % Build Ptent and coarse nullspace
   %rst: what is the testPerm option? It looks like it is basically
   %rst: commented out. No ... my mistake ... the if clause is the
   %rst: normal case where we produce a Ptent and a CNull . We might be doing
   %rst: some strange stuff to compress out the null space. Actually, it looks
   %rst: cnull is always injection of fnull. how does this work for elasticity
   %rst: oh ... I see it has to do with local or global svd ....
          if ~isfield(this.options_,'testPerm')
            coords=[FineLevel.Get2('xcoords'),FineLevel.Get2('ycoords'),FineLevel.Get2('xcoords')];
            [P,cnull] = this.CoarseNSFact_.Build(AggInfo, Amat, NS, this.options_,this.GetOutputLevel(),[],coords); %AggInfo/AggInfo

            % Switch around the nullspace signs for poor
            % convergence.  Question: Why?
            %cnull=[cnull(:,1),cnull(:,2),-cnull(:,3)];
          else
% rst: the else clause looks like it just always puts translations first
% rst: seems a little strange for 3D case as we would have 3 translations
            % Test of the effect initial guess vs. coarse
            % nullspace : hack code to use 2 translation
            % nullspace with translation/rotation Ptent
            % (set opt_global_svd=1)
            if (version('-release') == '2010a')
              %R2010a
              rot=1; % nullspace vectors order after global svd...
              tr1=2;
              tr2=3;
            else
              keyboard;
            end

            permtt=[tr1 tr2 rot];
            permrt=[rot tr1 tr2];
            % permrt=[rot tr2 tr1];

            [Ptt,cnulltt] = this.CoarseNSFact_(AggInfo, Amat, NS, this.options_,this.GetOutputLevel(), permtt);
            [Prt,cnullrt] = this.CoarseNSFact_(AggInfo, Amat, NS, this.options_,this.GetOutputLevel(), permrt);
            P=Ptt; cnull=cnulltt;
            % P=Prt; cnull=cnullrt;
            % P=Prt; cnull=cnulltt;
            % P=Ptt; cnull=cnullrt;
          end
        else
          %  A Simple F/C approach:
          %      1) Get root points
          %      2) Take Ptent to be zero except@CPts (its just an initial guess)
          %      3) Later define sparsity pattern via Schur complement dropping
          %          CptRows= EminFactory.BuildCoarseRows(FineLevel,CoarseLevel,...
          %                               this.CoalesceFact_, this.FCFact_,...
          %                               this.GetOutputLevel());
          %           cnull = NS(CptRows,:);
          %           P = sparse( Amat.RowMap.NDOFs(), length(CptRows) );
          %           P(CptRows,:) = speye( length(CptRows), length(CptRows));
        end
      end

% rst: Ptent should not appear in the function
% DEPRECATED
%      if Specs.TrueOrFalse('SavePtent',CoarseLevel.GetLevelId()),
%        CoarseLevel.Set('Ptent', P);
        Ptent = P.Copy();
        CoarseLevel.Set('Ptent',Ptent);
%      end

      % Get Ncpn
      if ~isfield(this.options_,'NCoarseDofPerNode') Ncpn = size(NS,2);
      else Ncpn = this.options_.NCoarseDofPerNode;end;

      % Pull matdata
      Amatrixdata = Amat.GetMatrixData();
      Pmatrixdata = P.GetMatrixData();

% rst: This looks like it should go in the sparsity pattern factory.
% rst: It basically takes the graph, de-coalesces, multiplies that
% rst: with Ptent and chooses that as the sparsity pattern.
% rst: By the way, this would best be done by having a GetAForSparsity().
       if(isfield(this.options_,'FilteredAP') && this.options_.FilteredAP==true...
          && FineLevel.IsAvailable('Graph'))
        % Use the Filtered A to build the AP pattern
        % The filtered A should come from the aggregation graph
        % from BuildAggregates, which is why we saved it.
        %
        % This means that we need a *nodal* Ptent, which we
        % should probably cache from somewhere else, but I'm lazy
        % at the moment
        fprintf('Emin: Using Filtered(A)*P for initial pattern\n');

        % Build the nodal pattern
        Graph=FineLevel.Get('Graph');
        CFact = CoalesceDropFactory();
        CFact.SetPostDropSpecifications([],[]);
        CFact.SetPreDropSpecifications([],[]);
        Pnode=CFact.Build_(P);
        ftemp = FineLevel.BuildMe();
        ftemp.Set('A', Graph);
        ctemp = CoarseLevel.BuildMe();
%        ctemp.Set('Ptent', Pnode);
        ctemp.Set('Ptent',Pnode);
        ctemp.Request('Ptent');
        ctemp.Request('Ptent');
        ctemp.Request('Ptent');
        ctemp.Request('Ptent');
        APnode=this.SparsityFact_.Build(ftemp, ctemp);

        % Expand back to dof-sized pattern
        %SparsityPattern=sparse(size(P.GetMatrixData(),1),size(P.GetMatrixData(),2));
        %for ii=1:Ncpn,
        %  for jj=1:Ncpn,
        %    SparsityPattern(ii:Ncpn:end,jj:Ncpn:end)=APnode;
        %  end
        %end
        SparsityPattern=this.DecoalesceMatrix(APnode,Ncpn,Ncpn);
      else
        % Get SparsityPattern the obvious way
%        CoarseLevel.Set('Ptent', P);
         if ~CoarseLevel.IsAvailable('Ptent')
            Ptent = P.Copy();
            CoarseLevel.Set('Ptent',Ptent);
         end
        SparsityPattern = this.SparsityFact_.Build(FineLevel, CoarseLevel);
%        CoarseLevel.Set('Ptent', []);   % not sure if we really want to store ptent
%        CoarseLevel.Save('Ptent',[]);
      end

      % OPTION: Smooth with pattern of A^2 P
      if(isfield(this.options_,'MultiSmoothing')),
        SparsityPattern=double((abs(Amatrixdata)>0)*SparsityPattern>0);
      end

      % else
      % % A Simple Schur Complement Sparsity Pattern
      %       NRows = P.RowMap.NDOFs();
      %       NCols = P.ColMap.NDOFs();
      %       temp  = ones(NRows,1); temp(CptRows) = 0;
      %       FptRows = find(temp);
      %       Aff = Amat.MatrixData(FptRows,FptRows)
      %       Afc = Amat.MatrixData(FptRows,CptRows)
      %       Pattern = sparse(NRows,NCols);
      %       Pattern(CptRows,,:) = speye(NCols,NCols);
      %       temp = Aff\Afc;
      %       Pattern(FptRows,:) = abs(temp) > .1;
      % end

      % Tests with optimal prolongator
% rst: I suspect that all of this within the if can be thrown out
      if (1==0)

        NRows = P.GetRowMap.NDOFs();
        NCols = P.GetColMap.NDOFs();

        CptRows = AggInfo.Roots;
        FptRows = setxor(1:NRows, CptRows);

        % Acc = Amatrixdata(CptRows,CptRows); %
        % Acf = Amatrixdata(CptRows,FptRows); %

        Aff = Amatrixdata(FptRows,FptRows);
        Afc = Amatrixdata(FptRows,CptRows);

        Pc = speye(NCols,NCols);
        Pf = -(Aff\Afc);

        OptP = sparse(NRows,NCols);
        OptP(CptRows,:) = Pc;
        OptP(FptRows,:) = Pf;

        % Schur = Acc+Acf*Pf;
        % full(Schur - OptP'*Amatrixdata*OptP)

        % TESTS:
        % SparsityPattern = spones(OptP); % use the pattern of the optimal prolongator
        % cnull(:)=1;               % force cnull=1 1 1 1
      end %if (1==0)

      residual = NS - P*cnull;

% rst: This can be put in BuildConstraints(). It accomplishes something
% rst: similar to Jacob's massage and so the two codes should be together.
% rst: The issue is that the near null space is not accurate near boundaries
% rst: as it corresponds to a true null space for the boundary-free PDE.
% rst: We might compensate by smoothing it or not satisfying near null space
% rst: constraints near boundaries. The trick is detecting boundaries. Option 1
% rst: checks how close the accuracy of the near null space locally. This
% rst: determines constraints to ignore.  Unfortunately, the null space is
% rst: inaccurate everywhere for transient problems due to the Mass. Option 2
% rst  looks for prolongator rows with a small number of nonzeros compared to
% rst: the number of null space modes.  These are assumed to be near boundaries
% rst: and so ignored (though interior rows could fall into this category if
% rst: there is not enough overlap between basis functions). We should really
% rst: check (for each prolongator row) whether we can find a P such that
% rst:   FineNull(row,:) = P(row,:)*CoarseNull(dofs,:) where dofs correspond
% rst: to nonzero locations within the row.  This is equivalent to
% rst:   FineNull(row,:)' = CoarseNull(dofs,:)'*P(row,:)'
% rst: which is a least-squares problem. If dofs is small and the number of
% rst: Null Space vectors is high, we have
% rst:       [ FineNull_1 ]  =   [ Coarse_11 Coarse_12 ]
% rst:       [ FineNull_2 ]  =   [ Coarse_21 Coarse_22 ] Prow'
% rst:       [ FineNull_3 ]  =   [ Coarse_31 Coarse_32 ]
% rst: and this might not be have a solution. Here, FineNull_i gives the
% rst: ith null space vector at row and Coarse_ij gives the ith coarse
% rst: null vector evaluated at the j^th nonzero location. It cannot be
% rst: solved if [FineNull_i] has a component in the left nullspace of Coarse_ij
% rst: We don't examine this, but do check for root point rows (which can be
% rst: solved when FineNull is injected to get CoarseNull.
      % Detection methods for Mandel-esque dropping
      DropWeight=zeros(size(SparsityPattern,1),1);
      if(isfield(this.options_,'DropConstraintsMethod')),
        % Ray-style dropping detection
        if(strcmp(this.options_.DropConstraintsMethod,'nullspace')),
          fprintf('Constraint Dropping: Nullspace\n');
          nn = Amat.GetRowMap().NDOFs();
          temp = spdiags( 1./diag(Amatrixdata),[0],nn,nn)*Amatrixdata*NS;
          temp = temp*spdiags( 1./max(NS)',[0],size(NS,2),size(NS,2));
          filter = abs(temp) > 1e-4;
          filter = filter*ones(size(NS,2),1);
          DropWeight(find(filter))=1;
          % make sure that all the root points are still constrained
          rrows = Node2DOF(FineLevel.Get('Aggregates').Roots,Amat.GetRowMap);
          DropWeight(rrows) = 0;
        elseif(strcmp(this.options_.DropConstraintsMethod,'skinny')),
          % Skinny non-root style dropping
          fprintf('Constraint Dropping: Skinny non-root\n');
          s=sum(abs(SparsityPattern)>0,2);
          isRoot = zeros(size(SparsityPattern,1),1);
          for jj=1:length(AggInfo.Roots)
            RIDX = (AggInfo.Roots(jj)-1)*Ncpn + [1:Ncpn];
            isRoot(RIDX) = 1;
          end
          DropWeight=(s<size(NS,2) & ~isRoot);

        elseif(strcmp(this.options_.DropConstraintsMethod,'fpoints'))
          fprintf('Constraint Dropping: Fpoints\n');
          DropWeight(1:end) = 1;
          rrows = Node2DOF(FineLevel.GetAggregates.Roots,Amat.GetRowMap);
          DropWeight(rrows) = 0;
        end
      end


      [B,ConstraintRhs]=EminPFactory.BuildConstraints(SparsityPattern,cnull,residual,this.ConstraintWgt_,DropWeight);

% for visualizing what gets weighted -- you also have to uncomment stuff in BuildConstraints
% also, dimensions are hardwired for the 120^2 problem
%global gwgt
%dbstack,keyboard
%figure; plot3(coords(:,1), coords(:,2), gwgt(1:2:28800),'r.'); box on; grid on;
%figure; plot3(coords(:,1), coords(:,2), gwgt(2:2:28800),'b.'); box on; grid on;
%figure; plot3(coords(:,1), coords(:,2), gwgt(28801:2:57600),'r.'); box on; grid on;
%figure; plot3(coords(:,1), coords(:,2), gwgt(28802:2:57600),'b.'); box on; grid on;
%figure; plot3(coords(:,1), coords(:,2), gwgt(57601:2:86400),'r.'); box on; grid on;
%figure; plot3(coords(:,1), coords(:,2), gwgt(57602:2:86400),'b.'); box on; grid on;

%rst: This is a modification of the sparsity pattern to essentially
%rst: put an I at the root nodes. We also remove constraints
%rst: that might be associated with these root point nodes and I believe
%rst: add constraints to force an I at the root nodes.
%rst: It seems like most of this should go in the sparsity pattern code.
%rst: It also seems that something should be flagged for the build
%rst: constraints. We might want to have something like ... if there
%rst: is only 1 nnz in a P row, we should either inject the fine null space
%rst: and only build one constraint which says that P_ij = 1 or if we
%rst: already have a fine and coarse null space ... P_ij = F_i/C_j
   % Ray's Secret Document, Item #4a:
   % Add constraints to force P = Ptent @ root nodes.
      if (isfield(this.options_,'PtentRootModifications')),
        if(strcmp(this.options_.PtentRootModifications,'4a'))
          fprintf('Running Ray''s Special Sauce, Option #4a\n');
       %
       %      (  P1  P2          )
       %      (  P3  P4          )
       %  P = (                  )
       %      (          P5  P6  )
       %      (          P7  P8  )
       %      (          P9  P10 )
       %
       % if P(2,:) corresponds to a root node:
       %
       %        P1   P2   P3  P4  ...
       %        |    |    |    |
       %      ( ..   ..   |    |  ... )
       %      (                |      )
       %  B = (           1    |      )
       %      (                       )
       %      (                1      )
       %      (                       )

       % 1) Finding the columns of B corresponding to root nodes.
       %    - There is one column in B for each nz of SparsityPattern
       %    - I mark root nodes in SparsityPattern (value -1)
       %    - nz of SparsityPattern are sorted as in B and I get the
          %      columns indices of root nodes in B
       SparsityPatternTmp = SparsityPattern;
       NAggregates = max(AggInfo.AggId);
       %mnulldim = this.options_.NCoarseDofPerNode;
       mnulldim = size(P.GetMatrixData(),2)/NAggregates;
       for i=1:NAggregates
         rnode = AggInfo.Roots(i);
         rows  = Node2DOF(rnode,P.GetRowMap());
         CDofs = ((i-1)*mnulldim+1:i*mnulldim);
         SparsityPatternTmp(rows,CDofs) = -1; % in SparsityPattern, we distinguish the P_i
                                        % corresponding to root nodes
       end
       [I,J,V] = find(SparsityPatternTmp);
       [I,perm] = sort(I); clear J; V = V(perm); % sort: same order as in BuildConstraints

       col = find(V==-1); % list of root node columns of B
       % col=[]

       % 2) Remove constraints involving root nodes.
       [I,J] = find(B(:,col));         % I: set of row in B involving root nodes
       toKeep = setxor(1:size(B,1),I); % toKeep: complement of the set I
       B = B(toKeep,:);                % Update B and ConstraintsRhs
       ConstraintRhs = ConstraintRhs(toKeep);

       % 3) Add new constraints for root nodes (at the end of B)
       nrowB = size(B,1);
       [I,J,V] = find(B);
       I = cat(1, I, (nrowB+1:nrowB+size(col))');
       J = cat(1, J, col);
       V = cat(1, V, repmat(1,size(col),1));
       B = sparse(I,J,V); clear I J V;
       ConstraintRhs = cat(1, ConstraintRhs, repmat(0,size(col),1));

     end % #4a
   end % #4a
   Ptentdata = Pmatrixdata; % bkp Ptent to test #4c afterward

      BBt = B * B';
      dBBt = diag(diag(BBt));
      fprintf('(before scaling) condest(BBt) = %g\n',condest(BBt,5));
      BBt = dBBt \ BBt;
      ConstraintRhs = dBBt \ ConstraintRhs;

% rst: This looks like we are making sure that the initial guess satisfies
% rst: the constraints. There should probably be an option to skip this.
% rst: if we are sure that the initial guess already satisfies the constraints.
      % Make sure P has the right sparsity pattern.
      % Then, force P NScoarse = NSfine (corresponding
      % to the update P = (NSfine-P*NScoarse)*inv(NScoarse'NScoarse)*NScoarse
      % (JG: useless here but I keep it in the case we get the
      % SparsityPattern and Ptent from the user)
      Pmatrixdata = Pmatrixdata .* SparsityPattern;

      % Satisfy the constraints (at least in a least-squares sense)
      % PInds is used by Flatten() and BlowUp() to
      % convert vector to sparse matrices and back
      % to vectors.
      %
      PInds = find(reshape(SparsityPattern',[],1));
      [nrows,ncols] = size(Pmatrixdata);
      Pmatrixdata =  Pmatrixdata + BlowUp(B'*(BBt\ConstraintRhs),nrows,ncols,PInds);
      dif = Pmatrixdata*cnull-NS;
      fprintf('column by column norm(FNull - P(noenergy)*CNull) = ');
      for i=1:size(NS,2), fprintf('%10.2e ',norm(dif(:,i))); end
      fprintf('\n');
      P.SetMatrixData(Pmatrixdata);

      % Diagnostic
      Pcons=P;

      % Energy Minimization

      if this.GetEnergyFunc_,
         temp = this.GetEnergyMatrix_(FineLevel);
         EminAmatrixdata = temp.GetMatrixData();
      else
         EminAmatrixdata = Amatrixdata;
      end

      if isa(this.Solver_,'function_handle') && this.SolverIterations_ > 0
        newP = this.Solver_(EminAmatrixdata, SparsityPattern, B, Pmatrixdata, this.SolverIterations_);
        P = Operator(newP, P.GetRowMap(), P.GetColMap(), P.GetApply());

% rst: I think we can just remove !!!
       % Test if Ray's Secret Document, Item #4a is working
       if (isfield(this.options_,'PtentRootModifications')),
           if(strcmp(this.options_.PtentRootModifications,'4a'))
               Test4a = this.CompareRootNodeParts(AggInfo, Ptentdata, newP, P.GetRowMap(),this.options_); % Test if Pfinal = Ptent @ root nodes
           end % #4a
       end % #4a
      end

      % P = Operator(OptP, P.GetRowMap(), P.GetColMap(), P.GetApply()); % use optimal prolongator

% rst: remove this
      opt_debug = 0;
      if (opt_debug)
        disp('P ='); disp(full(Pmatrixdata));
      end
% rst: These should be diagonally scaled.
      fprintf('condest(Pfinal''Pfinal) = %10.2e\n', condest(P.GetMatrixData()'*P.GetMatrixData()));
      fprintf('condest(P''AP) = %10.2e\n', condest(P.GetMatrixData()'*Amatrixdata*P.GetMatrixData()));
% rst: this should probably invoke something in the near null constraint class which
% rst: indicates how well the constraints are satisfied.
      fprintf('norm(P.cnull-fnull)=%f, nnzP = %d\n', norm(P*cnull-NS,inf), nnz(P.GetMatrixData()));
      dif = P*cnull-NS;
      fprintf('column by column norm(FNull - P(energy)*CNull) = ');
      for jj=1:size(NS,2), fprintf('%10.2e ',norm(dif(:,jj))); end
      fprintf('\n');

% rst: there is some scaling but somehow it seems that we should do some row-wise
% rst: scaling. (Note: it might be nice to use something other than diagonal scaling
% rst: within the energy minimization algorithms .. which reappears here). That is,
% rst:     trace(diag(PAP)/diag(P'*P))
% rst:
% rst: Note: Pcons is the P matrix before energy minimization
      D=diag(diag(Amatrixdata));
      PAP=P.GetMatrixData()'*(D\(Amatrixdata*P.GetMatrixData()));
      PtAPt=Ptentdata'*(D\(Amatrixdata*Ptentdata));
      PcAPc=Pcons.GetMatrixData()'*(D\(Amatrixdata*Pcons.GetMatrixData()));

      if ~isfield(this.options_,'NCoarseDofPerNode') Ncpn = size(NS,2);
      else Ncpn = this.options_.NCoarseDofPerNode;end

      % Total energy diagnostics
      energy=trace(PAP);
      tent_energy=trace(PtAPt);
      cons_energy=trace(PcAPc);
      fprintf('Initial energy (sum_i pt_i^T A pt_i)     = %6.4e\n',full(tent_energy));
      fprintf('Constrained energy (sum_i pc_i^T A pc_i) = %6.4e\n',full(cons_energy));
      fprintf('Final energy (sum_i p_i^T A p_i)        = %6.4e\n',full(energy));

% rst: local orthogonalization ... I wonder if we could stick this outside of Emin.
% rst: We could perhaps do something where we have a PostProc factory which first
% rst: invokes another factory ... in this case 'emin' and then does some post-processing
      % Ray's Secret Document, Item #3: Force linear independance
      % of P columns after energy minimization using an
      % aggregate-wise QR.
      %
      if (isfield(this.options_,'PtentRootModifications')),
        if(strcmp(this.options_.PtentRootModifications,'3')) % /!\ PtentRootModifications = string
          fprintf('Running Ray''s Special Sauce, Option #3\n');

          [P, cnull] = aggregateWiseQR(AggInfo, P, cnull);
          fprintf('condest(Pfinal''Pfinal) = %10.2e\n', condest(P.GetMatrixData()'*P.GetMatrixData()));
          fprintf('norm(P.cnull-fnull)=%f, nnzP = %d\n', norm(full(P*cnull-NS)), nnz(P.GetMatrixData()));
        end
      end

      CoarseLevel.Set('P', P);
      CoarseLevel.Set('NullSpace', cnull);

      if isa(this.SparsityFact_,'myPatternFactory')
          ArrNew = P'*EminAmatrixdata*P;
          ArrNew = ArrNew.GetMatrixData();
          this.options_.Arr = ArrNew;
          this.SparsityFact_.ArrNew_ = ArrNew;
      end

    end %Build()

   end % public methods
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

   methods (Static = true)
    function [B, ConstraintRhs] = BuildConstraints(Ppattern, cnull, rhs, wgtFactor, toWeight)
      %BUILDCONSTRAINTS Build the constrains matrix from the pattern of P and the coarse nullspace vectors.
      %
      %   SYNTAX   [B, ConstraintRhs] = obj.BuildConstraints(Ppattern, cnull, rhs, wgtFactor, toWeight);
      %
      %     Ppattern      - pattern of the prolongator
      %     cnull         - coarse nullspace vectors
      %     rhs           - right hand side
      %     wgtFactor     - scalar that multiplies too-skinny constraint equations
      %     toWeight      - A vector w/ zeros for equations that don't get weighted and ones for ones that do.
      %     B             - matrix of constraints
      %     constraintRhs - corresponding right hand side, reshaped and perhaps with element eliminated
      % (TODO: change the name of constraintRhs?)

      nulldim = size(cnull,2);
      finenulldim = size(rhs,2);

      % options: debug
      opt_debug = 0;

      % options: 2 methods to remove rows from B
      opt_local_svd  = 1; % default = 1
      opt_global_svd = 0; % default = 0

      if (opt_debug)
        disp('Ppattern ='); disp(''); disp(full(Ppattern));
        nulldim
      % disp('nullspace ='); disp(''); disp(full(nullspace));
        disp('cnull ='); disp(''); disp(full(cnull));
      end

      % Build B (using Ppattern and cnull)
      %
      % For each k in nullspace(:,k), fill in the values of
      % P with nullspace values at the same time shift
      % compress the nonzero matrix columns in each row and
      % shift the matrix columns to the right of the previous
      % rows matrix columns. This shifting is to mimic the
      % fact that each nonzero in the sparsity pattern is
      % actually its own unique unknown. The resulting
      % matrix will be m times as tall as P where
      % m is the number of nullspace vectors. The width of
      % the matrix B is nnz(Ppattern).

      % Get Ppattern in IJV format and sorted by row indices.
      [I,J] = find(Ppattern);
      [I,perm] = sort(I); J = J(perm);

      nnzP  = length(I);        % = nnz(Ppattern)
      nrowP = size(Ppattern,1); % = max(I)

      % Build B : B=(II,JJ,VV) with
      %  II = permutation([I 2I 3I ... (nulldim-1)*I])
      %  JJ = [J  J  J ... J]
      %  VV = nullspace values (if B(i,j)~=0, B(i,j)=cnull(j))
      II = zeros(nnzP*nulldim,1);
      JJ = zeros(nnzP*nulldim,1);
      for i=1:nulldim
        II((i-1)*nnzP+1:i*nnzP) = (I-1)*nulldim + i; % permute B so that rows involving same dofs are adjacent
        JJ((i-1)*nnzP+1:i*nnzP) = (1:nnzP); % so nnzP == size(B,2)
      end

      VV = reshape(double(full(cnull(J,:))), [],1);

      B = sparse(II,JJ,VV); clear I J II JJ VV;

      ConstraintRhs = reshape(rhs',[],1);

      if (opt_debug)
        disp('B (initial) ='); disp(''); disp(full(B));
      end

      % Remove rows from B:
      % Remove constraints which are identical to each other.
      % This can happen, for example, if we have two identical
      % nullspace vectors. It can also happen in more subtle
      % situations. Several methods for removing rows are implemented.

      % Method 1 : Local SVD (default)
      % Linear dependencies between rows in B only appears between rows
      % involving the same unknowns.
      %
      % For each row of Ppattern, the number of corresponding row in B is
      % nulldim. Between each row of Ppattern, there is a column-shift in
      % B. So, we only have to test linear dependencies of rows of B that
      % correspond to the  same line of Ppattern. Thanks to the previous
      % permutation, such rows are adjacent.
      % To resume, we extract and test the rank of subsystems:
      % ( aaaaaaa )  ( bbbbbbb )   ( ccccccc )
      % ( aaaaaaa )  ( bbbbbbb )   ( ccccccc )

      tol = 1e-8;
      fprintf('BuildConstraints: dropping rows of B w/ singular value < %d\n',tol);
      if (opt_local_svd)
        %
        icol=1;  % current col in B (shift)
        irow=1;  % current row in newB
        newnz=1; % nnz of newB (nnz-1)

        if size(Ppattern,2) == 1,  %handles an extreme case
           temp = Ppattern; temp(:,2) = 0; nnzs = sum(temp');
        else
           nnzs = sum(Ppattern'); % if we use nnz(Ppattern(i,:)) inside the
                                  % loop, it is too slow so we use
                                  % sum() to store nnz for each row.
        end

        % preallocate memory for newB
        I = ones(nnzP*nulldim,1);
        J = ones(nnzP*nulldim,1);
        V = zeros(nnzP*nulldim,1);
        % preallocate memory for rhs
        newConstraintRhs = zeros(length(ConstraintRhs),1);

        % We have to apply local SVD to the submatrix
        % M = B(rowstart:rowend,icol:icol+s-1);
        % but as this matlab operation is very slow, we use IJV/CSR format ...
        %
        % B in IJV format, sorted by rows
        [iB,jB,vB] = find(B);
        [iB,perm] = sort(iB); jB = jB(perm); vB = vB(perm);
        %
        % B in CSR format (row_ptr, jB, vB)
        % row_ptr vector stores the locations in the val vector 'vB' that start a row
        offset = zeros(size(B,1),1); % offset(i) = number of nz in row 'i'
        for i=1:length(iB), offset(iB(i)) = offset(iB(i))+1; end
        row_ptr = cat(1, 1, cumsum(offset)+1);

        %wgtVector = 100*ones(length(ConstraintRhs),1); %FIXME for visualization
        wgtVector = ones(length(ConstraintRhs),1);
        numUnderWeighted=0;
        %global gwgt
        %gwgt = 100*ones(length(ConstraintRhs),1);   % for later plotting purposes
        %gstep = 1;

        % Main loop
        for ii=1:size(Ppattern,1)
          s=nnzs(ii);           % number of unknowns for the current row in
                                % B ( == nnz(Ppattern(ii,:)) ).

          if s > 0,
             % Extract a submatrix M of B
             rowstart = nulldim*(ii-1)+1;
             rowend = nulldim*ii;
             % M = B(rowstart:rowend,icol:icol+s-1);         % -> slow
             ind    = row_ptr(rowstart):row_ptr(rowend+1)-1; % faster
             iBpart = iB(ind) - (rowstart-1);                %
             jBpart = jB(ind) - (icol-1);                    %
             vBpart = vB(ind);                               %
             M = sparse(iBpart, jBpart, vBpart,rowend-rowstart+1,s);

             [U,S,dummy]=svd(full(M));    % SVD on the submatrix
             z = abs(S);
             % Not sure what is right. Do we want to retain largest singular
             % values even if all are small (meaning a little submatrix in
             % BBt is small).  Consider, for example, scaling the null space
             % by 10^-5. Then, the constraint matrix is small, but we want
             % to keep it. On the other hand, consider a null space with SOME
             % small entries within a region. Interpolation in this region
             % should be unconstrained (because small coarse null space values
             % interpolate to small values without constraint).

             k = nnz(z > tol*max(max(abs(z))));       % k : rank of M

             % Deduce newB and newConstraintRhs
             M = U(:,1:k)' * M;
             % newB(irow:irow+k-1,icol:icol+s-1) = M;         % -> slow
             [iM,jM,vM] = find(M);                            % faster
             iMleng = length(iM);                             %
             I(newnz:newnz+iMleng-1) = iM(1:iMleng) + irow-1; %
             J(newnz:newnz+iMleng-1) = jM(1:iMleng) + icol-1; %
             V(newnz:newnz+iMleng-1) = vM(1:iMleng);          %
             newnz = newnz + iMleng;                          %

             newConstraintRhs(irow:irow+k-1) = U(:,1:k)' * ConstraintRhs(rowstart:rowend);

             if(toWeight(ii)),
               wgtVector(irow:irow+k-1) = wgtVector(irow:irow+k-1) * wgtFactor;
               numUnderWeighted = numUnderWeighted + k;
               %gwgt(gstep  : gstep+k-1)          = wgtVector(irow:irow+k-1);
             end
             %if k < size(S,1) gwgt(gstep+k: gstep+size(S,1)-1 ) = -100; end
             %gstep = gstep + size(S,1);

             irow = irow + k;
             icol=icol+s;
          end
        end %for ii=1:size(Ppattern,1)

        clear M U;

        fprintf('removing %d rows from constraints due to redundancy\n',size(B,1)-(irow-1));

        B = sparse(I,J,V);
        clear I J V
        ConstraintRhs = newConstraintRhs(1:irow-1);


        if wgtFactor ~= 1
          fprintf('Underweighting %d constraint equations\n',numUnderWeighted);
          wgtVector = diag( sparse( wgtVector(1:irow-1) ) );
          B = wgtVector*B;
        end
        clear wgtVector;

        % Now squeeze out any rows in B that may have been zeroed due to weighting.
        [iB,jB,vB] = find(B);
        [iB,perm] = sort(iB); jB = jB(perm); vB = vB(perm);
        uiB = unique(iB);  % used to pick out that part of the RHS corresponding to nonzero constraints
        pat = ones( size(B,1),1 );
        newrow = 1;
        last = iB(1);
        for ii = 1:length(iB)
          if iB(ii) == last
            iB(ii) = newrow;
          else
            last = iB(ii);
            newrow = newrow+1;
            iB(ii) = newrow;
          end
        end

        [m,n] = size(B);
        B = sparse(iB,jB,vB,newrow,n);
        clear iB jB vB
        ConstraintRhs = ConstraintRhs(uiB);
        clear uiB

        clear newConstraintRhs;
        clear icol irow l k;
      end %if (opt_local_svd)

      % Method 2 : Global SVD on B : simpler than 'local SVD' method but slower
      if (opt_global_svd)
        [U,S,dummy]=svd(full(B));
        clear dummy
        k = nnz(abs(S) > tol);
        clear S V;

        fprintf('removing %d rows from constraints due to redundancy\n',size(B,1)-k);

        B = U(:,1:k)'*B;
        ConstraintRhs = U(:,1:k)'*ConstraintRhs;

        clear k U;
      end

      if (opt_debug)
        disp('B (final) ='); disp(''); disp(full(B));
      end

    end %BuildConstraints()

  end %static methods

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Protected methods                                                         %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  methods (Access=protected)

    function [nDiff] = CompareRootNodeParts(this,AggInfo, P1, P2, RowMap, options)
      %COMPAREROOTNODEPARTS
      %
      %   SYNTAX   nDiff = obj.CompareRootNodeParts(AggInfo, P1, P2, RowMap, options);
      %
      %     AggInfo - Aggregate structure
      %     P1        - first prolongator to compare
      %     P2        - second prolongator to compare
      %     RowMap    - row map common to both prolongators
      %     options   - options array
      %     nDiff     - number of entries that differ in P1 and P2
      NAggregates = max(AggInfo.AggId);
      nulldim = options.NCoarseDofPerNode;

      nDiff=0;
      for i=1:NAggregates
        rnode=AggInfo.Roots(i);
        rows = Node2DOF(rnode,RowMap);
        CDofs   = ((i-1)*nulldim+1:i*nulldim);

        localP1 = P1(rows,CDofs);
        localP2 = P2(rows,CDofs);
        nDiff = nDiff + nnz(abs(localP1 - localP2)>1e-14);
      end
    end %function CompareRootNodeParts()

    function Copy_(this, src, mc)
      %COPY_
      %
      %   SYNTAX   obj.Copy_(src, mc);
      %
      %     src - Object to copy
      %     mc  - MATLAB Metaclass
      [cmd, data, mc] = this.CopyCmd_(src,mc);
      eval(cmd);
    end

    function Matrix2=DecoalesceMatrix(this,Matrix,Nr,Nc)
      % DECOALESCEMATRIX
      %
      % SYNTAX M=obj.DecoalesceMatrix(Matrix,Rb,Cb)
      %
      %   Matrix  - "nodal" matrix
      %   Rb      - number of rows per block
      %   Cb      - number of columns per block
      %   Matrix2 - A matrix with a RbxCb block for each non-zero in Matrix
      %              The entries of matrix2 will be all ones.
        [R,C,V]=find(Matrix);
        Blk=Nr*Nc;
        NNZ=length(R);
        Rblk=reshape(repmat(1:Nr,Nc,1),Blk,1);
        Cblk=reshape(repmat(1:Nc,Nr,1)',Blk,1);

        NewR=reshape(repmat((R-1)*Nr,1,Blk)',NNZ*Blk,1)+reshape(repmat(Rblk',NNZ,1)',NNZ*Blk,1);
        NewC=reshape(repmat((C-1)*Nc,1,Blk)',NNZ*Blk,1)+reshape(repmat(Cblk',NNZ,1)',NNZ*Blk,1);
        V=ones(NNZ*Blk,1);
        Matrix2=sparse(NewR,NewC,V);
      end
  end % protected methods

end %class EminPFactory
