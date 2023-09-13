%% MinDescentPFactory
% concrete implementation of PFactory that provides prolongator using some
% energy minimization ideas combined with a simple Descent minimization
% algorithm
%%
% the following steps occur
% # build fine level nullspace
% # create initial prolongator
% # create prescribed sparsity pattern for prolongator (PatternFactory)
% # prepare projection operators for nullspace preservation ans sparsity
% constraint
% # run minimization routine of |AP|_F using a MinDescent variant
%
% * MinIter_ minimization iterations are performed (constant, adaptive step
% length strategy)
% * StepLength is prescribed (constant step length strategy)
%
% MinDescentPFactory has support for prolongation and restriction mode (you
% can use it with GenericRFactory)



classdef MinDescentPFactory < PFactory
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties                                                                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (Access=private)
%         PatternFact_        = [];

        PatternFact_ = []       % aggregation
        InitPFact_   = []       % PFactory for initial guess

        QR_                  = true;

        % number of minimization iterations
        MinIter_             = [];

        % constant StepLengthStrategy_: step length
        % adaptive: initial step length
        % global/local: empty
        StepLength_          = [];

        % StepLengthStrategy_ can be constant, adaptive, global, local
        StepLengthStrategy_  = [];

        % can be true or false
        % only optimization for non dirichlet Dofs?
        % true: we do not prescribe a wrong nullspace approximation for the
        % dirichlet bcs (default)
        % false: just "optimize" prolongator/restrictor for all Dofs even
        % though our optimal prolongation operators are disturbing the
        % solver at dirichlet bc dofs.
        IgnoreDirichletDofs_ = true;

        % initial projection in restriction mode?
        InitialProjectionForRestrictor_ = false; %% FIXME!!!

        reUseP_          = false; %TODO option not available for
                                  %the moment
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Public methods                                                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
        %function [this] = MinDescentPFactory(PatternFact, AggregationFact, StepLengthStrategy)
        function [this] = MinDescentPFactory(PatternFact, InitPFact, StepLengthStrategy)
            %MinDescentPFactory Constructor
            %
            %   SYNTAX obj = MinDescentPRFactory(PatternFact, AggregationFact, StepLengthStrategy)
            %
            %    PatternFact        - pattern factory object (for prescribed
            %    sparsity pattern
            %    AggregationFact    - for aggregation
            %    StepLengthStrategy - type = string, can be either 'constant',
            %                         'adaptive', 'global' (default) or 'local',
            %                         defines step length strategy that is used
            %                         for minimization descent algorithm
            %
            if varexist('PatternFact')  this.PatternFact_  = PatternFact;
            else this.PatternFact_ = AP_PatternFactory(); end;
            if varexist('InitPFact'), this.InitPFact_ = InitPFact;
            else this.InitPFact_ = TentativePFactory(); end;

            % internal logic: check if pattern factory knows its Factory
            % for the initial prolongator
            % if not: set the same InitPFact that is used for this PFactory
            if isempty(this.PatternFact_.GetInitialPFactory())
                this.PatternFact_.SetInitialPFactory(this.InitPFact_);
            end

            if varexist('StepLengthStrategy')
               if strcmp(StepLengthStrategy, 'constant') || ...
                  strcmp(StepLengthStrategy, 'adaptive') || ...
                  strcmp(StepLengthStrategy, 'global')   || ...
                  strcmp(StepLengthStrategy, 'local')
                   this.StepLengthStrategy_ = StepLengthStrategy;
               else
                  error('unknown StepLengthStrategy %s, must be either constant, adaptive, global or local',StepLengthStrategy);
               end
            else
                this.StepLengthStrategy_ = 'global'; % default
            end
            this.MinIter_            = 5;    % default values
            this.StepLength_         = 0.1;

            %fprintf('constructor MinDescentPRFactory\n');
        end

        function [ToF] = SupportsRestrictionMode(this)
            ToF = true;
        end

        function [ToF] = InitialProjectionForRestrictor(this, ToF)
          if varexist('ToF'),
              ToFold = this.InitialProjectionForRestrictor_;
              this.InitialProjectionForRestrictor_ = ToF;
              ToF = ToFold;
          else ToF = this.InitialProjectionForRestrictor_; end
        end

        function [ToF] = ReUseP(this, ToF)
          if varexist('ToF'),
              ToFold = this.reUseP_;
              this.reUseP_ = ToF;
              ToF = ToFold;
          else ToF = this.reUseP_; end
        end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

        function SetNeeds(this, FineLevel, CoarseLevel)
        % Obtain any cross factory specifications

            FineLevel.Request('NullSpace'); % only needed for finest level???
            CoarseLevel.Request('NullSpace'); % needed all over the algorithm. freed in BuildP
            CoarseLevel.Request('P',this.InitPFact_); % only needed for GetInitialP

            CoarseLevel.Request('BtBFact'); % needed by SatisfyConstraints (quite often)
                                            % generated by CalculateBtBFact

            if this.prolongation_mode_ == true
                if ~CoarseLevel.IsAvailable('P', this.InitPFact_)
                    % prepare for generating initial prolongator
                    this.InitPFact_.SetNeeds(FineLevel,CoarseLevel);
                end

                this.PatternFact_.SetNeeds(FineLevel,CoarseLevel);

                CoarseLevel.Request('Ppattern',this.PatternFact_);

                if strcmp(this.StepLengthStrategy_,'global') || ...
                   strcmp(this.StepLengthStrategy_,'adaptive')
                    CoarseLevel.Request('MinDescent_StepLengthReduction');  % will be freed in restrictor
                    CoarseLevel.Request('MinDescent_MinIter');
                end
            end


        end

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

        function [P] = GetInitialP(this, FineLevel, CoarseLevel)
        % returns CoarseLevel.'Pinitial' and CoarseLevel.'NullSpace'
        % if there's no 'Pinitial' variable, create one...

          if(~CoarseLevel.IsAvailable('P', this.InitPFact_))
                if this.prolongation_mode_ == false
                   warning('no initial P information (generated by InitPFact_) in restriction mode?? very suspicious!\n');
                end
                this.InitPFact_.BuildP(FineLevel,CoarseLevel);   % bad: this overwrites CoarseLevel.'P'
          end
          P = CoarseLevel.Get('P',this.InitPFact_);
          CoarseLevel.Release('P',this.InitPFact_);

        end

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
        function SetMaxCoarseSize(this, MaxCoarseSize)
            this.MaxCoarseSize_ = MaxCoarseSize;
        end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
        function SetPatternFactory(this,PatternFact)
            % sets pattern factory
            %
            %   SYNTAX SetPatternFactory(PatternFact)
            %
            %    PatternFact  - PatternFactory for prescribed pattern
            this.PatternFact_ = PatternFact;
        end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
        function SetMinimizationIters(this, MinIters)
            % sets number of minimization iterations within descent method
            %
            %   SYNTAX SetMinimizationIters(MinIters)
            %
            %    MinIters  - number of minimization iterations performed
            %                within descent minimization routine
            this.MinIter_ = MinIters;
        end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
        function SetStepLength(this,steplength)
            % set fixed step length for minimization iteration (only for
            % constant steplength strategy)
            %
            %   SYNTAX SetStepLength(steplength)
            %
            %    steplength  - sets step length for use within descent
            %                  minimization routine (only for constant step
            %                  length strategy)

            this.StepLength_ = steplength;
        end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
        function SetStepLengthStrategy(this, strat)
            % set step length strategy
            %
            %   SYNTAX SetStepLengthStrategy(strat)
            %
            %    strat       - sets step length strategy

            this.StepLengthStrategy_ = strat;
        end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
        function TentativeWithQR(this, value)
            % Specify whether orthogonalization performed when forming tentative prolongator
            this.QR_ = value;
        end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
        function [ToF] = IgnoreDirichletDofs(this, ToF)
            % ignore DiricheletDofs within descent minimization routine?
            % (get/set function)
            %
            %   SYNTAX ToF = IgnoreDirichletDofs(ToF)
            %
            %    ToF  - ignore DirichletDofs within descent minimization
            %           routine (true/false)

            if varexist('ToF') this.IgnoreDirichletDofs_ = ToF; end;
            ToF = this.IgnoreDirichletDofs_;
        end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
        function flag = BuildP(this,FineLevel,CoarseLevel)
            % Build transfer operators using some descent mehtod for
            % optimizing the energy of transfer operators
            %
            % SYNTAX flag = Build(FineLevel, CoarseLevel)
            %
            %    FineLevel     - FineLevel (Level)
            %    CoarseLevel   - coarse level (Level)
            flag = true;

            %% get initial guess for prolongator/restrictor and coarse grid nullspace
            % runs e.g. aggregation process and builds the tentative
            % prolongator
            % Note: In GetInitialP also Ptent is saved on the coarse level
            % (-> needed for transfer operator pattern e.g. AP_Pattern)
            % that is, GetInitialP has to be called before generation of
            % transfer operator sparsity pattern...
            [Pinitial] = this.GetInitialP(FineLevel,CoarseLevel);

            %% build pattern for transfer operators
            % if there's no pattern stored as 'Ppattern' in CoarseLevel, we
            % generate a new transfer operator pattern using the
            % PatternFactory.
            % Note: The pattern factory itself uses information from
            % FineLevel and CoarseLevel for generating the transfer
            % operator pattern (independent of prolongation or restriction
            % mode!!!). We suppose that a pattern is only generated for the
            % prolongator and properly reused for the restrictor.
            if ~CoarseLevel.IsAvailable('Ppattern',this.PatternFact_)
                if this.prolongation_mode_ == false, warning('no transfer operator pattern available in restriction mode?? very suspicious!'); end;
                Ppattern = this.PatternFact_.Build(FineLevel,CoarseLevel);
            else
                % reuse existing transfer operator pattern (given by the
                % user or generated in the prolongation mode)
                if this.prolongation_mode_ == true, warning('pattern available but not in restriciton mode?'); end;
                Ppattern = CoarseLevel.Get('Ppattern',this.PatternFact_);
                CoarseLevel.Release('Ppattern',this.PatternFact_);
            end



            %% calculate BtB
            % needed for projection operator of constraint P * Bone = Bzero
            % is calculated only once in the prolongation mode and reused
            % in the restriction mode (needed by SatisfyConstraints)
            if this.prolongation_mode_ == true
                CalculateBtBFact(this,FineLevel,CoarseLevel,Ppattern);
            end

            %% prepare MinDescent transfer operator improvement method
            % copy initial guess for transfer operator
            % P stores the improved prolongation/restriction operator
            P = Pinitial.Copy();
            PP = P.GetMatrixData();

            %% initial projection (constraints)
            % project down the initial prolongator (restrictor?) to fulfill
            % the constraints for sparsity pattern and nullspace
            % preservation
            if this.prolongation_mode_ == true
                %% prescribe sparsity pattern
                PP = PP .* Ppattern;

                %% constraint P*Bone = Bzero
                % force P*Bone = Bzero (corresponding to the update updateP =
                % (Bzero-P*Bone)*inv(Bone'Bone)*Bone' we obtain P = P + updateP
                % = P + (Bzero-P*Bone)*inv(Bone'Bone)*Bone' = P *
                % (I-Bone*inv(Bone'Bone)*Bone') + Bzero*inv(Bone'Bone)*Bone' =
                % \Pi_B(P) in [Wiesner, "On multigrid transfer operators for
                % nonsymmetric positive definite systems", formula (4.38)]
                if(CoarseLevel.IsAvailable('NullSpace'))
                    Bzero = FineLevel.Get('NullSpace');
                    Bone = CoarseLevel.Get('NullSpace');  % we need coarse nullspace more than once!
                    NQ = PP * Bone - Bzero;
                    PP = this.SatisfyConstraints(FineLevel, CoarseLevel,PP, Ppattern, NQ);
                else
                    error('no coarse level nullspace?');
                end

                %% TODO: check this option
                % store PP that fulfills all constraints for reuse within
                % restriction mode (see below)
                CoarseLevel.Set('PP',PP);

            else
                %% TODO: this is a new option for restriction mode?????? check me!!!
                if this.InitialProjectionForRestrictor_ == true
                    PP = CoarseLevel.Get('PP');
                end
            end

            %% now do "optimization"
            switch lower(this.StepLengthStrategy_)
                case 'constant'
                    [PP] = this.MinDescentConstant(FineLevel,CoarseLevel,PP,Ppattern);
                case 'adaptive'
                    [PP] = this.MinDescentAdaptive(FineLevel,CoarseLevel,PP,Ppattern);
                case 'global'
                    [PP] = this.MinDescentGlobal(FineLevel,CoarseLevel,PP,Ppattern);
                %case 'local'
                    %[PP,RR] = this.MinDescentLocal(AA,PP,Ppattern,RR,Rpattern,Bzero,Bone,BtBFact,nPDE);
                otherwise
                    error('unknown StepLengthStrategy: %s, must be either constant, adaptive or global',this.StepLengthStratgy_);
            end

            % now we can free the FineLevel.NullSpace
            FineLevel.Release('NullSpace');
            CoarseLevel.Release('BtBFact');
            CoarseLevel.Release('NullSpace');

            %% catch spurious zero rows in PP
            % (robustness!)
            % fill these rows and columns with the values of Pinitial
            % This is the same as setting the omega values to zero for these
            % rows and columns (restriction mode)
            zerorows = find(max(abs(PP'))==0);
            if ~isempty (zerorows)
                if this.GetOutputLevel() > 5, fprintf('repair %i zero rows in P \n', length(zerorows)); end;
                PPqr = Pinitial.GetMatrixData(); %Ptent.GetMatrixData();
                PP(zerorows,:) = PPqr(zerorows,:);
                zerorows = find(max(abs(PP'))==0);
                if ~isempty(zerorows)
                    error('repair of zero rows in P failed?');
                end
            end

            %% write back PP to MueMat operator P
            P.SetMatrixData(PP);

            %% set prolongator/restrictor
            this.SetTransferOperator(CoarseLevel,P);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Static methods                                                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Static = true)
        function [dir] = SatisfyConstraints(FineLevel, CoarseLevel, dir, Ppattern, NQ)
            if ~CoarseLevel.IsAvailable('BtBFact')
                error ('no BtBFact data stored in CoarseLevel! Call CalculateBtBFact first!');
            end
            BtBFact = CoarseLevel.Get('BtBFact'); % we need BtBFact more than once!
            Bone = CoarseLevel.Get('NullSpace');  % we need coarse nullspace more than once!

            Amat  = FineLevel.Get('A'); % note: always "plain" A (only needed for nPDE)
            nPDE = Amat.GetRowMap().ConstBlkSize();
            if (Amat.GetRowMap().HasVariableBlkSize()) error('we only support constant block sizes'); end;

            if ~varexist('NQ')
                % default
                NQ = dir * Bone;
            end

            % check dimension of nullspace
            if(size(Bone,2) > 1)
                [nrows,ncols] = size(dir);
                % loop over all NDOFs
                % note: we have only NNodes BtBFact operators!
                % TODO: we should loop over all NNodes and handle these
                % things blockwise!!
                for i=1:nrows
                    [aaa,bbb,ccc] = find(Ppattern(i,:)); % get row pattern
                    % FIXME: this is the row to block information because
                    % we're currently looping over all rows
                    % Furthermore we're implicitly assuming that the we
                    % have constant block size.
                    Rtt = BtBFact(:,:,ceil(i/nPDE)); % this is Bone factorized for Bone'*Bone
                    % inv(Bone) * inv(Bone') * (Bzero - P * Bone)'
                    t = Rtt\(Rtt'\NQ(i,:)');
                    % note: we obtain formula (4.38) since
                    % t' * Bone' = (Bzero - P*Bone) inv(Bone) inv(Bone')
                    % * Bone' = (Bzero - P*Bone) inv(Bone' Bone) * Bone'
                    dir(i,bbb) = dir(i,bbb) - t'*Bone(bbb,:)';  % update
                end
            else
                % for 1dim nullspace BtBFact currently contains
                % Bone*inv(Bone'Bone)*Bone'
                % note: BtBFact is a tensor third order, but for 1dim
                % nsp saved as a simple matrix -> multiply diagonal
                % matrix with NQ as diagonal as replacement for "tensor
                % product"
                Nnodes = length(NQ);
                dir = dir - spdiags(NQ,[0],Nnodes,Nnodes) * BtBFact;
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Private methods                                                           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = private)
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
        function CalculateBtBFact(this,FineLevel,CoarseLevel,Ppattern)
            %CALCULATEBTBFACT
            %
            %   SYNTAX CalculateBtBFact(FineLevel,CoarseLevel,Ppattern)
            %
            %    FineLevel       - FineLevel
            %    CoarseLevel     - CoarseLevel
            %    Ppattern        - transfer operator pattern (prolongator
            %                      shape, MATLAB matrix)
            %
            % calculates BtBFact for transfer operator projection.
            % has to be calculated only once in prolongation mode.
            % result is stored as 'BtBFact' in CoarseLevel.BtBFact and is
            % loaded/used by SatisfyConstraints

            Bzero = FineLevel.Get('NullSpace');
            Bone  = CoarseLevel.Get('NullSpace');
            nulldim = size(Bzero,2);
            Nnodes = this.GetA(FineLevel).GetRowMap().NNodes();
            nPDE = this.GetA(FineLevel).GetRowMap().ConstBlkSize();

            if(nulldim>1)
                % factorize Nnodes (nulldim x nulldim) matrices needed to satisfy constraints
                BtBFact = zeros(nulldim,nulldim,Nnodes);  % stores compressed tensor
                for i=1:Nnodes    % loop over all Nnodes equality constraint sets
                    [aaa,bbb,ccc] = find(Ppattern((i-1)*nPDE+1:i*nPDE,:)); % pick out correct rows from Ppattern
                    if(~isempty(bbb))
                        RestrictedNull = Bone(bbb,:);    % = B^{k+1}_{nnz(i)} % pick out correct nnz rows from Bone
                        if length(find(eig(full(RestrictedNull'*RestrictedNull))<=0)>0)
                            fprintf('Problem with RestrictedNull\n');
                        end
                        BtBFact(:,:,i) = chol(RestrictedNull'*RestrictedNull);
                        % obtain (nulldim x nulldim) matrix from cholesky
                        % decomposition that belongs to the i-th set of eq
                        % constraints
                        % BtBFact (..,i) is a triangular matrix and
                        % represents a factor in the factorization of
                        % (RestrictedNull'*RestrictedNull), that we need for
                        % computing inv(Bone' Bone) = inv(Bone) inv(Bone')
                        % later
                    end
                end
            elseif(nulldim==1)
                [aaa,bbb,ccc] = find(Ppattern);
                z = sparse(aaa,bbb,Bone(bbb),Nnodes,size(Bone,1)); % = B^{k+1}_{nnz}
                BtBFact = full(diag(z*z'));
                BtBFact(find(BtBFact == 0)) = 1;
                BtBFact = spdiags(1./BtBFact,[0],Nnodes,Nnodes)*z;
                % note: z corresponds now to Bone*inv(B^t B)*Bone'!
            else
                BtBFact = [];
                error('upps. this is bad.');
            end

            CoarseLevel.Set('BtBFact',BtBFact);
        end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
        function [PP] = MinDescentConstant(this,FineLevel,CoarseLevel,PP,Ppattern)
            %MINDESCENTCONSTANT
            %
            %   SYNTAX PP = MinDescentConstant(FineLevel,CoarseLevel,PP,Ppattern)
            %
            %    FineLevel       - FineLevel
            %    CoarseLevel     - CoarseLevel
            %    PP              - MATLAB matrix with initial transfer
            %                      operator data for energy minimization
            %                      (prolongator shape)
            %    Ppattern        - transfer operator pattern (prolongator
            %                      shape, MATLAB matrix)
            %
            % use constant step length (this.StepLength_)

            Bzero = FineLevel.Get('NullSpace');
            Bone  = CoarseLevel.Get('NullSpace');
            nPDE = this.GetA(FineLevel).GetRowMap().ConstBlkSize();
            BtBFact = CoarseLevel.Get('BtBFact');

            AA  = this.GetA(FineLevel).GetMatrixData();
            ddd = spdiags(1./diag(AA),[0],size(AA,1),size(AA,2));    % inverse of point diagonal of A

            for k=1:this.MinIter_
                %% search direction for prolongator
                dir = -ddd * AA * PP;        % 1) search direction
                dir = dir .* Ppattern;       % 2) sparsity pattern
                dir2= this.SatisfyConstraints(FineLevel, CoarseLevel, dir,Ppattern); %) satisfy (O2) and (O3)

                %% correct search direction (prolongator)
                % ignore dirichlet dofs, i.e. no disturbance for dirichlet
                % dofs?
                if this.IgnoreDirichletDofs_ == true
                    correctnsp = find(abs(AA*Bzero)<0.2);
                    dir(correctnsp,:) = dir2(correctnsp,:);
                    %dir(correctnsp) = dir2(correctnsp); % TODO error!
                else
                    dir = dir2;
                end

                %% update
                PP = PP + this.StepLength_ * dir;

                if this.GetOutputLevel() > 5
                    if this.prolongation_mode_ == true
                        fprintf('MinDescent(%s) iter (prolongation): %d/%d steplength %10.5e, norm_Fro (A*P) %10.7e \n',this.StepLengthStrategy_,k,this.MinIter_,this.StepLength_,norm(AA*PP,'fro'));
                    else
                        fprintf('MinDescent(%s) iter (restriction ): %d/%d steplength %10.5e, norm_Fro (A*P) %10.7e \n',this.StepLengthStrategy_,k,this.MinIter_,this.StepLength_,norm(AA*PP,'fro'));
                    end
                end
            end

        end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
        function [PP] = MinDescentAdaptive(this,FineLevel,CoarseLevel,PP,Ppattern)
            %MINDESCENTADAPTIVE
            %
            %   SYNTAX PP = MinDescentAdaptive(FineLevel,CoarseLevel,PP,Ppattern)
            %
            %    FineLevel       - FineLevel
            %    CoarseLevel     - CoarseLevel
            %    PP              - MATLAB matrix with initial transfer
            %                      operator data for energy minimization
            %                      (prolongator shape)
            %    Ppattern        - transfer operator pattern (prolongator
            %                      shape, MATLAB matrix)
            %
            % use constant step length (this.StepLength_) with simple step
            % length control routine

            Bzero = FineLevel.Get('NullSpace');
            Bone  = CoarseLevel.Get('NullSpace');
            nPDE = this.GetA(FineLevel).GetRowMap().ConstBlkSize();
            BtBFact = CoarseLevel.Get('BtBFact');

            AA  = this.GetA(FineLevel).GetMatrixData();
            steplength = this.StepLength_;
            ddd = spdiags(1./diag(AA),[0],size(AA,1),size(AA,2));    % inverse of point diagonal of A

            lastAPnorm = inf; % only needed for adaptive StepLengthStrategy

            for k=1:this.MinIter_
                %% search direction for prolongator
                dir = -ddd * AA * PP;        % 1) search direction
                dir = dir .* Ppattern;       % 2) sparsity pattern
                dir2= this.SatisfyConstraints(FineLevel, CoarseLevel, dir,Ppattern); %) satisfy (O2) and (O3)

                %% correct search direction (prolongator)
                % ignore dirichlet dofs, i.e. no disturbance for dirichlet
                % dofs?
                if this.IgnoreDirichletDofs_ == true
                    correctnsp = find(abs(AA*Bzero)<0.2);
                    dir(correctnsp,:) = dir2(correctnsp,:);
                    %dir(correctnsp) = dir2(correctnsp); % TODO error
                else
                    dir = dir2;
                end

                %% simple step length control routine
                % reduce step length for reduction of norm
                if this.prolongation_mode_ == true
                    numStepLengthShortening = 0;
                    curAPnorm = norm (AA * (PP+steplength*dir),'fro');
                    while (curAPnorm > lastAPnorm)
                        steplength = 0.5 * steplength;
                        numStepLengthShortening = numStepLengthShortening + 1;
                        curAPnorm = norm (AA * (PP+steplength*dir),'fro');
                        if(steplength < 1e-10)
                            break;
                        end;
                    end
                    lastAPnorm = curAPnorm;
                    CoarseLevel.Set('MinDescent_StepLengthReduction',numStepLengthShortening);
                else
                    % iterations with full steplength until
                    % MinDescent_MinIter is exceeded
                    if CoarseLevel.IsAvailable('MinDescent_MinIter')
                        if k >= CoarseLevel.Get('MinDescent_MinIter')
                            % MinDescent_MinIter exceeded -> step length
                            % control routine -> reduce step length
                            numStepLengthShortening = CoarseLevel.Get('MinDescent_StepLengthReduction');
                            steplength = (0.5^numStepLengthShortening) * steplength;
                        end
                    end
                end

                %% update
                PP = PP + steplength * dir;

                %% step length break condition & output
                if this.prolongation_mode_ == true

                    CoarseLevel.Set('MinDescent_MinIter',this.MinIter_);
                    if(steplength < 1e-10)
                        CoarseLevel.Set('MinDescent_MinIter',k);
                        break;
                    end
                    if this.GetOutputLevel() > 5
                        fprintf('MinDescent(%s) iter (prolongation): %d/%d steplength %10.5e, norm_Fro (A*P) %10.7e \n',this.StepLengthStrategy_,k,this.MinIter_,steplength,norm(AA*PP,'fro'));
                    end
                else
                    if(steplength < 1e-10)
                       break;
                    end
                    if this.GetOutputLevel() > 5
                        fprintf('MinDescent(%s) iter (restriction ): %d/%d steplength %10.5e, norm_Fro (A*P) %10.7e \n',this.StepLengthStrategy_,k,this.MinIter_,steplength,norm(AA*PP,'fro'));
                    end

                end
            end

            if this.prolongation_mode_ == false,
                 CoarseLevel.Release('MinDescent_MinIter');
                 CoarseLevel.Release('MinDescent_StepLengthReduction');
            end
        end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
        function [PP] = MinDescentGlobal(this,FineLevel,CoarseLevel,PP,Ppattern)
            %MINDESCENTGLOBAL
            %
            %   SYNTAX PP = MinDescentAdaptive(FineLevel,CoarseLevel,PP,Ppattern)
            %
            %    FineLevel       - FineLevel
            %    CoarseLevel     - CoarseLevel
            %    PP              - MATLAB matrix with initial transfer
            %                      operator data for energy minimization
            %                      (prolongator shape)
            %    Ppattern        - transfer operator pattern (prolongator
            %                      shape, MATLAB matrix)
            %
            % uses global step length strategy with simple steplengh
            % control routine

            Bzero = FineLevel.Get('NullSpace');
            Bone  = CoarseLevel.Get('NullSpace');
            nPDE = this.GetA(FineLevel).GetRowMap().ConstBlkSize();
            BtBFact = CoarseLevel.Get('BtBFact');

            AA  = this.GetA(FineLevel).GetMatrixData();
            ddd = spdiags(1./diag(AA),[0],size(AA,1),size(AA,2));    % inverse of point diagonal of A

            lastAPnorm = inf; % only needed for adaptive StepLengthStrategy

            for k=1:this.MinIter_
                %% search direction for prolongator
                dir = -ddd * AA * PP;        % 1) search direction
                dir = dir .* Ppattern;       % 2) sparsity pattern
                dir2= this.SatisfyConstraints(FineLevel, CoarseLevel, dir,Ppattern); %) satisfy (O2) and (O3)

                %% correct search direction (prolongator)
                % ignore dirichlet dofs, i.e. no disturbance for dirichlet
                % dofs?
                % works only if nulldim == 1!!!
                if this.IgnoreDirichletDofs_ == true
                    correctnsp = find(abs(AA*Bzero)<0.2);
                    dir(correctnsp,:) = dir2(correctnsp,:);
                    %dir(correctnsp) = dir2(correctnsp); % TODO error
                else
                    dir = dir2;
                end

                %% recalculate global steplength within every iteration
                % same idea as PG-AMG (same formulas, too)
                steplength = full(sum(diag(PP'*AA'*AA*ddd*AA*PP)) / sum(diag(PP'*AA'*ddd*AA'*AA*ddd*AA*PP)));

                %% simple step length control routine
                % reduce step length for reduction of norm
                if this.prolongation_mode_ == true
                    numStepLengthShortening = 0;
                    curAPnorm = norm (AA * (PP+steplength*dir),'fro');
                    while (curAPnorm > lastAPnorm)
                        steplength = 0.5 * steplength;
                        numStepLengthShortening = numStepLengthShortening + 1;
                        curAPnorm = norm (AA * (PP+steplength*dir),'fro');
                        if(steplength < 1e-10)
                            break;
                        end;
                    end
                    lastAPnorm = curAPnorm;
                    CoarseLevel.Set('MinDescent_StepLengthReduction',numStepLengthShortening);
                else
                    % iterations with full steplength until
                    % MinDescent_MinIter is exceeded
                    if CoarseLevel.IsAvailable('MinDescent_MinIter')
                        if k >= CoarseLevel.Get('MinDescent_MinIter')
                            % MinDescent_MinIter exceeded -> step length
                            % control routine -> reduce step length
                            numStepLengthShortening = CoarseLevel.Get('MinDescent_StepLengthReduction');
                            steplength = (0.5^numStepLengthShortening) * steplength;
                        end
                    end
                end

                %% update
                PP = PP + steplength * dir;

                %% step length break condition & output
                if this.prolongation_mode_ == true
                    CoarseLevel.Set('MinDescent_MinIter',this.MinIter_);
                    if(steplength < 1e-10)
                        CoarseLevel.Set('MinDescent_MinIter',k);
                        break;
                    end
                    if this.GetOutputLevel() > 5
                        fprintf('MinDescent(%s) iter (prolongation): %d/%d steplength %10.5e, norm_Fro (A*P) %10.7e \n',this.StepLengthStrategy_,k,this.MinIter_,steplength,norm(AA*PP,'fro'));
                    end
                else
                    if k > CoarseLevel.Get('MinDescent_MinIter') % free MinDescent_MinIter
                        break;
                    end;
                    if this.GetOutputLevel() > 5
                        fprintf('MinDescent(%s) iter (restriction ): %d/%d steplength %10.5e, norm_Fro (A*P) %10.7e \n',this.StepLengthStrategy_,k,this.MinIter_,steplength,norm(AA*PP,'fro'));
                    end

                end
            end

            % free some variables
            if this.prolongation_mode_ == false
                 CoarseLevel.Release('MinDescent_MinIter');
                 CoarseLevel.Release('MinDescent_StepLengthReduction');
            end
        end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%         function [PP,RR] = MinDescentLocal(this,AA,PP,Ppattern,RR,Rpattern,Bzero,Bone,BtBFact,nPDE)
%             ddd = spdiags(1./diag(AA),[0],size(AA,1),size(AA,2));    % inverse of point diagonal of A
%
%             lastAPnorm = inf; % only needed for adaptive StepLengthStrategy
%
%             for k=1:this.MinIter_
%                 %% search direction for prolongator
%                 dir = -ddd * AA * PP;        % 1) search direction
%                 if(size(RR,1) > size(RR,2))
%                     RR = RR';
%                 end
%                 dir_r = -ddd * AA' * RR';
%
%                 %% calculate optimal step length for unconstrained directions!!
%                 nc = size(PP,2);
%                 omegajp = zeros(nc,1);
%                 omegajr = zeros(nc,1);
%                 for j=1:nc
%                     Pj = PP(:,j);
%                     Rj = RR(j,:);
%
%                     omegajp(j) = (Pj'*AA'*AA*ddd*AA*Pj) / (Pj'*AA'*ddd*AA'*AA*ddd*AA*Pj);
%                     omegajr(j) = (Rj *AA *AA'*ddd*AA'*Rj') / (Rj*AA*ddd*AA*AA'*ddd*AA'*Rj');
%
%                     if(isnan(omegajp(j))==1)
%                         omegajp(j) = full(sum(diag(PP'*AA'*AA*ddd*AA*PP)) / sum(diag(PP'*AA'*ddd*AA'*AA*ddd*AA*PP)));
%                     end
%
%                     if(isnan(omegajr(j))==1)
%                         omegajr(j) = full(sum(diag(RR*AA*AA'*ddd*AA'*RR'))/ sum(diag(RR*AA*ddd*AA*AA'*ddd*AA'*RR')));
%                     end
%
%                 end
%
%                 omegajp=spdiags(omegajp,0,nc,nc);
%                 omegajr=spdiags(omegajr,0,nc,nc);
%
%                 %% adapt unconstrained directions with steplength (~local steplength)
%                 dir = dir * omegajp;
%                 dir_r = dir_r * omegajr;
%
%                 %% satisfy constraints
%                 dir = dir .* Ppattern;
%                 dir2= this.SatisfyConstraints(FineLevel, CoarseLevel, dir,Ppattern); % satisfy (O2) and (O3)
%
%                 %% correct search direction (prolongator)
%                 % ignore dirichlet dofs, i.e. no disturbance for dirichlet
%                 % dofs?
%                 if this.IgnoreDirichletDofs_ == true
%                     correctnsp = find(abs(AA*Bzero)<0.2);
%                     %dir(correctnsp,:) = dir2(correctnsp,:);
%                     dir(correctnsp) = dir2(correctnsp);
%                 else
%                     dir = dir2;
%                 end
%
%                 %% search direction for restrictor
%                 dir_r = dir_r .* Rpattern';
%                 dir2_r= this.SatisfyConstraints(FineLevel, CoarseLevel, dir_r,Rpattern');
%
%                 %% correct search direction (restrictor)
%                 % ignore dirichlet dofs, i.e. no disturbance for dirichlet
%                 % dofs?
%                 if this.IgnoreDirichletDofs_ == true
%                     correctnsp = find(abs(AA'*Bzero)<0.2);
%                     %dir_r(correctnsp,:) = dir2_r(correctnsp,:);
%                     dir_r(correctnsp) = dir2_r(correctnsp);
%                 else
%                     dir_r = dir2_r;
%                 end
%
%                 %% break condition
%                 curAPnorm = norm (AA * (PP+dir),'fro');
%                 if (curAPnorm > lastAPnorm) && this.GetOutputLevel() > 5
%                     fprintf('MinDescent(local): steps: %i ||AP||_F = %10.5e\n',k,norm(AA*PP,'fro'));
%                     break;
%                 end
%
%                 % update
%                 PP = PP + dir;
%                 RR = RR + dir_r';
%
%                 lastAPnorm = curAPnorm;
%
%                 if this.GetOutputLevel() > 5
%                     fprintf('MinDescent(%s) iter: %d/%d, norm_Fro (A*P) %10.7e \n',this.StepLengthStrategy_,k,this.MinIter_,norm(AA*PP,'fro'));
%                 end
%             end
%         end
    end

    methods (Access = protected)

      % TODO: this functions should be moved to PFactory
      function [A] = GetA(this, FineLevel);
      %GETA
      %
      %   SYNTAX   [A] = GetA(FineLevel);
      %
      %     FineLevel - Level object for fine level
      %     A         - level system matrix for prolongation smoothing
      %
      % default implementation for PFactory derived prolongation operators
      % with support of "restriction" mode
      % The prolongation_mode_ flag controls what is used as system matrix
      % for smoothing the prolongation operator: A in "prolongation" mode
      % and the transposed of A in the "restriction" mode
          if this.prolongation_mode_ == true
              % PgPFactory is in prolongation mode
              % use system matrix A
              A = FineLevel.Get('A');
          else
              % PgPFactory is in restriction mode
              % use the transposed of A (downwinding)
              A = FineLevel.Get('A')';
          end
      end
      function SetTransferOperator(this,CoarseLevel,TransferOperator)
      %SETTRANSFEROPERATOR
      %
      %   SYNTAX   SetTransferOperator(CoarseLevel, TransferOperator;
      %
      %     CoarseLevel      - Level object for coarse level
      %     TransferOperator - prolongator (or transposed of restrictor in
      %                        "restriction" mode)
      %
      % stores the result of transfer operator smoothing in level data
      % structure as prolongator (if prolongation_mode_==true) or as
      % restrictor (if prolongation_mode_==false).
      % This is the default implementation for prolongation operator with
      % support of "restriction mode".
          if this.prolongation_mode_ == true
              % PgPFactory is in prolongation mode
              % set prolongation operator
              CoarseLevel.Set('P', TransferOperator, this);
          else
              % PgPFactory is in restriction mode
              % set restriction operator
              CoarseLevel.Set('R', TransferOperator', this);
          end
      end
        function Copy_(this,src,mc)
            % COPY_
            %  SYNTAX obj.Copy_(src,mc);
            %  src: Object to copy
            %  mc: MATLAB Metaclass
            [cmd, data, mc] = this.CopyCmd_(src,mc);
            eval(cmd);
        end
    end
end