classdef EminPFactory < PFactory
    properties (Access = private)
        EnergyMatrix_                = 'A';
        PatternFact_                 = [];
        InitPFact_                   = [];
        ConstraintFact_              = []
        EminSolver_                  = [];

        diagonalView_                = 'current';

        options_                     = []; % options for Emin and other factories
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Public methods                                                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods
        % TODO: should EnergyMatrix argument replaced by a factory?

        function [this] = EminPFactory(PatternFact, ConstraintFact, EminSolver, InitPFact, options, EnergyMatrix, diagonalView)
            %EMINPFACTORY Constructor
            %
            %   SYNTAX   obj = EminPFactory(PatternFact, ConstraintFact, EminSolver, ...
            %                               InitPFact, EnergyMatrix, diagonalView)
            %
            %     PatternFact     -  sparsity pattern factory
            %     ConstraintFact  -  interpolation constraint builder
            %     EminSolver      -  EminSolver handle for solving the 2x2 system
            %     InitPFact       -  initial prolongator guess
            %     EnergyMatrix    -  function which defines matrix used for energy (string)
            %     DiagonalView    -  view label

            % copy constructor
            if nargin == 1 && isa(PatternFact, class(this)), this.Copy_(PatternFact,[]); return; end

            if varexist('ConstraintFact'), this.ConstraintFact_ = ConstraintFact; end
            if varexist('EminSolver'),     this.EminSolver_     = EminSolver;     end
            if varexist('options'),        this.options_ = options;               end
            if varexist('EnergyMatrix'),   this.EnergyMatrix_   = EnergyMatrix;   end

            if varexist('PatternFact')  this.PatternFact_  = PatternFact;
            else error('EminPFactory: no pattern factory\n'); end;
            if varexist('InitPFact'), this.InitPFact_ = InitPFact;
            else this.InitPFact_ = 'default'; end;

            % internal logic: check if pattern factory knows its Factory
            % for the initial prolongator
            % if not: set the same InitPFact that is used for this PFactory
            if ~isempty(this.PatternFact_) && ...
                isempty(this.PatternFact_.GetInitialPFactory())
                this.PatternFact_.SetInitialPFactory(this.InitPFact_);
            end

            % internal logic: set pattern factory for constraints and
            % initial P factory to be the same as here
            if ~isempty(this.ConstraintFact_) && ...
                isempty(this.ConstraintFact_.GetPatternFactory())
                this.ConstraintFact_.SetPatternFactory(this.PatternFact_);
                this.ConstraintFact_.SetOptions(this.options_);
            end

            if varexist('diagonalView'),   this.SetDiagonalView(diagonalView);    end
        end

        function SetEminSolverFunction(this, EminSolver)
            this.EminSolver_         = EminSolver;
        end
        function SetConstraintFactory(this, ConstraintFact)
            this.ConstraintFact_ = ConstraintFact;
        end
        function [ToF] = SupportsRestrictionMode(this)
            ToF = true;
        end

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

        function SetNeeds(this, FineLevel, CoarseLevel)
            if this.CheckForReUsableP(FineLevel,CoarseLevel,this) == true, return; end

            FineLevel.Request(this.EnergyMatrix_);

            TwoLevel = TwoLevels(FineLevel,CoarseLevel);

            % 1) request for initial P
            TwoLevel.Request('P', this.InitPFact_, CoarseLevel);

            % 2) request for pattern constraint
            if ~isempty(this.PatternFact_), TwoLevel.Request('Ppattern', this.PatternFact_, CoarseLevel); end;


            if ~isempty(this.ConstraintFact_),
                % 3) request for constraints helper variables
                TwoLevel.Request('B',this.ConstraintFact_, CoarseLevel);
                CoarseLevel.Request('ConstraintRhs', this.ConstraintFact_);

                % 4) data generated locally by SatisfyConstraints
                CoarseLevel.Request('BBt');
                CoarseLevel.Request('PInds');
            end;


            % 5) request for EminSolver needs
            if ~isempty(this.EminSolver_)
                this.EminSolver_.SetNeeds(FineLevel,CoarseLevel);
            end

        end

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

        function SetDiagonalView(this, diagonalView)
            % indicates use of either point or block diagonal prolongator smoothing.
            %this.diagonalView_ = diagonalView;
            warning('SetDiagonalView() has no effect right now');
        end

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

        function SetEnergyMatrix(this, EnergyMatrix)
            this.EnergyMatrix_ = EnergyMatrix;
        end

        function [EnergyMatrix] = GetEnergyMatrix(this, FineLevel, CoarseLevel)

            if(this.prolongation_mode_ == true)
                EnergyMatrix = FineLevel.Get(this.EnergyMatrix_);
            else
                EnergyMatrix = FineLevel.Get(this.EnergyMatrix_)';
            end
        end

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
        function flag = Build(this, FineLevel, CoarseLevel)
            % Build transfer operators using some descent mehtod for
            % optimizing the energy of transfer operators
            %
            % SYNTAX flag = Build(FineLevel, CoarseLevel)
            %
            %    FineLevel     - FineLevel (Level)
            %    CoarseLevel   - coarse level (Level)
            flag = true;

            if this.CheckForReUsableP(FineLevel,CoarseLevel,this) == true, return; end

            TwoLevel = TwoLevels(FineLevel,CoarseLevel);

            % Initial Guess
            Pinitial = TwoLevel.Get('P', this.InitPFact_,CoarseLevel);
            this.PrintStats('Pinitial', Pinitial, FineLevel.Get('NullSpace'), CoarseLevel.Get('NullSpace'));

            % Define A for energy minimization
            EnergyMatrix = this.GetEnergyMatrix(FineLevel,CoarseLevel);

            % Constraint 1: Sparsity pattern
            if ~isempty(this.PatternFact_), Ppattern = TwoLevel.Get('Ppattern', this.PatternFact_, CoarseLevel); end;

            % Constraint 2: NullSpace (Constraint matrix)
            B = []; ConstraintRhs = [];
            if ~isempty(this.ConstraintFact_),
                B = TwoLevel.Get('B', this.ConstraintFact_, CoarseLevel);
                ConstraintRhs = TwoLevel.Get('ConstraintRhs', this.ConstraintFact_, CoarseLevel);
            end;

            % Enforce constraints for EminSolver initial guess
            Pinitial = this.SatisfyConstraints(Pinitial, ConstraintRhs, CoarseLevel, Ppattern, B);
            this.PrintStats('Pinitial after enforcing constraints', Pinitial, FineLevel.Get('NullSpace'), CoarseLevel.Get('NullSpace'));

            % Energy minimization
            P = this.EminSolver_.Iterate(EnergyMatrix, Ppattern, B, Pinitial, FineLevel, CoarseLevel, this.prolongation_mode_);
            this.PrintStats('P after emin', P, FineLevel.Get('NullSpace'), CoarseLevel.Get('NullSpace'));

            % Fix coarse nullspace
            %CoarseLevel.Set('NullSpace', P'*FineLevel.Get('NullSpace'));

            %% 9) catch spurious zero rows in P
            % (robustness!)
            % fill these rows and columns with the values of Pinitial
            % This is the same as setting the omega values to zero for these
            % rows and columns

            % only for test!
            %               zerorows = find(max(abs(P.GetMatrixData()'))==0);
            %               if ~isempty (zerorows)
            %                   fprintf('repair %i zero rows in P \n', length(zerorows));
            %                   PP = P.GetMatrixData();
            %
            %                   PPqr = Pinitial2.GetMatrixData();
            %
            %                   PP(zerorows,:) = PPqr(zerorows,:);
            %                   P = Operator(PP,P.GetRowMap(),P.GetColMap(),@MatlabApply);
            %                   zerorows = find(max(abs(PP'))==0);
            %                   if ~isempty(zerorows)
            %                       warning('repair of zero rows in P failed?'); % see Axels email
            %                       %2/5/11
            %                   end
            %               end


            %% only diagnostics
            % rst: These should be diagonally scaled.
            fprintf('condest(Pfinal''Pfinal) = %10.2e\n', condest(P.GetMatrixData()'*P.GetMatrixData()));
            % Need 'A': fprintf('condest(P''AP) = %10.2e\n', condest(P.GetMatrixData()'*FineLevel.Get('A').GetMatrixData()*P.GetMatrixData()));
            % rst: this should probably invoke something in the near null constraint class which
            % rst: indicates how well the constraints are satisfied.
            %             fprintf('norm(P.cnull-fnull)=%f, nnzP = %d\n', norm(P*cnull-NS,inf), nnz(P.GetMatrixData()));
            %             dif = P*cnull-NS;
            %             fprintf('column by column norm(FNull - P(energy)*CNull) = ');
            %             for jj=1:size(NS,2), fprintf('%10.2e ',norm(dif(:,jj))); end
            %             fprintf('\n');
            %%

            % Output
            this.SetTransferOperator(CoarseLevel,P);

            %             % rst: there is some scaling but somehow it seems that we should do some row-wise
            %             % rst: scaling. (Note: it might be nice to use something other than diagonal scaling
            %             % rst: within the energy minimization algorithms .. which reappears here). That is,
            %             % rst:     trace(diag(PAP)/diag(P'*P))
            %             % rst:
            %             % rst: Note: Pcons is the P matrix before energy minimization
            %             D=diag(diag(FineLevel.Get('A').GetMatrixData()));
            %             PtAPt=Pinitial.GetMatrixData()'*(D\(FineLevel.Get('A').GetMatrixData()*Pinitial.GetMatrixData));
            %             PAP=P.GetMatrixData()'*(D\(FineLevel.Get('A').GetMatrixData()*P.GetMatrixData()));
            %
            % %             if ~isfield(this.options_,'NCoarseDofPerNode') Ncpn = size(NS,2);
            % %             else Ncpn = this.options_.NCoarseDofPerNode;end
            %
            %             % Total energy diagnostics
            %             energy=trace(PAP);
            %             tent_energy = trace(PtAPt);
            %             fprintf('Initial energy (sum_i pt_i^T A pt_i)     = %6.4e\n',full(tent_energy));
            %             fprintf('Final energy (sum_i p_i^T A p_i)        = %6.4e\n',full(energy));
            %

            %% Release
            CoarseLevel.Release('P',this.InitPFact_)
            if ~isempty(this.PatternFact_), CoarseLevel.Release('Ppattern',this.PatternFact_); end;
            if ~isempty(this.ConstraintFact_),
                CoarseLevel.Release('B', this.ConstraintFact_);
                CoarseLevel.Release('ConstraintRhs', this.ConstraintFact_);
                CoarseLevel.Release('BBt');
                CoarseLevel.Release('PInds');
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
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Public methods (static)                                                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods (Static = true)

        function PrintStats(name, P, fnull, cnull)
            dif = P*cnull-fnull;
            fprintf('====== %s ======\n', name);
            fprintf('norm(P.cnull-fnull) = %f, nnzP = %d\n', norm(dif,inf), nnz(P.GetMatrixData()));
            fprintf('column by column norm(FNull - P*CNull) = ');
            for jj=1:size(fnull,2), fprintf('%10.2e ', norm(dif(:,jj))); end
            fprintf('\n');
            fprintf('===============\n', name);
            %keyboard
            %normalize(full([P*cnull fnull]))
        end

        function [Pvec] = Flatten(Pmat, PatternInds)
            Pvec  = reshape(Pmat',[],1);
            Pvec  = Pvec(PatternInds);
        end %Flatten()

        function [Unflattened] = BlowUp(Vec,nfine,ncoarse,PInds);
            nn = length(Vec);
            BigVec = sparse(PInds,ones(nn,1),Vec,nfine*ncoarse,1,nn);
            Unflattened = reshape(BigVec,ncoarse,nfine)';
        end %BlowUp()

        %function P = SatisfyConstraints(P_, rhs, CoarseLevel, PatternFact)
        function P = SatisfyConstraints(P_, rhs, CoarseLevel, Pattern, B)
            % satisfies constraints
            % requires:
            % 2 x Ppattern (from PatternFact)
            % 3 x B
            % sets:
            % 1 x BBt (if not already set)
            % 1 x Pinds (if not already set)
            if strcmp('Operator',class(P_)), P = P_.GetMatrixData(); else P=P_; end

            % 1) pattern constraint
            if ~isempty(Pattern), P = P .* Pattern; end;

            % 2) projection step for nullspace preservation
            if ~isempty(B),
                if ~CoarseLevel.IsAvailable('BBt')
                    BBt = B * B';
                    fprintf('(before scaling) condest(BBt) = %g\n',condest(BBt,5));
                    CoarseLevel.Set('BBt', BBt);
                end
                if ~CoarseLevel.IsAvailable('PInds')
                    PInds =  find(reshape(Pattern',[],1));
                    CoarseLevel.Set('PInds', PInds);
                    % took out scaling. Can't recall why we need it.
                end

                [nrows,ncols] = size(P);
                Pflat = EminPFactory.Flatten(P,CoarseLevel.Get('PInds'));
                temp = B*Pflat;
                if ~isempty(rhs), temp = temp - rhs; end
                P = EminPFactory.BlowUp(Pflat-B'*(CoarseLevel.Get('BBt')\temp),...
                    nrows,ncols,CoarseLevel.Get('PInds'));
            end;

            if strcmp('Operator',class(P_)), P = Operator(P, P_.GetRowMap(), P_.GetColMap(), P_.GetApply()); end
        end
    end

    methods (Access = protected)

        function Copy_(this, src, mc)
            %COPY_
            %
            %   SYNTAX   obj.Copy_(src, mc);
            %
            %     src - Object to copy
            %     mc  - MATLAB Metaclass
            [cmd, dummy, mc] = this.CopyCmd_(src,mc);
            eval(cmd);
        end

    end % methods
end
