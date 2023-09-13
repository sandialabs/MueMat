classdef SaPFactory < PFactory
    %SAPFACTORY Build a prolongator via the Smoothed Aggregation algorithm
    %
    % Specifically, the following occurs:
    %
    %     1) Produce graph from discretization matrix. This graph reflects dropping
    %        small entries and collapsing block matrices into point matrices.
    %
    %     2) Aggregate vertices of graph.
    %
    %     3) Make tentative prolongator from the aggregates and the null space.
    %
    %     4) Produce final prolongator via prolongator smoothing. That is,
    %           P_final = (I - omega Dinv A) P_tent
    %
    %        where
    %
    %            Dinv is either the inverse of A's point diagonal
    %                 or is the inverse of A's block diagonal.
    %
    %            omega is a damping parameter that is typically
    %                 1.333/lambdaMax(Dinv A)
    %
    properties (Access = private)
        InitPFact_
        diagonalView_ = 'current' % diagonal view label (default == current view)
        damping_
        AForSmoothingName_;
        SANullspace_ = false; % If true, 'fix' coarse nullspace so that Coarse = SmoothP' * Fine;
    end
    methods
        %function [this] = SaPFactory(CoalesceFact, AggFact, diagonalView)
        function [this] = SaPFactory(InitPFact, diagonalView)
            % Copy constructor
            if nargin == 1 && isa(InitPFact, class(this)), this.Copy_(InitPFact,[]); return; end
            %

            % constructor sets options
            if varexist('InitPFact'), this.InitPFact_ = InitPFact;
            else this.InitPFact_ = 'default'; end; %TentativePFactory(); end;

            if varexist('diagonalView'), this.diagonalView_ = diagonalView; end

            if this.GetOutputLevel() > 5
                fprintf('SaPFactory constructor\n');
            end

            this.damping_ = 4/3;
            this.AForSmoothingName_ = 'A';
        end

        function SetDampingFactor(this, value)
            % Set the damping factor used in the prolongator smoother
            this.damping_ = value;
        end

        function SetDiagonalView(this, diagonalView)
            % indicates use of either point or block diagonal prolongator smoothing.
            this.diagonalView_ = diagonalView;
        end

        function SetAForSmoothing(this, name)
            % indicates use of filtered version of A (small entries dropped) within prolognator smoothing
            this.AForSmoothingName_ = name;
        end

       function SetCoarseNullspace(this, TorF)
           % If true, set the coarse nullspace to enforce
           %  CoarseNullSpace = SAP' * FineNullspace;
           % Otherwise, use the coarse nullspace associated with InitP

           this.SANullspace_ = TorF
        end

        function SetNeeds(this, FineLevel, CoarseLevel)
            % check if prolongator is already built
            if this.CheckForReUsableP(FineLevel,CoarseLevel,this) == true
                if this.GetOutputLevel()>5, fprintf('SaPFactory: reUse P\n'); end;
                return;
            end

            TwoLevel = TwoLevels(FineLevel,CoarseLevel);
            TwoLevel.Request('P', this.InitPFact_, CoarseLevel);
            if this.damping_ ~= 0., FineLevel.Request(this.AForSmoothingName_); end
        end

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function flag = Build(this, FineLevel, CoarseLevel)
            % Construct a smoothed aggregation prolongator.
            % See also SaPFactory.
            flag = true;

            % check if prolongator is already built
            if this.CheckForReUsableP(FineLevel,CoarseLevel,this) == true
                if CoarseLevel.IsAvailable('P',this.InitPFact_),
                    CoarseLevel.Release('P', this.InitPFact_);
                end
                return;
            end

            %% Fetch initial guess for prolongator and coarse grid nullspace
            %       TwoLevels = TwoLevel(FineLevel, CoarseLevel);
            %       Ptent = TwoLevels.Get('P', this.InitPFact_, CoarseLevel);

            TwoLevel = TwoLevels(FineLevel,CoarseLevel);
            Ptent = TwoLevel.Get('P', this.InitPFact_, CoarseLevel);
%             if ~isempty(this.InitPFact_)
%                 % Needs InitPFact
%                 this.InitPFact_.Build(FineLevel,CoarseLevel);   % bad: this overwrites CoarseLevel.'P'
%             else
%                 error('SaPFactory::BuildP: no information about Ptent and no InitPFactory set\n');
%             end
%
%             Ptent = CoarseLevel.Get('P', this.InitPFact_);
            CoarseLevel.Release('P', this.InitPFact_);


            %% Smooth prolongator
            if this.damping_ == 0.
                P = Ptent;  % no smoothing
            else
                % Get and Release
                AForSmoothing = this.GetAForSmoothing(FineLevel);

                lambda  = AForSmoothing.GetDinvALambda(this.diagonalView_);
                BlkDiag = AForSmoothing.GetDiagonal([], this.diagonalView_);
                if isempty(BlkDiag.GetApplyInverse())
                    % Note: I prefer to keep this two lines of codes for the
                    % moment but this "if" case *never* appends because computation of
                    % lambda compute also inevitably diagonal factors.
                    BlkDiag.FactorBlkDiag();
                end

                APtent = AForSmoothing * Ptent;
                DinvAPtent = BlkDiag \ APtent;
                P = Ptent - (this.damping_/lambda) * DinvAPtent;
            end

            %% Set P
            this.SetTransferOperator(CoarseLevel,P); % == CoarseLevel.Set('P', P);

            % 'Fix' the coarse nullspace
            if (this.SANullspace_)
              Fine   = FineLevel.Get('NullSpace');
              Coarse = P' * Fine;
            else
              Coarse = CoarseLevel.Get('NullSpace');
            end
            % TODO: generating factory should be used in Level::Set()
            CoarseLevel.Set('NullSpace', Coarse);

        end %Build()

        %%
        function [ToF] = SupportsRestrictionMode(this)
            ToF = true;
        end

        function [A] = GetAForSmoothing(this, FineLevel);
            %GETA
            %
            %   SYNTAX   [A] = GetAForSmoothing(FineLevel);
            %
            %     FineLevel - Level object for fine level
            %     A         - level system matrix for prolongation smoothing
            %
            % default implementation for PFactory derived prolongation operators
            % with support of "restriction" mode
            % The prolongation_mode_ flag controls what is used as system matrix
            % for smoothing the prolongation operator: A in "prolongation" mode
            % and the transposed of A in the "restriction" mode
            if this.ProlongationMode() == true
                % SaPFactory is in prolongation mode
                % use system matrix A
                A = FineLevel.Get(this.AForSmoothingName_);
            else
                % SaPFactory is in restriction mode
                % use the transposed of A (downwinding)
                A = FineLevel.Get(this.AForSmoothingName_)';
            end

            FineLevel.Release(this.AForSmoothingName_);
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
            if this.ProlongationMode() == true
                % SaPFactory is in prolongation mode
                % set prolongation operator
                CoarseLevel.Set('P', TransferOperator, this);
            else
                % SaPFactory is in restriction mode
                % set restriction operator
                CoarseLevel.Set('R', TransferOperator', this);
            end
        end






    end %public methods

    methods (Access = protected)

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

    end % methods

end
