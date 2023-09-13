% TODO: remove argument AggFact - Not used inside of this (only there for compatibility reasons)
% TODO: merge with TentativePFactory. TentativePFactory has less option (NCoarseDofPerNode option not take into account, etc.)

classdef TentativePFactoryEx < PFactory
    properties (Access = private)
        % a little strange but Aggregation
        AggFact_        % factories now appear in the ApproxP class?
        CNullFact_      %
        options_        % all of the options here are pumped into the
        % null space factory. I wanted to make minimal changes
        % to the MueMat repository so I left this here.
        % However, options should definitely go away from
        % ApproxP and probably be replaced by something else
        % in CNullFact_.
        diagonalView_            = 'current'
    end

    methods
        function [this] = TentativePFactoryEx(CoalesceFact, AggFact, CNullFact, options, diagonalView)
            % Copy constructor
            if nargin == 1 && isa(CoalesceFact, class(this)), this.Copy_(CoalesceFact,[]); return; end

            if varexist('AggFact'), this.AggFact_ = AggFact;
            else                    this.AggFact_ = AggregationFactory(); end;

            % for compatibility reasons
            if varexist('CoalesceFact') && ~isempty(CoalesceFact), this.AggFact_.SetCoalesceFactory(CoalesceFact); end

            if varexist('CNullFact'), this.CNullFact_ = CNullFact;
            else                      this.CNullFact_ = CoarseNSFactory(); end;

            if varexist('options'), this.options_ = options;
            else                    this.options_ = []; end;

            if varexist('diagonalView'), this.SetDiagonalView(diagonalView); end

            if this.GetOutputLevel() > 5
                fprintf('TentativePFactoryEx constructor\n');
            end
        end

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

        function SetNeeds(this, FineLevel, CoarseLevel)
            if this.CheckForReUsableP(FineLevel,CoarseLevel,this) == true, return; end

            FineLevel.Request('Aggregates',this.AggFact_);
            FineLevel.Request('NullSpace');

            if ~isempty(this.AggFact_), this.AggFact_.SetNeeds(FineLevel); end
            if ~isempty(this.CNullFact_), this.CNullFact_.SetNeeds(FineLevel, CoarseLevel); end
        end

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

        function [ToF] = SupportsRestrictionMode(this) % should be static
            ToF = false;
        end

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

        function flag = Build(this,FineLevel,CoarseLevel)
            flag = true;

            if this.CheckForReUsableP(FineLevel,CoarseLevel,this) == true, return; end

            % Check out system matrix, that is used for aggregation
            A = FineLevel.Get('A');

            % Build aggregates
            this.AggFact_.Build(FineLevel);
            Aggregates = FineLevel.Get('Aggregates',this.AggFact_);

            % Check out fine level nullspace
            if ~FineLevel.IsAvailable('NullSpace');
                FineNullSpace = BuildNullSpace(A);
                FineLevel.Set('NullSpace', FineNullSpace);
            else
                FineNullSpace = FineLevel.Get('NullSpace');
            end

            %% smooth fine level nullspace
            % Jacob's nullspace massage to improve energy minimization
            % near the boundaries
            if isfield(this.options_,'SmoothNullspace'),
                fprintf('TentativePFactoryEx: Smoothing fine level nullspace \n');
                Adata = A.GetMatrixData();
                D  = spdiags(diag(Adata),0,size(Adata,1),size(Adata,2));
                DA = D\Adata;
                lambda_max = A.GetDinvALambda(this.diagonalView_);

                for iii=1:this.options_.SmoothNullspace
                    FineNullSpace = FineNullSpace - 4/(3*lambda_max)*(DA*FineNullSpace);
                end;
                FineNullSpace = normalize(FineNullSpace);
            end

            %% build coarse level nullspace and "tentative" prolongator
            % CoarseNSFactory allows for permutations of nullspace vectors
            % and grid coordinates as optional parameters
            perm = []; % no permutation %% TODO permutation from EminPFactory??
            if FineLevel.IsAvailable('xcoords') && FineLevel.IsAvailable('ycoords') && FineLevel.IsAvailable('zcoords')
                coords=[FineLevel.Get('xcoords'),FineLevel.Get('ycoords'),FineLevel.Get('zcoords')];
                [P,CoarseNullSpace] = this.CNullFact_.Build(Aggregates, A, FineNullSpace, this.options_,this.GetOutputLevel(),perm,coords);
            else
                [P,CoarseNullSpace] = this.CNullFact_.Build(Aggregates, A, FineNullSpace, this.options_,this.GetOutputLevel());
            end

            %% Output
            CoarseLevel.Set('P', P, this);
            CoarseLevel.Set('NullSpace', CoarseNullSpace);

            %% Release
            FineLevel.Release('A');
            FineLevel.Release('Aggregates',this.AggFact_);
            FineLevel.Release('NullSpace');

        end %Build

    end


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
