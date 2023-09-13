%% OptimalPFactory
% factory class for ideal prolongation (and restriction operator)
% which explicitely calculates
%             / -Aff^{-1} Afc \
%   P^{opt} = |               |
%             \      I        /
%
% This class is meant only for testing purposes and some studies
% The computational costs (computational time and memory footprint)
% extremely high!

classdef OptimalPFactory < PFactory
    % class implements optimal prolongation and restriction operator
    properties (Access = private)
        AggFact_ = [];
        ConstraintFact_ = [];
        reUseP_ = false;
    end

    methods
        function [this] = OptimalPFactory(AggregationFact,ConstraintFactory)
            if nargin == 1 && isa(AggregationFact, class(this)), this.Copy_(AggregationFact,[]); return; end
            if varexist('AggregationFact'), this.AggFact_ = AggregationFact;
            else this.AggFact_ = AggregationFactory(); end

            if varexist('ConstraintFactory'),
                this.ConstraintFact_ = ConstraintFactory;
                if isempty(this.ConstraintFact_.GetPatternFactory())
                    warning('OptimalPFactory with constraints: ConstraintFactory has no pattern! use standard AP pattern!');
                    this.ConstraintFact_.SetPatternFactory(AP_PatternFactory());
                end
            end

            if this.GetOutputLevel() > 5
                fprintf('TentativePFactory constructor\n');
            end
        end

        function SetNeeds(this, FineLevel, CoarseLevel)
            % Obtain any cross factory specifications
            %if ~this.ReUseP()
            this.AggFact_.SetNeeds(FineLevel);
            FineLevel.Request('Aggregates');

            FineLevel.Request('NullSpace');

            if ~isempty(this.ConstraintFact_)
                if ~CoarseLevel.IsAvailable('B') || ~CoarseLevel.IsAvailable('ConstraintRhs'),
                    this.ConstraintFact_.SetNeeds(FineLevel, CoarseLevel);
                end
                CoarseLevel.Request('ConstraintRhs');
                CoarseLevel.Request('Ppattern',this.ConstraintFact_.GetPatternFactory());
            end
        end

        function [ToF] = SupportsRestrictionMode(this) % should be static
            ToF = true;
        end

        function [ToF] = ReUseP(this, ToF)
            if varexist('ToF'),
                ToFold = this.reUseP_;
                this.reUseP_ = ToF;
                ToF = ToFold;
            else ToF = this.reUseP_;
            end
        end

        function flag = BuildP(this,FineLevel,CoarseLevel)
            flag = true;

            % check out system matrix A
            if this.prolongation_mode_ == true
                Amat = FineLevel.Get('A');
            else
                Amat = FineLevel.Get('A')';
            end
            % Check out fine level nullspace
            if ~FineLevel.IsAvailable('NullSpace');
                NS = BuildNullSpace(Amat);
                FineLevel.Set('NullSpace', NS);
            else
                NS = FineLevel.Get('NullSpace');
            end
            FineLevel.Release('NullSpace');

            % Build aggregates
            this.AggFact_.Build(FineLevel);
            AggInfo = FineLevel.Get('Aggregates');
            FineLevel.Release('Aggregates');

            roots = AggInfo.Roots;
            fpoints = ones(size(Amat.GetMatrixData(),1),1);
            fpoints(roots) = 0; fpoints = find(fpoints);

            Naggs   = max(AggInfo.AggId);
            nulldim = size(NS,2);
            RowMap = Amat.GetRowMap();
            rootDofs = Node2DOF(roots,RowMap);
            fpointDofs = Node2DOF(fpoints,RowMap);

            % todo point -> DOF information!
            AA = Amat.GetMatrixData();
            Aff = AA(fpointDofs,fpointDofs);
            Afc = AA(fpointDofs,rootDofs);
            Affinv = inv(Aff);

            % coarse level nullspace -> just injection
            cnull = NS(rootDofs,:);
            CoarseLevel.Set('NullSpace', cnull);

            % build prolongator
            PP = sparse(RowMap.NDOFs(),Naggs*nulldim);
            PP(fpointDofs,:) = -Affinv * Afc;
            PP(rootDofs,:) = eye(size(PP,2));

            P = Operator(PP, RowMap, Map(Naggs,nulldim), @MatlabApply);
            if ~isempty(this.ConstraintFact_)
                this.ConstraintFact_.Build(FineLevel, CoarseLevel); % fnull for rhs???
                PatFact = this.ConstraintFact_.GetPatternFactory();
                Rhs = CoarseLevel.Get('ConstraintRhs');
                Pattern = CoarseLevel.Get('Ppattern', this.ConstraintFact_.GetPatternFactory());
                CoarseLevel.Release('Ppattern', this.ConstraintFact_.GetPatternFactory());
%                 Rhs = sparse(NS - P*cnull);
                P = EminPFactory.SatisfyConstraints(P,Rhs,CoarseLevel,Pattern)
                CoarseLevel.Release('ConstraintRhs');
            end

            this.SetTransferOperator(CoarseLevel, P);

        end
    end

    methods (Access = private)
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

    methods (Access = protected)
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