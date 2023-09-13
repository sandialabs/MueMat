classdef Thresholding < PatternFilter
%TRESHOLDING implementation of a simple absolute thresholding pattern filter
% provides |Apply| function
% supports "filtering on fpoints only"

    properties (Access = private)
        FCSplitting_ = [];
        PThreshold_ = 0.0;
        fpointsonly_ = true;  % filtering only at the fpoints
    end

    methods
        function [this] = Thresholding(thresholdparam,FCSplitting)
            % copy constructor
            if nargin == 1 && isa(thresholdparam, class(this)), this.Copy_(thresholdparam,[]); return; end

            if varexist('thresholdparam') this.PThreshold_ = thresholdparam;
            else this.PThreshold_ = 0.0; end;

            if varexist('FCSplitting'), this.FCSplitting_ = FCSplitting;
            else error('Thresholding: no FCSplitting\n'); end;

            % DEPRECATED
            % Need.SavePtent ='true';
%             Need.SaveAggregates ='true';
%             if this.GetOutputLevel() > 5,
%                 fprintf('PatternFilter: Need Tentative prolongator (->AvoidZeroRows)\n');
%                 fprintf('PatternFilter: Need aggregates\n');
%             end
%             this.AddNeeds(Need);

            this.type_ = 'Thresholding';
        end

        function SetNeeds(this, FineLevel, CoarseLevel)
        % Obtain any cross factory specifications
            if ~FineLevel.IsRequested('FCSplitting',this.FCSplitting_),
                this.FCSplitting_.SetNeeds(FineLevel);
            end;
            FineLevel.Request('FCSplitting',this.FCSplitting_);
           CoarseLevel.Request('PtentForFilter');
        end

        function SetThresholdParameter(this, thresholdparam)
            % set threshold parameter
            this.PThreshold_ = thresholdparam;
        end

        function [tACPts] = FilterOnlyFineGridPoints(this, ToF)
            % get/set method: apply filter also at coarse grid nodes
            if varexist('ToF') this.fpointsonly_ = ToF; end;
            tACPts = this.fpointsonly_;
        end

        function [Pattern] = Apply(this, PatternMatrix, FineLevel, CoarseLevel)
        %Apply Apply filter for transfer operator pattern
        %
        % SYNTAX [Pattern] = Apply(PatternMatrix, FineLevel, Coarselevel)
        %
        %   PatternMatrix   - MatlabMatrix matrix with pattern data
        %   FineLevel       - FineLevel object (input)
        %                        provides access:
        %                           - aggregation info AggInfo
        %                           - fine level ID
        %   CoarseLevel      - CoarseLevel object (input)
        %                        provides access:
        %                           - tentative prolongator initialP
        %   Pattern          - binary pattern (output)

            % provide input information
            if ~FineLevel.IsAvailable('FCSplitting', this.FCSplitting_),
                this.FCSplitting_.Build(FineLevel);
            end
            AggInfo  = FineLevel.Get('FCSplitting',this.FCSplitting_); FineLevel.Release('FCSplitting',this.FCSplitting_);
            LevelID  = FineLevel.GetLevelId();
            initialP = CoarseLevel.Get('PtentForFilter'); CoarseLevel.Release('PtentForFilter');

            nFine = size(initialP.GetMatrixData(),1);   % use maps for this? %TODO DOF <-> NODE
            nCoarse = size(initialP.GetMatrixData(),2);
            roots = Node2DOF(AggInfo.cpoints,initialP.GetRowMap()); % determine DOFs of root nodes
            fpoints = ones(nFine,1); fpoints(roots) = 0; fpoints = find(fpoints); % DOFs of fine level nodes
            Pattern = zeros(nFine,nCoarse);
            if this.fpointsonly_
                Pattern(fpoints,:) = (abs(PatternMatrix(fpoints,:) ) > this.PThreshold_);
            else
                Pattern = (abs(PatternMatrix) > this.PThreshold_);
            end
            Pattern = spones(Pattern); % return binary pattern info

            Pattern = this.AvoidZeroRows(Pattern, initialP);

            this.RatioNNZpatternVStentativepattern(Pattern,initialP);
            %if (this.RatioNNZpatternVStentativepattern(Pattern,initialP) < 1.05) && (this.verbose_ > 5)
            %    warning('ratio nnz(pattern)/nnz(tentativepattern) = %f. filter maybe too strong!\n',this.RatioNNZpatternVStentativepattern(Pattern,initialP));
            %end
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
