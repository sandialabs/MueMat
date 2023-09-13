classdef AvgNNZperRowFilter < PatternFilter
%AvgNNZperRowFilter
% provides |Apply| function

    properties (Access = private)
        avgnnz_  = 0;
        fpointsonly_ = true;
    end

    methods
        function [this] = AvgNNZperRowFilter(avgnnz)
        %AvgNNZperRowFilter constructor
        %
        %  SYNTAX: [obj] = AvgNNZperRowFilter(avgnnz)
        %
        %    avgnnz  - input: number of average nnz per row
            if nargin == 1 && isa(avgnnz, class(this)), this.Copy_(avgnnz,[]); end;

            if varexist('avgnnz') this.avgnnz_ = avgnnz;
            else this.avgnnz_ = 1; end;

% DEPRECATED
%             Need.SavePtent ='true';
%             Need.SaveAggregates ='true';
%             if this.GetOutputLevel() > 5,
%                 fprintf('PatternFilter: Need Tentative prolongator (->AvoidZeroRows)\n');
%                 fprintf('PatternFilter: Need aggregates\n');
%             end
%             this.AddNeeds(Need);

            this.type_ = 'AvgNNZperRowFilter';
        end

        function [tACPts] = FilterOnlyFineGridPoints(this, ToF)
        %FilterOnlyFineGridPoints get/set function for fpointsonly_
        %property. if true, the filter is only applied for fine grid
        %nodes and the coarse grid dofs are not touched
        %
        %  SYNTAX: [tACPts] = FilterOnlyFineGridPoints(ToF)
        %
        %   ToF    - input: true or false (can be empty)
        %   tACPts - output: true or false
            if varexist('ToF') this.fpointsonly_ = ToF; end;
            tACPts = this.fpointsonly_;
        end

        function SetAvgNNZParameter(this, avgnnz)
        %SetAvgNNZParameter set internal relative thresholding parameter
        %
        %  SYNTAX: SetAvgNNZParameter(avgnnz)
        %
        %    avgnnz  - input: average number of nnz per row
            this.avgnnz_ = avgnnz;
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
            AggInfo  = FineLevel.Get('Aggregates');
            LevelID  = FineLevel.GetLevelId();
            initialP = CoarseLevel.Get('Ptent');

            nFine = size(initialP.GetMatrixData(),1);   % number of dofs
            nCoarse = size(initialP.GetMatrixData(),2); % number of dofs
            roots = Node2DOF(AggInfo.Roots, initialP.GetRowMap());  % DOFs of root nodes
            fpoints=ones(nFine,1); fpoints(roots) = 0; fpoints = find(fpoints); % DOFs of fine grid nodes

            Pattern = zeros(nFine,nCoarse);

            if this.fpointsonly_
                % we adapt only the pattern of the fine level nodes!
                PatternMatrix = PatternMatrix(fpoints,:);   % restrict me to fpoints!
                allowednnz = ceil(this.avgnnz_ * length(fpoints));
            else
                allowednnz = ceil(this.avgnnz_ * nFine);
            end

            % sort all nnz entries of PatternMatrix
            patData = sort(abs(PatternMatrix(find(PatternMatrix))),'descend');
            if length(patData) < nFine, warning('pattern data suspicious! pattern nnz = %i < nFine = %i\n',length(patData),nFine); end;
            allowednnz = min(allowednnz,length(patData));

            absoluteThreshold = patData(allowednnz);

            PatternMatrix = (abs(PatternMatrix) >= absoluteThreshold);

            PatternMatrix = spones(PatternMatrix);

            if this.fpointsonly_
               Pattern(fpoints,:) = PatternMatrix;
            else
               Pattern = PatternMatrix;
            end

            Pattern = this.AvoidZeroRows(Pattern, initialP);

            this.RatioNNZpatternVStentativepattern(Pattern,initialP);
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
