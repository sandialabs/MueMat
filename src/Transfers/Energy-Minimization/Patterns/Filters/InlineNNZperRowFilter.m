classdef InlineNNZperRowFilter < PatternFilter
    properties (Access = private)
        inlinefct_  = inline('1');
        fpointsonly_ = true;
        FCSplitting_ = [];
    end

    methods
        function [this] = InlineNNZperRowFilter(inlinefct,FCSplitting)
        %InlineNNZperRowFilter constructor
        %
        %  SYNTAX: [obj] = InlineNNZperRowFilter(inlinefct)
        %
        %    inlinefct  - input: number of average nnz per row
            if nargin == 1 && isa(inlinefct, class(this)), this.Copy_(inlinefct,[]); end;

            if varexist('inlinefct') this.inlinefct_ = inline(inlinefct,'LevelID');
            else this.inlinefct_ = inline('1'); end;

            if varexist('FCSplitting'), this.FCSplitting_ = FCSplitting;
            else error('Thresholding: no FCSplitting\n'); end;

            this.type_ = 'InlineNNZperRowFilter';
        end

        function SetNeeds(this, FineLevel, CoarseLevel)
        % Obtain any cross factory specifications
            if ~FineLevel.IsRequested('FCSplitting',this.FCSplitting_),
                this.FCSplitting_.SetNeeds(FineLevel);
            end;
            FineLevel.Request('FCSplitting',this.FCSplitting_);

           CoarseLevel.Request('PtentForFilter');
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

        function SetInlineNNZFunction(this, inlinefct)
        %SetInlineNNZFunction set internal relative thresholding parameter
        %
        %  SYNTAX: SetInlineNNZFunction(inlinefct)
        %
        %    inlinefct  - input: average number of nnz per row
            this.inlinefct_ = inline(inlinefct,'LevelID');
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
            if ~FineLevel.IsAvailable('FCSplitting',this.FCSplitting_),
                this.FCSplitting_.Build(FineLevel);
            end
            AggInfo  = FineLevel.Get('FCSplitting',this.FCSplitting_); FineLevel.Release('FCSplitting',this.FCSplitting_);
            LevelID  = FineLevel.GetLevelId();
            initialP = CoarseLevel.Get('PtentForFilter');
            CoarseLevel.Release('PtentForFilter');

            nFine = size(initialP.GetMatrixData(),1);   % number of dofs
            nCoarse = size(initialP.GetMatrixData(),2); % number of dofs
            roots = Node2DOF(AggInfo.cpoints, initialP.GetRowMap());  % DOFs of root nodes
            fpoints=ones(nFine,1); fpoints(roots) = 0; fpoints = find(fpoints); % DOFs of fine grid nodes

            Pattern = sparse(nFine,nCoarse); %zeros(nFine,nCoarse);

            if this.fpointsonly_
                % we adapt only the pattern of the fine level nodes!
                PatternMatrix = PatternMatrix(fpoints,:);   % restrict me to fpoints!
                allowednnz = ceil(this.inlinefct_(LevelID) * length(fpoints));
            else
                allowednnz = ceil(this.inlinefct_(LevelID) * nFine);
            end

            % find for every row at least 1 entry!
            %             maxabsval = norm(PatternMatrix,inf);
            %             for r=1:size(PatternMatrix,1)   % loop over all lines
            %                 rowvals = PatternMatrix(r,:);
            %                 if(length(find(rowvals) > 1))
            %                     rowData = sort(abs(rowvals(find(PatternMatrix(r,:)))),'descend');
            %                     PatternMatrix(r,find(rowvals>=rowData(1))) = maxabsval + 1;
            %                 else
            %                     % empty row
            %                     Ptent = initialP.GetMatrixData();
            %                     PatternMatrix(r,:) = (maxabsval+1) * spones(Ptent(r,:));
            %                 end
            %             end

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
