classdef ConstNNZperRowFilter < PatternFilter
    properties (Access = private)
        numEntries_ = 0;
    end

    methods
        function [this] = ConstNNZperRowFilter(numEntries)
            % copy constructor
            if nargin == 1 && isa(numEntries, class(this)), this.Copy_(numEntries,[]); return; end

            if varexist('numEntries') this.numEntries_ = numEntries;
            else this.numEntries_ = 5; end;

% DEPRECATED
%             Need.SavePtent ='true';
%             Need.SaveAggregates ='true';
%             if this.GetOutputLevel() > 5,
%                 fprintf('PatternFilter: Need Tentative prolongator (->AvoidZeroRows)\n');
%                 fprintf('PatternFilter: Need aggregates\n');
%             end
%             this.AddNeeds(Need);

            this.type_ = 'ConstNumberOfEntriesPerRow';
        end

        function SetNumEntriesParameter(this, numEntries)
            this.numEntries_ = numEntries;
        end

        function [tACPts] = FilteringAtCoarsePoints(this, ToF)
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
            AggInfo  = FineLevel.Get('Aggregates');
            LevelID  = FineLevel.GetLevelId();
            initialP = CoarseLevel.Get('Ptent');

            nFine = size(initialP.GetMatrixData(),1);   % use maps for this? %TODO DOF <-> NODE
            nCoarse = size(initialP.GetMatrixData(),2);
            Pattern = zeros(nFine,nCoarse);
            roots = Node2DOF(AggInfo.Roots,initialP.GetRowMap()); % determine DOFs of root nodes
            fpoints = ones(nFine,1); fpoints(roots) = 0; fpoints = find(fpoints); % DOFs of fine level nodes

            % strategy that prescribes # entries per row
            for i=1:length(fpoints)
                currow = PatternMatrix(fpoints(i),:);
                rowidx = find(currow); % number of entries in row
                rowval = currow(rowidx);
                rowval = sort(abs(rowval),'descend');

                if(length(rowval) > this.numEntries_ )
                    %rowval(3)

                    rowidxelim = find(abs(currow) < rowval(this.numEntries_));
                    PatternMatrix(fpoints(i),rowidxelim) = 0;

                    if length(find(PatternMatrix(fpoints(i),:))) ~= this.numEntries_;
                        fprintf('pattern in row %i has %i entries?\n',fpoints(i),length(find(PatternMatrix(fpoints(i),:))));
                    end
                end
            end

            Pattern = spones(PatternMatrix); % return binary pattern info

            Pattern = this.AvoidZeroRows(Pattern, initialP);
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
