classdef PatternFilter < VerboseObject
%PATTERNFILTER base class for a transfer operator pattern filter
%
%  provides virtual functions for pattern filters (Apply) and common service
%  functions (e.g. AvoidZeroRows)

    properties (Access = private)
    end
    properties (Access = protected)
        type_ = [];
    end

    methods
        function [this] = PatternFilter(obj)
            % pattern filter

            % copy constructor
            if nargin==1 && isa(obj, class(this)), this.Copy_(obj,[]); return; end

            % declare needs: Ptent needed for fallback method
% DEPRECATED
%             Need.SavePtent ='true';
%             if this.GetOutputLevel() > 5,
%                 fprintf('PatternFilter: Need Tentative prolongator (->AvoidZeroRows)\n');
%             end
%             this.AddNeeds(Need);
        end

        function [type] = GetFilterType(this)
            % returns filter type (e.g. Thresholding)
            type = this.type_;
        end

        function [Pattern] = AvoidZeroRows(this, Pattern, initialP)
        % AvoidZeroRows search for zero rows in given prolongator pattern
        % and fills zero rows with corresponding information from initialP
        % operator
        %
        % SYNTAX [Pattern] = AvoidZeroRows(Pattern, initialP)
        %
        %   Pattern         - Matlab matrix with pattern data
        %   initialP        - tentative prolongator of type Operator
        % output: binary matlab matrix for pattern

            % make sure, that there are no zero rows in Pattern
            % fill up zero rows with information from Ptent!
            % for most problems this should be ok (Ptent = 1 entry per row)
            % maybe prolematic for nulldim > 1? -> elasticity?
            Ptent = initialP.GetMatrixData();

            if this.GetOutputLevel() > 5
                fprintf('PatternFilter(%s) zerorows(Ptent)=%i zerorows(Ppattern)=%i \n',this.type_,length(find(max(abs(Ptent'))==0)),length(find(max(abs(Pattern'))==0)));
                fprintf('PatternFilter(%s) nnz(Ptent)=%i nnz(Ppattern)=%i ',this.type_,nnz(Ptent),nnz(Pattern));
            end

            % avoid zero rows
            zerorows = find(full(sum(Pattern'))==0.0);
            Pattern(zerorows,:) = abs(Ptent(zerorows,:));

            if this.GetOutputLevel() > 5
                fprintf('PatternFilter(%s)nnz(Ptent+Ppattern)=%i\n',this.type_,nnz(Pattern));
                fprintf('PatternFilter(%s)zerorows(Ptent+Ppattern)=%i\n',this.type_,length(find(max(abs(Pattern'))==0)));
            end

            Pattern = spones(Pattern); % all entries in Pattern 0 or 1
        end

        function [ratio] = RatioNNZpatternVStentativepattern(this, Pattern, initialP)
        % RatioNNZpatternVStentativepattern calculates ratio
        % nnz(Pattern)/nnz(initialP) for current level
        % if ratio = 1 the filter is too strong -> tentative prolongator
        % pattern
        %
        %   SYNTAX [ratio] = RatioNNZpatternVStentativepattern(Pattern, initialP)
        %
        %       Pattern      - pattern after filtering
        %       initialP     - tentative prolongator of type Operator
            ratio = nnz(Pattern) / nnz(initialP.GetMatrixData());

            if ratio < 1.05 && this.GetOutputLevel() > 2
               warning('PatternFilter(%s) filtered pattern (almost) the same as tentative pattern (nnz ratio = %f)! filter too strong!\n',this.type_,ratio);
            end
        end
    end

    methods(Abstract = true)
        [Pattern] = Apply(PatternMatrix, FineLevel, CoarseLevel); %Apply Apply filter for transfer operator pattern (virtual)
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
