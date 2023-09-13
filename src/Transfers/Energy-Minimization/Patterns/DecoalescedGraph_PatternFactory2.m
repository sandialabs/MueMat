%% DecoalescedGraph_PatternFactory
% Pattern = initialP * AfilteredDecoalesced

classdef DecoalescedGraph_PatternFactory2 < PatternFactory
    properties (Access = private)

    end

    methods
        function [this] = DecoalescedGraph_PatternFactory2(InitPFact)
            % copy constructor
            if nargin == 1 && isa(InitPFact, class(this)), this.Copy_(InitPFact,[]); return; end

            if varexist('InitPFact'), this.InitPFact_ = InitPFact;
            else error('DecoalescedGraph_Pattern: no InitPFact given!\n'); end

            this.type_ = 'DecoalescedGraph_PatternFactory2'; % What is that??
        end

        function SetNeeds(this, FineLevel, CoarseLevel)
           %TODO
        end

        function [Pattern] = BuildPattern(this, FineLevel, CoarseLevel)

            initialP = CoarseLevel.Get('P', this.InitPFact_);
            Amatrix = FineLevel.Get('AfilteredDecoalesced');

            % MATLAB only:
            % Get the pattern of initialP using either GetPattern or
            % GetMatrixData
            initialP_Pattern = initialP.GetPattern();
            if (isempty(initialP_Pattern))
                initialP_Pattern = spones(initialP.GetMatrixData());

                % Check if initialP might have some true zeros
                % due to zeros in the nullspace vectors and issue a warning
                Nullspace = FineLevel.Get('NullSpace');
                if (nnz(Nullspace==0) > 0)
                    warning('DecoalescedGraph_PatternFactory2: Pattern might be shrunk due to true zeros in the matrix initialP');
                end
            end

            % Pattern = A*P
            Pattern = spones(Amatrix.GetMatrixData()) * initialP_Pattern;
        end
    end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Protected methods                                                         %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    end % protected methods
end