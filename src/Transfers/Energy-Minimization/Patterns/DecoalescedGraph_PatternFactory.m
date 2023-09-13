%% DecoalescedGraph_PatternFactory
% concrete implementation of sparsity pattern factory that uses a
% decoalesced version of the graph of matrix A

% FIXME:
% Even if this class is suppose to be the old 'FilteredAP' option, it does not
% seems to use FilteredA :(. See also DecoalescedGraph_PatternFactory2

classdef DecoalescedGraph_PatternFactory < PatternFactory
    properties (Access = private)
        options_ = [];
        CoalesceFact_ = [];
    end

    methods
        function [this] = DecoalescedGraph_PatternFactory(filterType, CoalesceFact, InitPFact, options)
            % DecoalescedGraph_PatternFactory constructor
            %
            %   SYNTAX obj = DecoalescedGraph_PatternFactory(filterType,options)
            %
            %     filterType   - optional filtering (default = [])
            %     CoalesceFact - CoalesceDropFactory (that generated Graph
            %                    of matrix)
            %     InitPFact    - PFactory object for initial guess 'Ptent'
            %     options      - options stored in options_ (used for
            %                    fattening)
            %

            % copy constructor
            if nargin == 1 && isa(filterType, class(this)), this.Copy_(filterType,[]); return; end

            if varexist('options'), this.options_ = options; end

            if varexist('CoalesceFact'), this.CoalesceFact_ = CoalesceFact;
            else error('DecoalescedGraph_Pattern: no CoalesceDropFactory given!\n'); end

            if varexist('InitPFact'), this.InitPFact_ = InitPFact;
            else error('DecoalescedGraph_Pattern: no InitPFact given!\n'); end

            if varexist('filterType') this.filter_ = filterType;
            else this.filter_ = []; end

            this.type_ = 'DecoalescedGraph';
        end
        function SetNeeds(this, FineLevel, CoarseLevel)
            % FIXME: does not seems right
            % Ex1:  if isempty(PatternFact.GetInitialPFactory()), Factory
            %       set in BuildPattern() but no need requested.
            % Ex2:  PtentForFilter not requested if needed
            % Obtain any cross factory specifications

            SetNeeds@PatternFactory(this,FineLevel,CoarseLevel);

            FineLevel.Request('Graph',this.CoalesceFact_);
            CoarseLevel.Request('P', this.InitPFact_);

            if ~isfield(this.options_,'NCoarseDofPerNode')
                FineLevel.Request('NullSpace');
            end
        end
        function [Pattern] = BuildPattern(this, FineLevel, CoarseLevel)
            %BUILD Build a pattern based upon smoothed aggregation.
            %
            %   SYNTAX   Pattern = obj.Build(FineLevel, CoarseLevel)
            %
            %     FineLevel        - FineLevel object (input)
            %                        provides access:
            %                           - matrix A
            %                           - fine level nullspace Bzero
            %                           - aggregation info AggInfo
            %                           - fine level ID
            %     CoarseLevel      - CoarseLevel object (input)
            %                        provides access:
            %                           - coarse level nullspace Bone
            %                           - tentative prolongator initialP
            %     Pattern          - binary sparsity pattern (output)


            if nargout() ~= 1, error('Build method returns 1 argument'); end

            initialP = this.GetInitialP(FineLevel,CoarseLevel);

            if ~FineLevel.IsAvailable('Graph',this.CoalesceFact_), error('no graph available for DecoalescedGraph_PatternFactory'); end;

            Graph = FineLevel.Get('Graph',this.CoalesceFact_);
            FineLevel.Release('Graph',this.CoalesceFact_);
            CFact = CoalesceDropFactory();
            CFact.SetPostDropSpecifications([],[]);
            CFact.SetPreDropSpecifications ([],[]);
            Pnode = CFact.Build_(initialP);
            ftemp = FineLevel.BuildMe();
            ftemp.Set('A', Graph);
            ctemp = CoarseLevel.BuildMe();
            ctemp.Request('P', this);
            ctemp.Request('P', this);
            ctemp.Set('P', Pnode, this);

            PatternFact = AP_PatternFactory();
            ctemp.Request('Ppattern', PatternFact);
            % link pattern factories together!
            if isempty(PatternFact.GetInitialPFactory())
                PatternFact.SetInitialPFactory(this);
            end

            PatternFact.Build(ftemp,ctemp);
            APnode = ctemp.Get('Ppattern',PatternFact);
            ctemp.Release('Ppattern',PatternFact);

            if ~isfield(this.options_,'NCoarseDofPerNode')
                NS = FineLevel.Get('NullSpace');
                FineLevel.Release('NullSpace');
                Ncpn = size(NS,2);
            else Ncpn = this.options_.NCoarseDofPerNode; end;

            Pattern = this.DecoalesceMatrix(APnode,Ncpn,Ncpn);
        end

        function [initialP] = GetInitialP(this,FineLevel,CoarseLevel)
            %GETINITIALP communication layer function between levels and
            %this factory
            %
            %   SYNTAX   Pinitial = obj.GetInitialP(FineLevel, CoarseLevel)
            %
            %     FineLevel        - FineLevel object (input)
            %                        provides access:
            %                           - matrix A
            %                           - fine level nullspace Bzero
            %                           - aggregation info AggInfo
            %                           - fine level ID
            %     CoarseLevel      - CoarseLevel object (input)
            %                        provides access:
            %                           - coarse level nullspace Bone
            %                           - tentative prolongator initialP

            if ~CoarseLevel.IsAvailable('P', this.InitPFact_)
                this.InitPFact_.BuildP(FineLevel,CoarseLevel);   % bad: this overwrites CoarseLevel.'P'
            end
            initialP = CoarseLevel.Get('P', this.InitPFact_);
            CoarseLevel.Release('P', this.InitPFact_);

            if ~isempty(this.filter_)
                % provide information for filter
                CoarseLevel.Set('PtentForFilter',initialP);
            end
        end %GetInitialP()
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
        function Matrix2=DecoalesceMatrix(this,Matrix,Nr,Nc)
            % DECOALESCEMATRIX
            %
            % SYNTAX M=obj.DecoalesceMatrix(Matrix,Rb,Cb)
            %
            %   Matrix  - "nodal" matrix
            %   Rb      - number of rows per block
            %   Cb      - number of columns per block
            %   Matrix2 - A matrix with a RbxCb block for each non-zero in Matrix
            %              The entries of matrix2 will be all ones.
            [R,C,V]=find(Matrix);
            Blk=Nr*Nc;
            NNZ=length(R);
            Rblk=reshape(repmat(1:Nr,Nc,1),Blk,1);
            Cblk=reshape(repmat(1:Nc,Nr,1)',Blk,1);

            NewR=reshape(repmat((R-1)*Nr,1,Blk)',NNZ*Blk,1)+reshape(repmat(Rblk',NNZ,1)',NNZ*Blk,1);
            NewC=reshape(repmat((C-1)*Nc,1,Blk)',NNZ*Blk,1)+reshape(repmat(Cblk',NNZ,1)',NNZ*Blk,1);
            V=ones(NNZ*Blk,1);
            Matrix2=sparse(NewR,NewC,V);
        end
    end % protected methods
end