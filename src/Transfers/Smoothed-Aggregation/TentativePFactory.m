%% TentativePFactory
% class implements PFactory interface and creates a tentative prolongation
% operator.
%%

classdef TentativePFactory < PFactory
    % class implements PFactory interface and creates a tentative prolongation
    % operator.
    % class also provides static functions BuildAggregates and
    % MakeTentative that can be used by other prolongation operator classes
    properties (Access=protected)
        AggFact_             = [];
        QR_                  = true; % use QR decomposition for improving nullspace information per default
    end

    methods
        % default parameters for constructor:
        % AggregationFactory = AggregationFactory() for standard aggregation
        % techniques
        function [this] = TentativePFactory(CoalesceFact, AggregationFact)
            % TODO: copy constructor

            if nargin == 1,
                if nargin == 1 && isa(CoalesceFact, class(this)), this.Copy_(CoalesceFact,[]); return; end
                warning('TentativePFactory with only one parameter = copy constructor');
            end % TMP

            if varexist('AggregationFact'), this.AggFact_ = AggregationFact;
            %else this.AggFact_ = 'default'; end this is not working
            else this.AggFact_ = AggregationFactory(); end

            % for compatibility reasons
            if varexist('CoalesceFact'), this.AggFact_.SetCoalesceFactory(CoalesceFact); end

            if this.GetOutputLevel() > 5
                fprintf('TentativePFactory constructor\n');
            end
        end

        function TentativeWithQR(this, value)
            % Specify whether orthogonalization performed when forming tentative prolongator
            this.QR_ = value;
        end

        function SetNeeds(this, FineLevel, CoarseLevel)
          if this.CheckForReUsableP(FineLevel,CoarseLevel,this) == true, return; end

            % translate AggFact_ (e.g. 'default' to default factory)
            this.AggFact_ = FineLevel.InterpretHandle('Aggregates',this.AggFact_);
            if ~FineLevel.IsRequested('Aggregates', this.AggFact_) & ...
               ~FineLevel.IsAvailable('Aggregates', this.AggFact_)
                if ~ismethod(this.AggFact_,'SetNeeds'),
                    error('TentativePFactory.SetNeeds: AggFact_ must be a factory\n');
                end
                this.AggFact_.SetNeeds(FineLevel);
            end
            FineLevel.Request('Aggregates',this.AggFact_);

            FineLevel.Request('NullSpace');
        end

        function [ToF] = SupportsRestrictionMode(this) % should be static
            ToF = false;
        end

        function flag = Build(this,FineLevel,CoarseLevel)
            flag = true;

            %% 0) check if prolongator is already built -> reuse it if possible
            if this.CheckForReUsableP(FineLevel,CoarseLevel,this) == true, return; end

            % Check out system matrix, that is used for aggregation
            Amat = FineLevel.Get('A');

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
            AggInfo = FineLevel.Get('Aggregates',this.AggFact_);
            FineLevel.Release('Aggregates',this.AggFact_);

            % Build tentative prolongator P0
            [P,cnull] = TentativePFactory.MakeTentative(AggInfo,Amat,NS,this.QR_,this.GetOutputLevel());

% FUTURE??
%            % instead of above line maybe we should use something like
%            CoarseNSFact = SaCoarseNSFactory();
%            CoarseNSFact.TentativeWithQR(this.QR_);
%            [P,cnull] = CoarseNSFact.Build(AggInfo,Amat,NS,[],this.GetOutputLevel());
%            clear CoarseNSFact;

            % Store information

            CoarseLevel.Set('P', P, this);
            CoarseLevel.Set('NullSpace', cnull);

        end
    end

    methods (Static = true)


        function  [P,CoarseNull] = MakeTentative(AggInfo,Amat,nullspace,QROnOrOff,OutputLevel)
        % Make a tentative prolongator optionally use QR to orthogonalize columns.

            %
            % Currently making an assumption
            %    1) The line 'LocalCoarseNull = R(1:nulldim,1:nulldim);' is
            %       basically assuming that localP has at least nulldim rows.
            %
            % I'm also not completely sure what happens if there are not
            % at least nulldim linearly independent columns in localP. I
            % think the QR still gives a  square full rank matrix Q ... so
            % probably everything is okay. The R is rank deficient and this
            % might eventually cause a rank deficient P on the next coarsest
            % level, but I haven't checked this out. I'm adding some random
            % stuff in empty columns (in MakeNoQRTentative.m) to hopefully
            % avoid any serious problems.
            %
            % It might also be interesting to see if 'VarBlk' column
            % maps could be generated for P in the case that we have
            % some empty (or degenerate columns). Right now MakeNoQRTentative.m
            % jams in some random numbers when it detects that we have an
            % empty column so that P remains a 'ConstBlk' column map matrix.
            %

            P = TentativePFactory.MakeNoQRTentative(AggInfo,Amat,nullspace,OutputLevel);

            AggId      = AggInfo.AggId;
            NodesInAgg = AggInfo.NodesInAgg;
            Naggs   = max(AggId);
            nulldim = size(nullspace,2);
            MatData = P.GetMatrixData();
            CoarseNull = sparse(Naggs*nulldim,nulldim);

            % set up vectors for faster lookup of nodes in a given aggregate
            [Aggs,Nodes] = find(NodesInAgg); % Nodes are in ascending order
            [Aggs,s] = sort(Aggs);
            Nodes = Nodes(s);

            indAgg = zeros(Naggs+1,1);
            indAgg(1) = 1;
            cnt=1;
            for ii=1:length(Aggs)
              if Aggs(ii) ~= cnt
                cnt=cnt+1;
                indAgg(cnt) = ii;
              end
            end
            indAgg(end) = length(Aggs)+1;

            count = 1;
            RowMap = Amat.GetRowMap();
            for i=1:Naggs
              BlkRows = Nodes(indAgg(i):indAgg(i+1)-1);
              PtRowList = Node2DOF(BlkRows,RowMap);
              localP = MatData(PtRowList,count:count+nulldim-1);
              if QROnOrOff,
                [Q,R] = qr(full(localP));
              else
                Q = localP;
                R = speye(size(localP,2),size(localP,2));
              end
              if size(localP,1) < nulldim,
                fprintf('MakeTentative:: Not enough DOFs for QR\n');
                keyboard;
              end
              MatData(PtRowList,count:count+nulldim-1) = Q(:,1:nulldim);
              LocalCoarseNull = R(1:nulldim,1:nulldim);
              CoarseNull((i-1)*nulldim+1:i*nulldim,1:nulldim) = LocalCoarseNull;
              count = count + nulldim;
            end

            %fprintf('An expensive check on coarse nullspace => %e\n',...
            %norm(nullspace-MatData*pinv(full(MatData'*MatData))*MatData'*nullspace,'fro'));

            P.SetMatrixData(MatData);
        end %MakeTentative()

        function  [P] = MakeNoQRTentative(AggInfo,Amat,nullspace,OutputLevel)
        % Make a tentative prolongator without using QR to orthogonalize columns.

            %
            % Note: we might also have a constant block size in one dimension
            % while a variable block size in the other.
            %
            % It might be nice if we could generate variable block (in
            % column direction) P's in the case where we have a column
            % of all zeros (or if it is linearly dependent). That is,
            % if designed properly, perhaps it is not so hard to
            % have variable block stuff floating around on all levels
            % in the MG hierarchy as opposed to just the finest level.
            %
            AggregateId      = AggInfo.AggId;
            NodesInAgg = AggInfo.NodesInAgg;
            RowMap = Amat.GetRowMap();
            Naggs   = max(AggregateId);

            % set up vectors for faster lookup of nodes in a given aggregate
            [Aggs,Nodes] = find(NodesInAgg); % Nodes are in ascending order
            [Aggs,s] = sort(Aggs);
            Nodes = Nodes(s);

            indAgg = zeros(Naggs+1,1);
            indAgg(1) = 1;
            cnt=1;
            for ii=1:length(Aggs)
              if Aggs(ii) ~= cnt
                cnt=cnt+1;
                indAgg(cnt) = ii;
              end
            end
            indAgg(end) = length(Aggs)+1;

            nulldim = size(nullspace,2);
            MatData    = sparse(RowMap.NDOFs(),Naggs*nulldim);
            MatPattern = sparse(RowMap.NDOFs(),Naggs*nulldim); % MATLAB specific: pattern with "true" zeros

            warnings = 0;
            count = 1;
            PrevStream = RandStream.setGlobalStream(RandStream.create('mrg32k3a','NumStreams',1));
            for i=1:Naggs
             BlkRows = Nodes(indAgg(i):indAgg(i+1)-1);
             PtRowList = Node2DOF(BlkRows,RowMap);
             for j=1:nulldim
                MatData(PtRowList,count)    = nullspace(PtRowList,j);
                MatPattern(PtRowList,count) = 1;
                if nnz(MatData(PtRowList,count)) == 0,
                   MatData(PtRowList,count) = rand(length(PtRowList),1);
                   if warnings == 0,
                      if OutputLevel > 5,
                         fprintf('MakeNoQRTentative Warning: Inserting random numbers to avoid empty columns\n');
                      end
                   end
                   warnings = 1;
                end
                count = count + 1;
             end
            end

            P = Operator(MatData, RowMap, Map(Naggs,nulldim), @MatlabApply);
            P.SetPattern(MatPattern);  % MATLAB specific: pattern with "true" zeros

            RandStream.setGlobalStream(PrevStream);
        end %MakeNoQRTentative()

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
