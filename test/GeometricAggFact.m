%TODO: create a real AggInfo structure

% A factory which creates aggregates in the structure AggInfo
%
classdef GeometricAggFact < NeedsObject
    properties (Access = private)
        GraphName_ = 'Graph'
%         reUseAggregates_ = false        
    end
    methods
        function [this] = GeometricAggFact(obj)
            % Create aggregation factory which builds geometric 2D aggregates
            % corresponding to rectangular subdomains. In particular,
            % aggregates are centered at domain cross-points. Points within
            % subdomains are assigned to the closest cross-point.  Thus, away
            % from boundaries, each aggregate takes 1/4 of the 4 subdomains
            % which meet at the cross-point while a subdomain in the corner of the
            % entire domain has only 1 cross-point and so all of its points are
            % assigned to this cross-point.
            %
            % Copy constructor
            if nargin == 1 && isa(obj, class(this)), this.Copy_(obj,[]); return; end
            %
        end
        function SetGraphName(this, name)
            this.GraphName_ = name;
        end
        
        
        %       function SetNeeds(this, FineLevel, CoarseLevel)
        %         FineLevel.Request(this.GraphName_)
        %         FineLevel.Request('xcoords');
        %         FineLevel.Request('ycoords');
        %         FineLevel.Request('DomainList');
        %         FineLevel.Request('CrossPointList');
        %         FineLevel.Request('Bcs');
        %       end
        
%         function [ToF] = ReUseAggregates(this, ToF)
%             if varexist('ToF'),
%                 ToFold = this.reUseAggregates_;
%                 this.reUseAggregates_ = ToF;
%                 ToF = ToFold;
%             else ToF = this.reUseAggregates_;
%             end
%         end
        
        function SetNeeds(this, CurrentLevel)
            % Obtain any cross factory specifications
               if CurrentLevel.IsAvailable('Aggregates',this), return; end; % TODO: fix me for use with factory
%             if ~this.ReUseAggregates()
                CurrentLevel.Request(this.GraphName_);
                CurrentLevel.Request('xcoords');
                CurrentLevel.Request('ycoords');
                CurrentLevel.Request('DomainList');
                CurrentLevel.Request('CrossPointList');
                CurrentLevel.Request('Bcs');
%             end
        end
        
        function flag = Build(this, CurrentLevel, Specs)
            % take a matrix graph and build a set of aggregates
            %
            % On output:
            %     AggInfo.AggId(i)        gives the aggregate id corresponding to the
            %                             aggregate containing the ith vertex of Graph
            %
            %     AggInfo.NodesInAgg(k,j) = 1 indicates that the jth vertex is
            %                             assigned to the kth aggregate.
            %
            fprintf('===============> Entering GeometricAggFact\n');
            totalID=tic;
            flag = true;
            
            if CurrentLevel.IsAvailable('Aggregates',this), return; end; % TODO: fix me for use with factory
            
            Graph             = CurrentLevel.Get(this.GraphName_);
            xcoords           = CurrentLevel.Get('xcoords');
            ycoords           = CurrentLevel.Get('ycoords');
            DomainList        = CurrentLevel.Get('DomainList');
            CrossPointList    = CurrentLevel.Get('CrossPointList');
            BoundaryConditions= CurrentLevel.Get('Bcs');
            
            % figure out the bounding box of the aggregate regions. Normally, this
            % would be defined using the midpoint of each subdomain connected
            % to a cross-point. However, domains which touch boundaries are not
            % shared between aggregates, thus we must include all the domains points
            % in that coordinate direction within an entire aggregate.
            %
            %   minvalue(i,k) should correspond to the lowest bounding box value
            %   of the ith crosspont in the kth coordinate direction.
            %
            minvalue = ones(CrossPointList.NCrossPoints,DomainList.dim);
            maxvalue = ones(CrossPointList.NCrossPoints,DomainList.dim);
            minvalue(:,1) = max(xcoords) + 1;
            maxvalue(:,1) = min(xcoords) - 1;
            if DomainList.dim > 1,
                minvalue(:,2) = max(ycoords) + 1;
                maxvalue(:,2) = min(ycoords) - 1;
            end
            for i=1:CrossPointList.NCrossPoints
                mydomains = CrossPointList.Domains(i,:);
                for k=1:DomainList.dim
                    for j=1:length(mydomains)
                        mymin=(DomainList.LowerLeftCorner(mydomains(j),k) + ...
                            DomainList.UpperRightCorner(mydomains(j),k))/2;
                        mymax = mymin;
                        
                        if DomainList.Bcs(mydomains(j),2*k-1) == 'b',
                            mymin = DomainList.LowerLeftCorner(mydomains(j),k);
                        end
                        if DomainList.Bcs(mydomains(j),2*k) == 'b',
                            mymax = DomainList.UpperRightCorner(mydomains(j),k);
                        end
                        minvalue(i,k) = min(minvalue(i,k),mymin);
                        maxvalue(i,k) = max(maxvalue(i,k),mymax);
                    end
                end
            end
            % which aggregate owns each point
            AggregateId = zeros(length(xcoords),1);
            %  could be slow
            %Naggs = CrossPointList.NCrossPoints;
            n     = length(xcoords);
            droppedCrossPoints = sparse(CrossPointList.NCrossPoints,1);
            Naggs = 0;
            ii=1;
            emptyAggs = 0;
            tic
            fprintf('Setting up AggregateId structure:    ');
            aggSizes = zeros(CrossPointList.NCrossPoints,1);
            for j=1:CrossPointList.NCrossPoints
                inside = ones(n,1);
                inside(xcoords  <= minvalue(j,1)) = 0;
                inside(xcoords  >  maxvalue(j,1)) = 0;
                if DomainList.dim > 1,
                    inside(ycoords <= minvalue(j,2)) = 0;
                    inside(ycoords >  maxvalue(j,2)) = 0;
                end
                myAggSize = nnz(inside);
                if myAggSize > 0
                  Naggs = Naggs+1;
                  aggSizes(Naggs) = myAggSize;
                  AggregateId(inside>0) = Naggs;
                else
                  emptyAggs = emptyAggs + 1;
                  droppedCrossPoints(j) = 1;
                end
            end %for j=1:CrossPointList.NCrossPoints
            toc
            if emptyAggs > 0
              fprintf('WARNING: found %d empty aggregates\n',emptyAggs);
            end
            aggSizes = aggSizes(1:Naggs,:);
            maxAggSize = max(aggSizes);
            minAggSize = min(aggSizes);
            nbins = 10;
            delta = floor((maxAggSize - minAggSize) / nbins);
            if delta == 0, delta = 1; end
            bins = zeros(1:nbins);
            minmax = zeros(nbins,2);
            for ii=1:nbins,
              minmax(ii,:) = [minAggSize+delta*(ii-1)   minAggSize+delta*ii];
            end
            minmax(nbins,2) = maxAggSize;
            for ii=1:Naggs
               mybin = ceil((aggSizes(ii) - minAggSize) / delta);
               if mybin > nbins , mybin = nbins; end
               if mybin == 0, mybin = 1; end
               if aggSizes(ii) < minmax(mybin,1) || aggSizes(ii) > minmax(mybin,2)
                 error('wrong bin!');
               end
               bins(mybin) = bins(mybin) + 1;
            end
            fprintf('#bins = %d, min agg size = %d, max agg = %d\n',nbins,minAggSize,maxAggSize);
            fprintf('bin        minPts             maxPts     #aggs\n');
            for ii=1:nbins,
              fprintf('%2d   (%9d      to %9d)  :   %5d \n',ii,minmax(ii,1),minmax(ii,2),bins(ii));
            end

            %dbstack; keyboard
            RowMap = Graph.GetRowMap();

            IJV = zeros(RowMap.NDOFs(),3);
            jj=1;
            for i=1:RowMap.NDOFs()
              if (AggregateId(i) > -1)
                IJV(jj,:) = [AggregateId(i) i 1];
                jj=jj+1;
              end
            end
            NodesInAggregate = spconvert(IJV);
            clear IJV
            AggInfo.AggId      = AggregateId;
            AggInfo.NodesInAgg = NodesInAggregate;

            % set up vectors for faster lookup of nodes in a given aggregate
            [AggregateId,s] = sort(AggregateId);
            NodesInAggregate = [1:RowMap.NDOFs()];
            NodesInAggregate = NodesInAggregate(s);
            indAgg = zeros(Naggs+1,1);
            indAgg(1) = 1;
            cnt=1;
            for ii=1:length(AggregateId)
              if AggregateId(ii) ~= cnt
                cnt=cnt+1;
                indAgg(cnt) = ii;
              end
            end
            indAgg(end) = length(AggregateId)+1;
            
            % pick a point which seems far from aggregate boundaries (both
            % internal and real domain boundaries).
            roots = zeros(Naggs,1);
            A = Graph.GetMatrixData;
            n = RowMap.NDOFs();
            tvec = zeros(n,1);
            tvec(BoundaryConditions) = 1.;
            junk = zeros(n,1);

            tic

            fprintf('Finding aggregate centers:  ');

            % set up fast access to rows of A
            % sI and sJ are row/col indices of nonzeros in A
            % sI is in contiguous ascending order
            % indI(k):indI(k+1)-1 are the indices of sI and sJ for row k
            [I,J] = find(A);
            [sI,s] = sort(I); sJ = J(s);
            clear I J
            indI = zeros(size(A,1)+1,1);
            ctr=1;
            indI(1)=1;
            for ii=1:length(sI)
              if sI(ii) ~= ctr  
                while ctr<sI(ii)
                  ctr = ctr+1;
                  indI(ctr) = ii;
                end  
              end  
            end  
            nnzPerRow = spones(A) * ones(size(A,2),1);
            
            if ctr < length(indI)
              for ii=ctr+1:length(indI)
                indI(ii) = length(sI)+1;
              end
            end

            for i=1:Naggs

                %rows in current aggregate
                RowInds = NodesInAggregate(indAgg(i):indAgg(i+1)-1);
                aggSize = length(RowInds);

                %preallocate space for submatrix
                nnzInSubMatrix=0;
                for jj=1:aggSize
                  nnzInSubMatrix = nnzInSubMatrix + nnzPerRow(RowInds(jj));
                end
                subI = zeros(nnzInSubMatrix,1);
                subJ = zeros(nnzInSubMatrix,1);

                %do the moral equivalent of subA=spones(A(RowInds,:))
                cnt = 1;
                for jj=1:aggSize
                  s = indI(RowInds(jj));
                  e = indI(RowInds(jj)+1);
                  nnzInMyRow = e - s;
                  if nnzInMyRow>0
                    subI(cnt:cnt+nnzInMyRow-1) = jj;
                    subJ(cnt:cnt+nnzInMyRow-1) = sJ(s:e-1);
                  end
                  cnt = cnt + nnzInMyRow;
                end
                subA = sparse(subI,subJ,1,aggSize,size(A,2));

                ColInds = unique(subJ);
                junk(ColInds) = 1; junk(RowInds) = 0;
                border = ColInds(junk(ColInds)>0);
                junk(border) = 0;
                tvec(border) = 1;
                At = subA*tvec;
                [ii,jj] = min(At);
                % keep doing matvecs until 1's on boundary have propagated to all nodes
                while ii == 0,
                    tvec(RowInds) = At;
                    At = subA*tvec;
                    [ii,jj] = min(At); %sparse(tvec)',
                end
                roots(i) = RowInds(jj);
                tvec(ColInds) = 0;
            end;
            toc
            AggInfo.Roots = roots;
            CurrentLevel.Set('Aggregates', AggInfo, this);
            CurrentLevel.Keep('DroppedCrossPoints', droppedCrossPoints );
            CurrentLevel.Set('DroppedCrossPoints', droppedCrossPoints);
            
            CurrentLevel.Release(this.GraphName_)
            CurrentLevel.Release('xcoords');
            CurrentLevel.Release('ycoords');
            CurrentLevel.Release('DomainList');
            CurrentLevel.Release('CrossPointList');
            CurrentLevel.Release('Bcs');

            fprintf('Done building aggregates:    '); toc(totalID)
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
            [cmd, data, mc] = this.CopyCmd_(src,mc);
            eval(cmd);
        end
        
    end % methods
    
end
