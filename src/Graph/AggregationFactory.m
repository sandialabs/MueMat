%TODO: create a real AggInfo structure

% A factory which creates aggregates in the structure AggInfo
%
classdef AggregationFactory < VerboseObject
   properties (Access = private)
      Algorithm_
      TargetSize_
      DomainPointsPerDim_
      AggPointsPerDim_
      chop_
      mloptions_
      RemoveAggsThisSizeOrLess_

      CoalesceFact_
      reUseAggregates_ = false
      GraphName_ = 'Graph'
   end
   methods
      function [this] = AggregationFactory(CoalesceFact)
       % Copy constructor
       if nargin == 1 && isa(CoalesceFact, class(this)), this.Copy_(CoalesceFact,[]); return; end
       %

       % CoalesceFactory = CoalesceDropFactory() for basic amalgamation and
       % dropping of coarse grid level problems
       if varexist('CoalesceFact'), this.CoalesceFact_ = CoalesceFact;
       else this.CoalesceFact_ = CoalesceDropFactory(); end

         % create aggregation factory
         this.Algorithm_ = 'graph';
         this.TargetSize_         = 14;
         this.DomainPointsPerDim_ = [30 30];
         this.AggPointsPerDim_    = [3 3];
         this.RemoveAggsThisSizeOrLess_ = 0;
      end

      function SetCoalesceFactory(this,CoalesceFact)
        this.CoalesceFact_ = CoalesceFact;
      end

      function SetMLOptions(this,varargin)
      % set options for ML's aggregation
        this.mloptions_=varargin;
      end
      function SetAlgorithm(this,type)
      % set aggregation algorithm: graph, rcm, square,cubic, ml, order, rectangle, hypergraph
         if ( ~strcmp(type,'graph' ) && ~strcmp(type,'rcm')   && ...
              ~strcmp(type,'square') && ~strcmp(type,'cubic') && ...
              ~strcmp(type,'ml')     && ~strcmp(type,'Uniform1D') && ...
              ~strcmp(type,'order')  && ~strcmp(type,'rectangle')) && ...
              ~strcmp(type,'hypergraph')
            fprintf('Aggregation Factory: Unknown Aggregation Algorithm\n');
            keyboard;
         end
         if strcmp(type,'ml') && ~exist('ml','file')
           warning('AggregationFactory.SetAlgorithm(): MLMEX not found. Cannot use ML Aggregation.');
           return;
         end

         this.Algorithm_ = type;
      end
      function SetTargetSize(this, size)
      % set target aggregate size for rcm and order algorithms
         this.TargetSize_ = size;
         if ( ~strcmp(this.Algorithm_,'rcm' ) && ...
              ~strcmp(this.Algorithm_,'order') && ...
              ~strcmp(this.Algorithm_,'hypergraph') )
            fprintf('Aggregation Factory Warning: TargetSize will be ignored by this algorithm!\n');
         end
      end
      function SetDomainPtsPerDim(this, PpdArray)
      % set # of points in each coordinate dimension for order & rectangle algorithms. 1D, 2D, or 3D is decided by length(PpdArray).
         this.DomainPointsPerDim_ = PpdArray;
         if ~strcmp(this.Algorithm_,'order')&&~strcmp(this.Algorithm_,'rectangle' )&&~strcmp(this.Algorithm_,'Uniform1D' ),
            fprintf('Aggregation Factory Warning: DomainPtsPerDim only used with order or rectangle algorithms!\n');
         end
      end
      function SetAggPtsPerDim(this, PpdArray)
      % set # of agg points in each coordinate dimension for rectangle algorithm
         this.AggPointsPerDim_ = PpdArray;
         if ~strcmp(this.Algorithm_,'rectangle' ),
            fprintf('Aggregation Factory Warning: AggPtsPerDim only used with rectangle algorithm!\n');
         end
      end
      function SetChop(this,chop)
        this.chop_=chop;
      end
      function SetRemoveAggsThisSizeOrLess(this, MinSize)
         this.RemoveAggsThisSizeOrLess_ = MinSize;
      end

%       function [ToF] = ReUseAggregates(this, ToF)
%         if varexist('ToF'),
%           ToFold = this.reUseAggregates_;
%           this.reUseAggregates_ = ToF;
%           ToF = ToFold;
%         else ToF = this.reUseAggregates_;
%         end
%       end

      function [mloptions] = GetMLOptions(this)
        mloptions = this.mloptions_;
      end
      function [type] = GetAlgorithm(this)
        type = this.Algorithm_;
      end
      function [size] = GetTargetSize(this)
         size = this.TargetSize_;
        if ( ~strcmp(this.Algorithm_,'rcm' ) && ...
              ~strcmp(this.Algorithm_,'order') && ...
              ~strcmp(this.Algorithm_,'hypergraph') )
          fprintf('Aggregation Factory Warning: TargetSize will be ignored by this algorithm!\n');
         end
      end
      function [PpdArray] = GetDomainPtsPerDim(this)
      % get # of points in each coordinate dimension for order & rectangle algorithms. 1D, 2D, or 3D is decided by length(PpdArray).
         PpdArray = this.DomainPointsPerDim_;
         if ~strcmp(this.Algorithm_,'order')&&~strcmp(this.Algorithm_,'rectangle' ),
            fprintf('Aggregation Factory Warning: DomainPtsPerDim only used with order or rectangle algorithms!\n');
         end
      end
      function [PpdArray] = GetAggPtsPerDim(this)
      % get # of agg points in each coordinate dimension for rectangle algorithm
         PpdArray = this.AggPointsPerDim_;
         if ~strcmp(this.Algorithm_,'rectangle' ),
            fprintf('Aggregation Factory Warning: AggPtsPerDim only used with rectangle algorithm!\n');
         end
      end
      function [chop] = GetChop(this)
        chop=this.chop_;
      end
      function [MinSize] = GetRemoveAggsThisSizeOrLess(this)
         MinSize = this.RemoveAggsThisSizeOrLess_;
      end


      function SetGraphName(this, name)
       % indicates use of filtered version of A (small entries dropped) within prolognator smoothing
       this.GraphName_ = name;
      end

      function SetNeeds(this, CurrentLevel)
        % Obtain any cross factory specifications
          if CurrentLevel.IsAvailable('Aggregates',this), return; end;

          if ~isempty(this.CoalesceFact_), this.CoalesceFact_.SetNeeds(CurrentLevel); end;
          CurrentLevel.Request(this.GraphName_,this.CoalesceFact_);
        %end
      end

      function flag = Build(this,CurrentLevel,Specs)
        flag = true;

        if CurrentLevel.IsAvailable('Aggregates',this), return; end;

        if ~isempty(this.CoalesceFact_)
          this.CoalesceFact_.Build(CurrentLevel);
        elseif CurrentLevel.IsAvailable(this.GraphName_,this.CoalesceFact_)
          fprintf('Aggregation: using a user defined graph (s) (skipped call to CoalesceFact)\n', this.GraphName_);
        end

        Graph = CurrentLevel.Get(this.GraphName_,this.CoalesceFact_);
        CurrentLevel.Release(this.GraphName_,this.CoalesceFact_);
        AggInfo = this.Build_(Graph);

        CurrentLevel.Set('Aggregates', AggInfo, this);

      end

      function [AggInfo] = Build_(this, Graph, Specs)
      % take a matrix graph and build a set of aggregates
      %
      % On output:
      %     AggInfo.AggId(i)        gives the aggregate id corresponding to the
      %                             aggregate containing the ith vertex of Graph
      %
      %     AggInfo.NodesInAgg(k,j) = 1 indicates that the jth vertex is
      %                             assigned to the kth aggregate.
      %
        AggInfo = this.Basic(Graph);
      end

      function [AggInfo] = Basic(this,Amat)
      %    Aggregate Amat using an algorithm specified by this.Algorithm_
      %
      %       Algorithm='graph'   Smoothed aggregation graph based aggregation.
      %
      %       Algorithm='rcm'     Simple RCM aggregation. Default size
      %                           of aggregates is 14. Can be altered by
      %                           SetTargetSize().
      %
      %       Algorithm='square'  Make square aggregates. I believe this makes
      %                           3x3's and assumes domain is square.
      %
      %       Algorithm='cubic'   Make cubic aggregates. I believe this makes
      %                           3x3x3's and assumes cube domain.
      %
      %       Algorithm='rectangle' Make rectangular aggregates assuming domain
      %                             is rectangular. Defaults: 30x30 2d domain
      %                             to be tiled with 3x3 aggregates. Defaults
      %                             can be changed via
      %                                SetDomainPtsPerDim(PpdArray)
      %                                SetAggPtsPerDim(PpdArray)
      %
      %       Algorithm='order'   Ordering-based aggregation for logical
      %                           rectangular meshes in an arbitrary number of
      %                           dimensions.  Aggressiveness of coarsening
      %                           set by the SetTargetSize(). Domain size
      %                           set via SetDomainPtsPerDim(PpdArray)
      %
      %       Algorithm='ml'      Use ML's aggregation routine,
      %                           faking root points by grabbing
      %                           the highest degree node in the aggregate.
      %
      %       Algorithm='hypergraph' Aggregate by hypergraph partitioning.
      %                              Either PaToH or Zoltan is required.
      %                              Use SetTargetSize() to set target size.
      %
      %
      %  See this.Build() documentation for more information on output.
      %
      %    Note: There is reduncancy between rectangle, square, cubic, & order.
      %    These are different aggregation algorithms from the past. They have
      %    different options/features so all are included.

      A = Amat.GetMatrixData();
      RowMap = Amat.GetRowMap();

      if strcmp(this.Algorithm_,'rectangle') || strcmp(this.Algorithm_,'order')
        dim   = length(this.DomainPointsPerDim_);
        Nnodes= Amat.GetRowMap().NNodes();
        while ( prod(this.DomainPointsPerDim_) ~= Nnodes)
          if Nnodes > prod(this.DomainPointsPerDim_),
             error('# points per dimension via SetDomainPtsPerDim() does not match total number of nodes\n');
          end
          for i=1:dim
             this.DomainPointsPerDim_(i)=this.DomainPointsPerDim_(i)/this.AggPointsPerDim_(i);
          end
        end
      end

      if strcmp(this.Algorithm_,'rcm')
         %  Simple aggregation based on rcm reordering and then a linear
         %  decompsition.
         Naggs = round(RowMap.NDOFs()/this.TargetSize_ );
         if Naggs < 1, Naggs = 1; end
         [AggregateId, perm] = rcmpart(A,Naggs);
         AggregateId= AggregateId+1; % Make 1-based.

         % find roots nodes
         nn = length(AggregateId);
         mm = max(AggregateId);
         Ptent = sparse(1:nn, AggregateId, ones(nn,1), nn, mm);
         Apattern = (A~=0);
         AP = Ptent .* (Apattern * Ptent);
         [C,I] = max(AP);
         roots = I;
      elseif strcmp(this.Algorithm_,'rectangle')
         RectAgg.Nx    = 1; RectAgg.Ny    = 1; RectAgg.Nz    = 1;
         RectAgg.Diamx = 1; RectAgg.Diamy = 1; RectAgg.Diamz = 1;
         dim   = length(this.DomainPointsPerDim_);
         if length(this.AggPointsPerDim_) ~= dim,
            fprintf('Mismatch between dimensions of domain and aggregates\n');
            keyboard;
         end
         RectAgg.Dim   = dim;
         if dim > 0, RectAgg.Nx    = this.DomainPointsPerDim_(1); end
         if dim > 1, RectAgg.Ny    = this.DomainPointsPerDim_(2); end
         if dim > 2, RectAgg.Nz    = this.DomainPointsPerDim_(3); end
         if dim > 0, RectAgg.Diamx = this.AggPointsPerDim_(1);    end
         if dim > 1, RectAgg.Diamy = this.AggPointsPerDim_(2);    end
         if dim > 2, RectAgg.Diamz = this.AggPointsPerDim_(3);    end

         [AggregateId, roots] = SquareAggs(Amat, RectAgg);
      elseif strcmp(this.Algorithm_,'square')
         [AggregateId, roots] = AggregationFactory.GraphAggregation(A,2);
      elseif strcmp(this.Algorithm_,'cubic')
         [AggregateId, roots] = AggregationFactory.GraphAggregation(A,3);
      elseif strcmp(this.Algorithm_,'order')
         [AggregateId, roots] = AggregationFactory.OrderAggregation(this.DomainPointsPerDim_, this.TargetSize_,zeros(dim,1));
      elseif strcmp(this.Algorithm_,'ml')
         if(isempty(this.mloptions_)) opts={}; else opts=this.mloptions_;end
         [AggregateId, roots] = AggregationFactory.MLAggregation(A,1,opts);
      % this stuff was ripped from something else and is not
      % tested in its current form
      elseif strcmp(this.Algorithm_,'Uniform1D')
         AggDiam = this.AggPointsPerDim_;
         nnn = RowMap.NDOFs();
         starter = (AggDiam+1)/2;
         AggregateId = ceil( (1:nnn)/AggDiam);
         roots = (starter:AggDiam:nnn);
      elseif strcmp(this.Algorithm_,'hypergraph')
         % Check if PaToH is installed
         if ~exist('PaToH','file')
           error('The PaToH library and mex file must be installed to use hypergraph aggregation');
         end
         %  Aggregation based on hypergraph partitioning.
         Naggs = round(RowMap.NDOFs()/this.TargetSize_ );
         if Naggs < 1, Naggs = 1; end
         AggregateId = PaToH(A,Naggs) + 1;
         % find roots
         roots = AggregationFactory.FindRoots(A, AggregateId);
      else % this.Algorithm_ == graph
           [AggregateId, roots] = AggregationFactory.GraphAggregation(A);
      end

      % Optional aggregate chopper - HAQ only in 2D
      if(this.chop_==true) [AggregateId,roots]=AggregationFactory.SplitAggregates(A,2,AggregateId,roots);end

      Naggs = max(AggregateId);

      IJV = zeros(RowMap.NDOFs(),3);
      jj=1;
      for i=1:RowMap.NDOFs()
        if (AggregateId(i) > -1), IJV(jj,:) = [AggregateId(i) i 1]; end
        jj=jj+1;
      end
      NodesInAggregate = spconvert(IJV);
      clear IJV

      if this.RemoveAggsThisSizeOrLess_ > 0,
         AggSizes = NodesInAggregate*ones(RowMap.NDOFs(),1);
         ThrowAway = find(AggSizes <= this.RemoveAggsThisSizeOrLess_);
         if length(ThrowAway) > 0,
            fprintf('Removing %d aggregates with |Agg_i| <= %d\n',length(ThrowAway),this.RemoveAggsThisSizeOrLess_);
            keep = ones(Naggs,1); keep(ThrowAway) = 0; keep = find(keep);
            NodesInAggregate = NodesInAggregate(keep,:);
            AggregateId = (1:length(keep))*NodesInAggregate;
            roots = roots(keep);
         end
      end

      AggInfo.AggId      = AggregateId;
      AggInfo.NodesInAgg = NodesInAggregate;
      AggInfo.Roots      = roots;
      end
   end
   methods (Static = true)

      function [agg,roots] = MLAggregation(A,ppn,mlopts)
      % Uses ML's aggregation routine to generate aggregates.  We
      % fake root nodes by finding the node associated with the
      % highest degree dof.
        fprintf('ML agg: Using %d equations per node\n',ppn);

        if ~exist('ml','file'), error('MLMEX not found'); end

        agg=double(1+ml('aggregate',A,'PDE equations',ppn,...
                        'coarse: max size',1,...
                        mlopts{:}));
        roots = AggregationFactory.FindRoots(A, agg);
      end

      function [roots] = FindRoots(A, agg, steps)
      % Find root notes given an aggregation (partition) of A.
      % By default (steps=1), pick a vertex of maximum degree.
      % For large aggregates, steps should be larger, too.
        if (nargin<3)
          steps=1;
        end
        Nagg=max(agg);
        roots = zeros(Nagg,1);
        max_agg_size=-1;
        for i=1:Nagg,
          idx=find(agg==i);
          degree= ones(length(idx),1);
          for s=1:steps
             degree=(abs(A(idx,idx))>0)*degree;
          end
          [v,j]=max(degree);
          roots(i)=idx(j);
          % roots(i)=floor(idx(j)/ppn); % ML takes ppn into account?
          if(length(idx)>max_agg_size), max_agg_size=length(idx);end
        end
        fprintf('FindRoots: Maximum aggregate size %d\n', max_agg_size);
      end


      function [myaggregate, roots] = GraphAggregation(Graph, makesquare)
      %
      % Given a graph this function returns the aggregates in the
      % usual aggregate data structure (with agg.length giving the number
      % of aggregates, and agg.info(aggnum).nodelist giving a list of the
      % nodes in aggregate number aggnum.  The algorithm used in this
      % function tries to minimize the size of each element - it does not
      % take into account how strongly connected an orphan element is in
      % deciding which neighboring aggregate to add it to.

      if ~varexist('makesquare'), makesquare = 0; end
      nnode=size(Graph,2);


      roots       = zeros(size(Graph,1),1);
      nodesInAgg  = [0];              % list of nodes already assigned to some aggregate
      myaggregate = -1*ones(nnode,1); % keeps track of which aggregate each node is in
      aggtmp.length=0;                % current number of aggregates


      if (makesquare == 0),

        didDot = false;
        dotCount = 0;
        lineCount = ceil(nnode/8000);
        fprintf('Aggregating\n%2d ',lineCount);
        lineCount = lineCount-1;
        for inode = 1:nnode
          if ( mod(inode,100) == 0)
            fprintf(1,'.');
            dotCount = dotCount+1;
            if dotCount == 80,
              if lineCount > 0, fprintf('\n%2d ',lineCount);
              else              fprintf('\n');               end
              dotCount = 0;
              lineCount = lineCount-1;
            end
            didDot = true;
          end
          if (myaggregate(inode)==-1)
            % Node inode is not in an aggregate yet.
            % Check if inode is adjacent to a current aggregate:  step through each
            % node jj of each aggregate to see if jj is adjacent to inode.

            adj=0; ctr=1;
            Gi_nnz = find(Graph(:,inode));

            while (~adj && (ctr<=length(Gi_nnz)))
              % adj = 1 if node Gi_nnz(ctr) is already in an aggregate
              adj=~(myaggregate(Gi_nnz(ctr))==-1);
              ctr=ctr+1;
            end

            if ~adj
              % no aggregate adjacent to inode, so make new aggregate containing inode
              aggtmp.length=aggtmp.length+1;

              % add this new aggregate to agg:
              roots(aggtmp.length) = inode;
              aggtmp.info(aggtmp.length).nodelist=Gi_nnz;

              %myaggregate(Gi_nnz) = aggtmp.length;
              myaggregate(aggtmp.info(aggtmp.length).nodelist) = aggtmp.length;

              nodesInAgg = [nodesInAgg; aggtmp.info(aggtmp.length).nodelist];
            end

          end %if isempty(find(nodesInAgg==inode))
        end %for inode
        if didDot, fprintf('\n'); end

      elseif (makesquare == 2)
         %%%%%%%%%% Square aggs for 2d %%%%%%%%%%%%%%%
         nx = sqrt(size(Graph,1));
         interior = [ -1-nx   -nx    -nx+1 ...
                       -1       0       1  ...
                      -1+nx    nx     nx+1 ];
         top      = [ -1-nx   -nx    -nx+1 ...
                       -1       0       1  ];
         topright = [ -1-nx   -nx    ...
                       -1       0     ];
         right    = [ -1-nx   -nx          ...
                       -1       0          ...
                      -1+nx    nx          ];
         AggId = 0;
         for nn=2:3:nx
            for mm=2:3:nx
               AggId = AggId + 1;
               mid = (nn-1)*nx+mm;
               nodes = interior;
               if nn==nx,
                  if mm==nx, nodes = topright;
                  else       nodes = top; end
               elseif mm==nx, nodes = right; end
               myaggregate(mid+nodes) = AggId;
               aggtmp.length = aggtmp.length + 1;
               roots(aggtmp.length) = mid;
               aggtmp.info(aggtmp.length).nodelist= (mid+nodes)';
            end
         end

      elseif (makesquare == 3)
         %%%%%%%%%% Cubic aggs for 3d %%%%%%%%%%%%%%%
         nx = floor(power(size(Graph,1),1/3)+1e-5);
         nx2 = nx*nx;
         interior2d = [ -1-nx   -nx    -nx+1 ...
                         -1       0       1  ...
                        -1+nx    nx     nx+1 ];
         top2d      = [ -1-nx   -nx    -nx+1 ...
                         -1       0       1  ];
         topright2d = [ -1-nx   -nx    ...
                         -1       0     ];
         right2d    = [ -1-nx   -nx          ...
                         -1       0          ...
                        -1+nx    nx          ];
         interior = [interior2d-nx2 interior2d interior2d+nx2];
         top      = [     top2d-nx2      top2d      top2d+nx2];
         right    = [   right2d-nx2    right2d    right2d+nx2];
         topright = [topright2d-nx2 topright2d topright2d+nx2];
         front         = [interior2d-nx2 interior2d];
         frontright    = [   right2d-nx2    right2d];
         fronttop      = [     top2d-nx2      top2d];
         fronttopright = [topright2d-nx2 topright2d];
         AggId = 0;
         for nn=2:3:nx
            for mm=2:3:nx
               for oo=2:3:nx
                  AggId = AggId + 1;
                  mid = (oo-1)*nx*nx + (nn-1)*nx+mm;
                  nodes = interior;
                  if nn==nx,
                     if mm==nx,
                        if oo==nx, nodes = fronttopright;
                        else       nodes =      topright; end
                     else
                        if oo==nx, nodes = fronttop;
                        else       nodes =      top; end
                     end
                  else
                     if mm==nx,
                        if oo==nx, nodes = frontright;
                        else       nodes =      right; end
                     else
                        if oo==nx, nodes = front; end
                     end
                  end
                  myaggregate(mid+nodes) = AggId;
                  aggtmp.length=aggtmp.length+1;
                  roots(aggtmp.length) = mid;
                  aggtmp.info(aggtmp.length).nodelist= (mid+nodes)';
               end
            end
         end
      else
        error('unknowns option for makesquare');
      end %if ~makesquare



      % So now we've tried to put all elements into aggregates, but there may be
      % some left over that are adjacent to aggregates, but that are not in an
      % aggregate

      for i=1:nnode

        if (myaggregate(i) < 0)
          % if isempty(find(nodesInAgg==i))  %i is not in any aggregate
          % Put i into the smallest aggregate that i is adjacent to.

          Gi_nnz=find(Graph(:,i));

          tnbrs = [-1];
          tcount = 1;
          for k=1:length(Gi_nnz)
            jjj = myaggregate(Gi_nnz(k));
            if ( isempty(find(jjj == tnbrs)))
              tcount = tcount +1;
              tnbrs(tcount) = jjj;
            end
          end

          % "nbrs" is a list of aggregates that i is adjacent to.
          nnbrs = tnbrs(2:tcount);

          % Determine size of each neighboring aggregate.
          nbrsizes=[];
          for k=1:length(nnbrs)
            nbrsizes = [nbrsizes;length(aggtmp.info(nnbrs(k)).nodelist)];
          end
          %smnbri is index into nnbrs of smallest neighbor
          [dummy,smnbri]=min(nbrsizes);

          if size(nnbrs,2) ~= 0
            smnbr=nnbrs(smnbri);
            aggtmp.info(smnbr).nodelist=[aggtmp.info(smnbr).nodelist;i];
            myaggregate(i) = smnbr;
          end
        end
      end

      % DBC
%       [myaggregate, aggtmp] = ...
%          AggregationFactory.RemoveAggsWithOnlyOneBlk(myaggregate, aggtmp, nnode);

      roots = roots(1:aggtmp.length);
      end %GraphAggregation()

      function [myaggregate, aggtmp] = RemoveAggsWithOnlyOneBlk(...
            myaggregate, aggtmp, nnode)
      % Build new aggtmp that has only aggregates with more than one blk

         % set myaggregate to -1, if my aggregate has only me as a member
         for i=1:nnode
            myagg = myaggregate(i);
            if myagg ~= -1 % there could be a node that is not part of any aggregate
               if size(aggtmp.info(myagg).nodelist,1) == 1
                  myaggregate(i) = -1;
               end
            end
         end

         % create new aggtmp, where aggs with only one blk are ommited
         length2 = 1;
         aggtmp2.info = [];
         for iagg=1:aggtmp.length
            numnode = size(aggtmp.info(iagg).nodelist,1);
            if (numnode > 1)
               aggtmp2.info(length2).nodelist = aggtmp.info(iagg).nodelist;
               length2 = length2 + 1;
            end
         end
         aggtmp = aggtmp2;
         aggtmp.length = size(aggtmp.info,2);

         % correct aggregate number, since we just ommited some aggs
         % somewhere in the middle
         newaggid = 1;
         for iagg=1:aggtmp.length
            numnode = size(aggtmp.info(iagg).nodelist,1);
            for inode = 1:numnode
               nodeid = aggtmp.info(iagg).nodelist(inode);
               myaggregate(nodeid) = newaggid;
            end
            newaggid = newaggid + 1;
         end

      end % RemoveAggsWithOnlyOneBlk

      function [agg, roots] = OrderAggregation(ppd,aggpatt,offset)
      % Generate imperfectly perfect aggregates, as specified by the
      % pattern aggpatt which gets tiled time and time again.  A
      % parameter is included to create an offset in the first aggregate
      % on the left side --- so it doesn't need to line up directly with
      % the mesh.  An offset of n means that we skip the first n nodes in
      % the aggpatt.

        N=prod(ppd);
        Nplane=cumprod(ppd);
        dim=length(ppd);
        agg=ones(N,1);

        % Calculate the stencils
        for I=1:dim,
          AP=aggpatt{I};
          stencil{I}=[];
          for J=1:length(AP),
            stencil{I}=[stencil{I},(J-1)*ones(1,AP(J))];
          end
          ssize(I)=length(stencil{I});
          aps(I)=length(AP);

        end

        % Generate the per-dimension aggregates
        for I=1:dim,
          STC=stencil{I};
          % Generate the header stuff
          header=STC(1+offset(I):length(STC));
          base=length(header);

          % Cleanup for missing aggregates
          if(base > 0), header=header-min(header); max_header=max(header);
          else max_header=0; end

          % Generate the main chunk of stuff
          nresp=ceil( (ppd(I)-base) / ssize(I));
          agg_up=reshape(aps(I)*repmat(0:nresp-1,ssize(I),1),1,ssize(I)*nresp)+ repmat(STC,1,nresp);
          agg_d{I}=[header';1+max_header+agg_up(1:(ppd(I)-base))'];
          nagg(I)=1+max(agg_d{I});

          % Compute root nodes for each dimension
          BIDX=[0,find([0;diff(agg_d{I});length(agg_d{I})]')];
          rrt=BIDX+[floor(diff(BIDX)/2),0];
          rrt_d{I}=rrt(1:length(rrt)-1)';
        end
        Na=prod(nagg);
        Nap=cumprod(nagg);

        % Get the plane indices right
        for I=2:dim,
          agg_d{I} = agg_d{I} * prod(nagg(1:I-1));
          rrt_d{I} = ((rrt_d{I}-1)*Nplane(I-1)+1);
        end

        % Aggregate, calculating root nodes
        roots=zeros(Na,1);
        for I=1:dim,
          agg = agg + reshape(repmat(agg_d{I},N/Nplane(I),Nplane(I)/ppd(I))',N,1);
          roots = roots + reshape(repmat(rrt_d{I},Na/Nap(I),Nap(I)/nagg(I))',Na,1)-1;
        end
        roots=roots+1;

        % This is a HAQ.  I should really fix this.
        roots=sort(roots);
      end



      function [agg,roots]=SplitAggregates(Amat,dim,agg,roots)
      % Make sure no aggregate is bigger than 3^d nodes by
      % splitting offensive aggregates in half.
        aold=agg;
        rold=roots;

      % NTS - Only works in 2D right now
       NumSplit=0;
       Na=max(agg);
       Nanew=Na+1;
       for I=1:Na,
         AIDX=find(agg==I);
         MULT = length(AIDX) / 3;

         if(MULT>3 & MULT<6),
           % Grab the matrix & find minimum degree nodes
           NumSplit=NumSplit+1;
           Alocal=Amat(AIDX,AIDX);
           N=size(Alocal,1);
           DEGREE=(abs(Alocal)>0)*ones(N,1);
           IDX=find(DEGREE <=4);

           if(length(IDX)>4), fprintf('Too many corners\n');keyboard;end
           % Split into four pieces - each corner grabs its
           % neighbors and we hope this is sufficient
           REMAIN=1:N;
           for J=1:length(IDX),
             NEWAGG{J}=find(Alocal(IDX(J),:));
             REMAIN=setdiff(REMAIN,NEWAGG{J});
           end
           if(length(REMAIN)~=0), fprintf('4-way corner re-aggregation did not work.\n');keyboard;end
           roots(I)=AIDX(IDX(1));
           for J=2:length(IDX),
             agg(AIDX(NEWAGG{J}))=Nanew;
             roots(Nanew)=AIDX(IDX(J));
             Nanew=Nanew+1;
           end
         elseif(MULT>6),
             fprintf('That''s one fat aggregate...\n');keyboard;
         end
       end
       fprintf('Split: Splitting %d (%4.1f%%) aggregates\n',NumSplit,100*NumSplit/Na);
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
