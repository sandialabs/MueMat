classdef CoarseNSFactory
% Factory for generating a coarse grid representation of the nullspace

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Public methods                                                            %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  methods

    function [this] = CoarseNSFactory(arg)
      %COARSENSFACTORY
      %
      %   SYNTAX   obj = CoarseNSFactory();
      %
      %     obj -

       % Copy constructor
       if nargin == 1 && isa(arg, class(this)), this.Copy_(arg,[]); return; end
       %
    end

    function SetNeeds(this, FineLevel, CoarseLevel)
    end

    function  [Ptent2,CNull] = Build(this,AggInfo,Amat,FNull, options,OutputLevel,varargin)
      %BUILD Build a coarse grid representation of the near nullspace.
      %
      %   SYNTAX   [Ptent2, CNull] = obj.Build(AggInfo, Amat, FNull, AggInfo, options, OutputLevel, varargin);
      %
      %     AggInfo   - Aggregation object
      %     Amat        - an Operator object
      %     FNull       - fine nullspace
      %     options     - option structure
      %     OutputLevel - verbosity level
      %     varargin    - optional parameter
      %     Ptent2      - optional tentative prolongator
      %     CNull       - coarse grid nullspace

      opt_local_svd  = 1; % default = 1
      opt_global_svd = 0; % default = 0
      if isfield(options,'GlobalSvd')
        opt_global_svd = 1;
        opt_local_svd  = 0;
      end

      opt_local_scaling  = 0; opt_global_scaling = 0;
      if (opt_local_svd),  opt_local_scaling  = 1; end % default = 1 for local svd
      if (opt_global_svd), opt_global_scaling = 1; end % default = 1 for global svd

      % Get the null space dimension and choose the number of coarse dofs
      % per node which should be less than or equal to the null space
      % dimension.
      NullDim = size(FNull,2);
      NFine  = Amat.GetRowMap().NDOFs();

      % MakeTentative can take an extra argument that allow to permute
      % columns of the nullspace vectors for experiments.
      if (length(varargin)==0),
        %perm = 1:NullDim;
      elseif (length(varargin)==1),
        if(~isempty(varargin{1})), perm = varargin{1};end
      elseif (length(varargin)>1),
        if(~isempty(varargin{1})), perm = varargin{1};end
        coords=varargin{2};
      end
      if varexist('perm'), iperm(perm)=1:NullDim; end

      if ~isfield(options,'NCoarseDofPerNode')
        NCoarseDofPerNode = NullDim; % Defaults to the same as the
                                     % fine nullspace
      else
        NCoarseDofPerNode = options.NCoarseDofPerNode;  %Note: cannot exceed size of smallest
                                                        %      aggregate due to QR issues.
      end

      % The scaling of null space vectors affects what gets emphasized in the tentative
      % prolongator when NCoarseDofPerNode < NullDim.  A Weight is introduced to reflect
      % the relative importance of modes. For example, we might weight the function '1'
      % heavily if we have 3 vectors corresponding to '1', 'x', 'y'.  Similarly, we might
      % heavily weight translations if we have 3 vectors corresponding to x translation,
      % y translation, and xy rotation.  Here, the nullspace vectors are first normalized
      % and then weighted according to user-provided weights.
      globalWeighting=1;
      if (opt_global_scaling)
        globalWeighting=diag(1./sqrt(diag(FNull'*FNull))); % global scaling
      end
      if isfield(options,'CoarseNullspaceWeightingMultiplier'), %user weighting
          % global scaling before user weighting
        globalWeighting = globalWeighting * diag(options.CoarseNullspaceWeightingMultiplier);
      end
     % CNmode='inject JHU';
     % CNmode='inject CMS';

      %keyboard
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %
      %  Local SVD
      %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (opt_local_svd)
        % For Elasticiy2D, NCoarseDofPerNode=2 and nulldim=3:
        % nullspace vectors are orthogonal (ss=diag([1,1,1])) so
        % 'global SVD' on FNull do nothing and we just throw away
        % one nullspace vector. Depending on the matlab version,
        % the two first vectors of uu ([uu,ss,vv]=svd())
        % corresponds to the two translations mode (R2009bSP1) or
        % the rotation and one translation (R2010a). Results are
        % then matlab version dependant.

        % Ray's idea:
        % Depending on the aggregate, it might be interesting to
        % keep translations/translation vectors or
        % rotation/translation vectors: importance of modes might
        % depend on where we are on the grid.

        % 'local SVD' is an aggregate-wise version of the 'global
        % SVD' code. Take a look to 'global SVD' for more information.

        NAggregates  = length(AggInfo.Roots);
        NCoarse      = NAggregates*NCoarseDofPerNode;
        Ptent        = spalloc(NFine,NCoarse,NFine*NCoarseDofPerNode);
        CNull        = zeros(NAggregates*NCoarseDofPerNode,NullDim);

        tol = 1e-12;
        fprintf('dropping all CNull entries < %g\n',tol);

        %Set up global injection operator.  This is used in the case that
        %part or all of the fine grid nullspace is injected to the coarse grid.
        injectionR = sparse(NCoarse,NFine);
        Arowmap = Amat.GetRowMap();
        for ii=1:NAggregates
          rootDofs = Node2DOF(AggInfo.Roots(ii),Arowmap);
          CDofs    = (ii-1)*NCoarseDofPerNode+1 : ii*NCoarseDofPerNode;
          %injectionR(CDofs,rootDofs) = eye(NullDim - NCoarseDofPerNode);
          injectionR(CDofs,rootDofs) = eye(NCoarseDofPerNode);
        end
        % In case we want to append part of the global nullspace
        appendedCNull = injectionR * FNull;
        for ii=1:size(appendedCNull,2), appendedCNull(:,ii) = appendedCNull(:,ii) / norm(appendedCNull(:,ii)); end


        % Choose nullspace mode
        if isfield(options,'CNullMode'), CNmode=options.CNullMode;
        else CNmode='inject CMS';end

        for ii=1:NAggregates
          nodes = find(AggInfo.NodesInAgg(ii,:));
          rows  = Node2DOF(nodes,Amat.GetRowMap());
          nRows = length(rows);
          CDofs = ((ii-1)*NCoarseDofPerNode+1:ii* NCoarseDofPerNode);

          localWeighting = globalWeighting;
          if (opt_local_scaling)
            colNorms = sum(FNull(rows,:).^2).^(1/2); localScaling = diag(1./colNorms);
            localWeighting = localScaling * localWeighting; % local scaling before user weighting
          end %if (opt_local_scaling)

          % SVD null space to throw away small components if NCoarseDofPerNode < NullDim.
          if(~strcmp(CNmode,'inject CMS')), % svd not used by option 'inject CMS'
            localFNull = full(FNull(rows,:)*localWeighting);
            [uu,ss,vv] = svd( localFNull, 0);
            %fprintf('agg %d\n',ii); disp(diag(ss)')
            %if varexist('perm'), uu=uu(:,perm); end

            % Do usual QR-smoothed aggregation business to get Ptent and the coarse null space
            % FIXME JJH shouldn't be necessary, as SVD gives us orthogonality
            %[Ptent(rows,CDofs),CNull(CDofs,1:NCoarseDofPerNode)]=qr(uu(1:nRows,1:NCoarseDofPerNode),0);

            % The first NCoarseDofPerNode columns of Ptent are just the corresponding columns of the local SVD,
            % i.e., the most locally significant fine grid nullspace components.
            % The coarse grid nullspace is then simply the identity.
            Ptent(rows,CDofs) = sparse(uu(:,1:NCoarseDofPerNode));
            CNull(CDofs,1:NCoarseDofPerNode) = eye(NCoarseDofPerNode);
            %FIXME experimental
            %[Q,R]=qr(localFNull,0);
            %Ptent(rows,CDofs) = Q(:,1:NCoarseDofPerNode);
            %% dimensions are screwy -- don't get full rotation
            %CNull(CDofs,:) = R(1:NCoarseDofPerNode,:);
          end

          % Now figure out how we want to represent the remaining nullspace vectors on the coarse level.

          % Append extra columns to coarse null space to reflect all columns of uu.
          % choose one of the next three
          %%appendedCNull = drop(Ptent(rows, CDofs)'*uu(1:nRows,NCoarseDofPerNode+1:NullDim),tol,0); %local w/ dropping
          %%appendedCNull = Ptent(rows, CDofs)'*uu(1:nRows,NCoarseDofPerNode+1:NullDim);             %local w/ no dropping
          %%appendedCNull = drop(Ptent(rows, CDofs)'*FNull(rows,NCoarseDofPerNode+1:NullDim),tol,0); %global
          %CNull(CDofs,NCoarseDofPerNode+1:NullDim)= appendedCNull;
          %if ~varexist('flag3a'), flag3a = 1; fprintf('  ==> Applying Ptent to remaining NS components\n'); end

          %Inject remainder of fine grid nullspace
          if(strcmp(CNmode,'inject partialFN')),
            if ~varexist('flag3b'), flag3b = 1; fprintf('  ==> Injecting remaining %d fine grid NS components\n', NullDim-NCoarseDofPerNode); end
            %for jj=1:size(appendedCNull,2), appendedCNull(:,jj) = appendedCNull(:,jj) / norm(appendedCNull(:,jj)); end
            CNull(CDofs,NCoarseDofPerNode+1:NullDim)= appendedCNull(CDofs,NCoarseDofPerNode+1:NullDim);
          end

          %Inject remainder of local SVD
          if(strcmp(CNmode,'inject SVD')),
            if ~varexist('flag3c'), flag3c = 1;
              fprintf('  ==> Injecting remaining %d local SVD components\n',NullDim-NCoarseDofPerNode);
            end
            appendedCNull = injectionR(CDofs,rows)*uu(1:nRows,NCoarseDofPerNode+1:NullDim);
            appendedCNull = drop(appendedCNull,tol,0);
            for jj=1:size(appendedCNull,2)
              if norm(appendedCNull(:,jj)) > 1e-12, appendedCNull(:,jj) = appendedCNull(:,jj) / norm(appendedCNull(:,jj));end
            end
            %appendedCNull = injectionR(CDofs,rows) * Q(:,NCoarseDofPerNode+1:NullDim); %FIXME experimental
            CNull(CDofs,NCoarseDofPerNode+1:NullDim)= appendedCNull;
          end

          % Inject entire fine grid nullspace
          if(strcmp(CNmode,'inject JHU')),
            if ~varexist('flag3d'), flag3d = 1; fprintf('  ==> Injecting entire fine grid NS JHU style\n');end
            CNull(CDofs,:) = appendedCNull(CDofs,:);
          end

          % Injection CMS style
          if(strcmp(CNmode,'inject CMS')),
            if ~varexist('flag3e'), flag3e = 1; fprintf('  ==> Injecting entire fine grid NS CMS style\n'); end
            % Ptent(rows,CDofs) = sparse(uu(:,1:NCoarseDofPerNode));
            rootDofs = Node2DOF(AggInfo.Roots(ii),Amat.GetRowMap());
            CNull(CDofs,:)=FNull(rootDofs,:);
% I think this is what we need so that Ptent is consistent with the injected null space.
%
% That is,    FNull(:,1:NCoarseDofPerNode) = Ptent * CNull(:,1:NCoarseDofPerNode)
%
            Ptent(rows,CDofs) = FNull(rows,1:NCoarseDofPerNode)*inv(CNull(CDofs,1:NCoarseDofPerNode));
          end

          % Average CMS style
          if(strcmp(CNmode,'average CMS')),
            if ~varexist('flag3f'), flag3f = 1; fprintf('  ==> Averaging entire fine grid NS CMS style\n'); end
            Ptent(rows,CDofs) = sparse(uu(:,1:NCoarseDofPerNode));
            for jj=1:NCoarseDofPerNode,
              CNull(CDofs(jj),:)=mean(FNull(rows(jj:NCoarseDofPerNode:end),:));
            end
          end

          % Linear interpolation CMS - HAQ for 2D & pretty aggregates
          if(strcmp(CNmode,'linear CMS')),
            if ~varexist('flag3g'), flag3g = 1; fprintf('  ==> Linearly interpolating entire fine grid NS CMS style\n'); end
            if ~varexist('coords'), fprintf('ERROR: Linear interpolation REQUIRES coordinates\n');keyboard;end
            root=AggInfo.Roots(ii);
            NodeLoc=coords(nodes,:);
            Lf=prod(1 - (coords(nodes,:)- repmat(coords(root,:),length(nodes),1))/2,2);
            Lf=Lf/sum(Lf);
            for jj=1:NCoarseDofPerNode,
              CNull(CDofs(jj),:)=Lf'*FNull(rows(jj:NCoarseDofPerNode:end),:);
            end
            Ptent(rows,CDofs) = ones(length(rows),NCoarseDofPerNode);
          end

          % Reverse SVD
          if(~strcmp(CNmode,'inject CMS')),
            if varexist('perm'), CNull(CDofs,:)=CNull(CDofs,iperm); uu=uu(:,iperm); end
            if ~varexist('fooBar'), fooBar = 1; fprintf('\n****** skipping reverse SVD scaling ******\n\n'); end
          end
          %CNull(CDofs,:) = (CNull(CDofs,:)*ss*vv')/localWeighting;
          %dif = FNull(rows,:) - Ptent(rows,CDofs)*CNull(CDofs,:);
          %fprintf('agg %d: column by column norm(FNull - Ptent*CNull) = ',ii);
          %for jj=1:NullDim, fprintf('%10.2e ',norm(dif(:,jj))); end fprintf('\n');
        end %for ii=1:NAggregates
      end % opt_local_svd


      % Manually regnerate the nullspace from the coordinates - HAQ
      % for 2D/3D elasticity
      if(strcmp(CNmode,'regen CMS')),
        if ~varexist('coords'), fprintf('ERROR: Regeneration REQUIRES coordinates\n');keyboard;end
        for ii=1:NAggregates,
          ccoord(ii,:)=coords(AggInfo.Roots(ii),:);
        end
        CNull=build_elastic_rbm(ccoord);
      end


      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %
      %  Global SVD
      %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (opt_global_svd)
        %
        % SVD null space to throw away small components if NCoarseDofPerNode < NullDim.  No
        % attempt is made to satisfy these small components in the tentative prolongator.
        %
        [uu,ss,vv] = svd( full(FNull*globalWeighting), 0); % only uu(:,1:NCoarseDofPerNode)
                                                     % considered when building Ptent
                                                     %
        %
        if varexist('perm'), uu=uu(:,perm); end

        % Allocate tentative prolongator and coarse null space.
        NAggregates  = length(AggInfo.Roots);
        NCoarse      = NAggregates*NCoarseDofPerNode;
        Ptent        = spalloc(NFine,NCoarse,NFine*NCoarseDofPerNode);
        CNull        = zeros(NAggregates*NCoarseDofPerNode,NullDim);

        % Do usual QR-smoothed aggregation business to get Ptent and the coarse null
        % space using uu(:,1:NCoarseDofPerNode) as the fine null space.
        for i=1:NAggregates
          nodes    = find(AggInfo.NodesInAgg(i,:));
          rows = Node2DOF(nodes,Amat.GetRowMap());
          CDofs   = ((i-1)*NCoarseDofPerNode+1:i*NCoarseDofPerNode);
          [Ptent(rows,CDofs),CNull(CDofs,1:NCoarseDofPerNode)]=qr(uu(rows,1:NCoarseDofPerNode),0);
        end

        % Append extra columns to coarse null space to reflect all columns of uu.  Specifically,
        %
        %          Ptent CNull ~ uu
        %   Ptent' Ptent CNull = Ptent'*uu
        %                CNull = Ptent'*uu                (as Ptent'*Ptent = I)
        %                      = [ CNull(:,1:NullDim)  Ptent'uu(:,NCoarseDofPerNode+1:NullDim)];
        %
        % Note: if we ignoring extra columns in recursive calls, they are perhaps not needed?
        %
        CNull(:,NCoarseDofPerNode+1:NullDim)= Ptent'*uu(:,NCoarseDofPerNode+1:NullDim);
        %
        % CNull is appropriate for uu. That is, Ptent*CNull ~ uu. We apply right svd matrices
        % and reverse the scaling to instead obtain Ptent*CNull ~ FNull.
        %
        if varexist('perm') CNull=CNull(:,iperm); uu=uu(:,iperm); end
        CNull = CNull*ss*vv'*inv(globalWeighting);
      end % opt_global_svd

      if (opt_local_svd + opt_global_svd ~=1)
        disp('EminFactory: No SVD or global + local SVD activated');
        %TODO: takes the NCoarseDofPerNode first vectors of FNull if no SVD.
        keyboard;
      end

      %
      % DEBUG: Check results
      %
      temp = Ptent'*Ptent - speye(NCoarse,NCoarse);
      fprintf('Nnz(Ptent^T*Ptent - eye) = %d\n',nnz( abs(temp) > 1e-12));
      dif = FNull - Ptent*CNull;
      fprintf('column by column norm(FNull - Ptent*CNull) = ');
      for i=1:NullDim, fprintf('%10.2e ',norm(dif(:,i))); end
      fprintf('\n');

      fprintf('condest(Ptent''Ptent) = %10.2e\n', condest(Ptent'*Ptent));
      fprintf('norm(P.cnull-fnull)=%f, nnzP = %d\n', norm(full(Ptent*CNull-NFine)), nnz(Ptent));

      % Create col map
      colMap=Map(NAggregates,NCoarseDofPerNode);


      % Ray's Secret Document, Item #4 - Forcing interpolation at
      % root nodes to make sure energy minimization doesn't mess
      % stuff up too badly.
      if (isfield(options,'PtentRootModifications'))
        if(strcmp(options.PtentRootModifications,'4b'))
          fprintf('Running Ray''s Special Sauce, Option #4b\n');

          P = Operator(Ptent,Amat.GetRowMap(),colMap,@MatlabApply); % to be fixed ....

          [P, CNull] = this.aggregateWiseQR(AggInfo, P, CNull);

          fprintf('  condest(Ptent''Ptent) = %10.2e\n', condest(P.GetMatrixData()'*P.GetMatrixData()));
          fprintf('  norm(P.cnull-fnull)=%f, nnzP = %d\n', norm(full(P*CNull-NFine)), nnz(P.GetMatrixData()));

          %
          Ptent=P.GetMatrixData(); % to be fixed ...
          clear P;

        elseif(strcmp(options.PtentRootModifications,'4c')),
          fprintf('Running Ray''s Special Sauce, Option #4c\n');

          % Find the coarse unknowns corresponding to the root
          % nodes and replace that chunk of Ptent with the identity
          NAggregates  = length(AggInfo.Roots);
          for I=1:NAggregates,
            rnode=AggInfo.Roots(I);
            rows = Node2DOF(rnode,Amat.GetRowMap());
            CDofs   = ((I-1)*NCoarseDofPerNode+1:I*NCoarseDofPerNode);
            if(length(rows)~=length(CDofs)), fprintf('Identity Injection does not make sense here!\n');keyboard;end
            Ptent(rows,CDofs)=speye(NCoarseDofPerNode,NCoarseDofPerNode);
          end

          % Update CNull
          CNull=Ptent'*FNull;
        end %elseif

      end  %if (isfield(options,'PtentRootModifications'))

      % Wrap Ptent for output
      Ptent2=Operator(Ptent,Amat.GetRowMap(),colMap,@MatlabApply);
      fprintf('condest(final Ptent''Ptent) = %10.2e\n', condest(Ptent'*Ptent));

    end %Build()

  end %public methods

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Private methods                                                           %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  methods (Access=private)

    function [newP, CNull] = aggregateWiseQR(this,AggInfo, P, CNull)
      %AGGREGATEWISEQR
      %
      %   SYNTAX   [newP, CNull] = obj.aggregateWiseQR(AggInfo, P, CNull);
      %
      %     AggInfo - aggregate structure
      %     P         - prolongator Operator
      %     CNull     - coarse nullspace
      %     newP      - new prolongator
      Pmatrixdata = P.GetMatrixData();

      NodesInAgg  = AggInfo.NodesInAgg;
      NAggregates = max(AggInfo.AggId);

      NFine     = P.GetRowMap().NDOFs();  %size(P,1);
      NCoarse   = size(Pmatrixdata,2);    %NAggregates*nulldim;
      nulldim   = NCoarse/NAggregates;

      RowMap = P.GetRowMap();
      for i=1:NAggregates
        BlkRows = find(NodesInAgg(i,:));
        rnode=AggInfo.Roots(i);
        rows = Node2DOF(rnode,RowMap); % rows = Node2DOF(BlkRows,RowMap);

        CDofs   = ((i-1)*nulldim+1:i*nulldim);

        localP = Pmatrixdata(rows,CDofs);
        [localQ,localR] = qr(localP,0);
        %if ( size(localR,1) < length(CDofs) )
        if ( rank(full(localR)) < length(CDofs) )
          localR(length(CDofs),length(CDofs)) = 1.;
        end

        % debug: R(CDofs,CDofs) = localR;

        Pmatrixdata(:,CDofs) = Pmatrixdata(:,CDofs) / localR;
        CNull(CDofs,:)       = localR * CNull(CDofs,:);
      end

      % debug: Pmatrixdata = Pmatrixdata/R;
      % debug: CNull = R * CNull;

      newP = Operator(Pmatrixdata, P.GetRowMap(), P.GetColMap(), P.GetApply());
    end %function aggregateWiseQR()

  end %private methods

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

end %class CoarseNSFactory
