classdef ConstraintFactory < TwoLevelFactoryBase

    properties(Access = private)
        PatternFact_                 = [];
        options_                     = [];
        ConstraintWgt_               = 1;
        FCSplitting_                 = [];
    end

    methods
        %         function [this] = ConstraintFactory(PatternFact, InitPFact, options)
        function [this] = ConstraintFactory(PatternFact, FCSplitting, options)
            % Copy constructor
            if nargin == 1 && isa(PatternFact, class(this)), this.Copy_(PatternFact,[]); return; end

            if varexist('PatternFact')  this.PatternFact_  = PatternFact; end;
            if varexist('FCSplitting')  this.FCSplitting_  = FCSplitting; end;
            if varexist('options')      this.options_      = options;     end;

            this.ConstraintWgt_ = 1;
        end

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

        function SetNeeds(this, FineLevel, CoarseLevel)
            if isempty(this.PatternFact_),
                error('no PatternFactory in ConstraintFactory defined. Please define a PatternFactory and provide it to the constraint factory. We do not support no pattern factory for ConstraintFactory.');
            end;
            TwoLevel = TwoLevels(FineLevel,CoarseLevel);

            TwoLevel.Request('Ppattern', this.PatternFact_, CoarseLevel);
            CoarseLevel.Request('NullSpace');
            FineLevel.Request('NullSpace');

            % for constraint dropping we need access to the aggregates
            if(isfield(this.options_,'DropConstraintsMethod')),
                if ~FineLevel.IsRequested('FCSplitting', this.FCSplitting_),
                    this.FCSplitting_.SetNeeds(FineLevel);
                end;
                FineLevel.Request('FCSplitting',this.FCSplitting_);
            end
        end

        function SetOptions(this, options)
          this.options_ = options;
        end

        function SetConstraintWeight(this, alpha)
            %SETCONSTRAINTWEIGHT Sets a scalar weighting factor that scales all constraint equations
            % corresponding rows of P that have too few nonzeros to exactly interpolate
            % the fine nullspace.
            %
            %   SYNTAX   obj = obj.SetConstraintWeight(alpha);
            %
            %     alpha -
            %     obj   -
            this.ConstraintWgt_ = alpha;
        end

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

        function SetPatternFactory(this, PatternFact)
            this.PatternFact_ = PatternFact;
        end

        function [PatternFact] = GetPatternFactory(this)
            PatternFact = this.PatternFact_;
        end

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

        function Build(this, FineLevel, CoarseLevel, Specs)
            TwoLevel = TwoLevels(FineLevel,CoarseLevel);
            % build requested pattern, if not available
            Ppattern = TwoLevel.Get('Ppattern', this.PatternFact_, CoarseLevel);
            cnull    = CoarseLevel.Get('NullSpace'); % what should we do if fnull or cnull is not
            fnull    = FineLevel.Get('NullSpace');   % yet defined? We could call a default method
            % which puts translations for fnull. We could
            % also invoke something to inject fnull to
            % define cnull.

            %% drop constraints
            [toWeight] = this.DropConstraints(FineLevel,CoarseLevel);

            %% build constraints
            [B,ConstraintRhs] = ConstraintFactory.BuildConstraints(Ppattern,cnull,fnull,this.ConstraintWgt_,toWeight);
            %ConstraintFactory.BuildConstraints(Ppattern,cnull,residual,0,toWeight);
            %%% TODO residual (->EminPfactory) or fnull (ConstraintFactory)??

            %% Output
            CoarseLevel.Set('B', B, this);
            CoarseLevel.Set('ConstraintRhs', ConstraintRhs, this);

            %% Release
            CoarseLevel.Release('Ppattern',this.PatternFact_);
            CoarseLevel.Release('NullSpace');
            FineLevel.Release('NullSpace');

        end %BuildConstraints()
    end

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

    methods (Access = private)
        function [toWeight] = DropConstraints(this,FineLevel,CoarseLevel)
            TwoLevel = TwoLevels(FineLevel,CoarseLevel);
            % fixme: TwoLevel should not be used for FCSplitting!
            Ppattern = TwoLevel.Get('Ppattern', this.PatternFact_, CoarseLevel); % released in BuildConstraints
            toWeight = zeros(size(Ppattern,1),1);
            if(isfield(this.options_,'DropConstraintsMethod')),
                % Ray-style dropping detection
                if(strcmp(this.options_.DropConstraintsMethod,'nullspace')),
                    fprintf('Constraint Dropping: Nullspace\n');
                    Amat = FineLevel.Get('A');
                    Amatrixdata = Amat.GetMatrixData();
                    NS = FineLevel.Get('NullSpace'); % released in BuildConstraints
                    nn = Amat.GetRowMap().NDOFs();
                    temp = spdiags( 1./diag(Amatrixdata),[0],nn,nn)*Amatrixdata*NS;
                    temp = temp*spdiags( 1./max(NS)',[0],size(NS,2),size(NS,2));
                    filter = abs(temp) > 1e-4;
                    filter = filter*ones(size(NS,2),1);
                    toWeight(find(filter))=1;
                    % make sure that all the root points are still constrained
                    rrows = Node2DOF(TwoLevel.Get('FCSplitting',this.FCSplitting_,FineLevel).cpoints,Amat.GetRowMap);
                    toWeight(rrows) = 0;
                elseif(strcmp(this.options_.DropConstraintsMethod,'skinny')),
                    % Skinny non-root style dropping
                    fprintf('Constraint Dropping: Skinny non-root\n');
                    FC = TwoLevel.Get('FCSplitting',this.FCSplitting_,FineLevel);
                    if ~isfield(this.options_,'NCoarseDofPerNode') % released in BuildConstraints
                        NS = FineLevel.Get('NullSpace');
                        Ncpn = size(NS,2);
                    else Ncpn = this.options_.NCoarseDofPerNode; end;
                    s=sum(abs(Ppattern)>0,2);
                    isRoot = zeros(size(Ppattern,1),1);
                    for jj=1:length(FC.cpoints)
                        RIDX = (FC.cpoints(jj)-1)*Ncpn + [1:Ncpn];
                        isRoot(RIDX) = 1;
                    end
                    toWeight=(s<Ncpn & ~isRoot);
                elseif(strcmp(this.options_.DropConstraintsMethod,'fpoints'))
                    fprintf('Constraint Dropping: Fpoints\n');
                    Amat = FineLevel.Get('A');
                    toWeight(1:end) = 1;
                    rrows = Node2DOF(TwoLevel.Get('FCSplitting',this.FCSplitting_,FineLevel).cpoints,Amat.GetRowMap());
                    toWeight(rrows) = 0;
                end

                % release variables if "ordered" in SetNeeds
                FineLevel.Release('FCSplitting',this.FCSplitting_);    % release aggregates variable (only needed in this function)
            end

        end % DropConstraints
    end % end methods(private)

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

    methods (Static = true)
        function [B, ConstraintRhs] = BuildConstraints(Ppattern, cnull, rhs, wgtFactor, toWeight)
            %BUILDCONSTRAINTS Build the constrains matrix from the pattern of P and the coarse nullspace vectors.
            %
            %   SYNTAX   [B, ConstraintRhs] = obj.BuildConstraints(Ppattern, cnull, rhs, wgtFactor, toWeight);
            %
            %     Ppattern      - pattern of the prolongator
            %     cnull         - coarse nullspace vectors
            %     rhs           - right hand side
            %     wgtFactor     - scalar that multiplies too-skinny constraint equations
            %     toWeight      - A vector w/ zeros for equations that don't get weighted and ones for ones that do.
            %     B             - matrix of constraints
            %     constraintRhs - corresponding right hand side, reshaped and perhaps with element eliminated
            % (TODO: change the name of constraintRhs?)
            %
            % builds matrix B using cnull and Ppattern
            % and then removes rows, that are not "necessary"

            nulldim = size(cnull,2);
            finenulldim = size(rhs,2);

            % options: debug
            opt_debug = 0;

            % options: 2 methods to remove rows from B
            opt_local_svd  = 1; % default = 1
            opt_global_svd = 0; % default = 0

            if (opt_debug)
                disp('Ppattern ='); disp(''); disp(full(Ppattern));
                nulldim
                % disp('nullspace ='); disp(''); disp(full(nullspace));
                disp('cnull ='); disp(''); disp(full(cnull));
            end

            % Build B (using Ppattern and cnull)
            %
            % For each k in nullspace(:,k), fill in the values of
            % P with nullspace values at the same time shift
            % compress the nonzero matrix columns in each row and
            % shift the matrix columns to the right of the previous
            % rows matrix columns. This shifting is to mimic the
            % fact that each nonzero in the sparsity pattern is
            % actually its own unique unknown. The resulting
            % matrix will be m times as tall as P where
            % m is the number of nullspace vectors. The width of
            % the matrix B is nnz(Ppattern).

            % Get Ppattern in IJV format and sorted by row indices.
            [I,J] = find(Ppattern);
            [I,perm] = sort(I); J = J(perm);

            nnzP  = length(I);        % = nnz(Ppattern)
            nrowP = size(Ppattern,1); % = max(I)

            % Build B : B=(II,JJ,VV) with
            %  II = permutation([I 2I 3I ... (nulldim-1)*I])
            %  JJ = [J  J  J ... J]
            %  VV = nullspace values (if B(i,j)~=0, B(i,j)=cnull(j))
            II = zeros(nnzP*nulldim,1);
            JJ = zeros(nnzP*nulldim,1);
            for i=1:nulldim
                II((i-1)*nnzP+1:i*nnzP) = (I-1)*nulldim + i; % permute B so that rows involving same dofs are adjacent
                JJ((i-1)*nnzP+1:i*nnzP) = (1:nnzP); % so nnzP == size(B,2)
            end

            VV = reshape(double(full(cnull(J,:))), [],1);

            B = sparse(II,JJ,VV); clear I J II JJ VV;

            ConstraintRhs = reshape(rhs',[],1);

            if (opt_debug)
                disp('B (initial) ='); disp(''); disp(full(B));
            end

            % Remove rows from B:
            % Remove constraints which are identical to each other.
            % This can happen, for example, if we have two identical
            % nullspace vectors. It can also happen in more subtle
            % situations. Several methods for removing rows are implemented.

            % Method 1 : Local SVD (default)
            % Linear dependencies between rows in B only appears between rows
            % involving the same unknowns.
            %
            % For each row of Ppattern, the number of corresponding row in B is
            % nulldim. Between each row of Ppattern, there is a column-shift in
            % B. So, we only have to test linear dependencies of rows of B that
            % correspond to the  same line of Ppattern. Thanks to the previous
            % permutation, such rows are adjacent.
            % To resume, we extract and test the rank of subsystems:
            % ( aaaaaaa )  ( bbbbbbb )   ( ccccccc )
            % ( aaaaaaa )  ( bbbbbbb )   ( ccccccc )

            tol = 1e-8;
            fprintf('BuildConstraints: dropping rows of B w/ singular value < %d\n',tol);
            if (opt_local_svd)
                %
                icol=1;  % current col in B (shift)
                irow=1;  % current row in newB
                newnz=1; % nnz of newB (nnz-1)

                if size(Ppattern,2) == 1,  %handles an extreme case
                    temp = Ppattern; temp(:,2) = 0; nnzs = sum(temp');
                else
                    nnzs = sum(Ppattern'); % if we use nnz(Ppattern(i,:)) inside the
                    % loop, it is too slow so we use
                    % sum() to store nnz for each row.
                end

                % preallocate memory for newB
                I = ones(nnzP*nulldim,1);
                J = ones(nnzP*nulldim,1);
                V = zeros(nnzP*nulldim,1);
                % preallocate memory for rhs
                newConstraintRhs = zeros(length(ConstraintRhs),1);

                % We have to apply local SVD to the submatrix
                % M = B(rowstart:rowend,icol:icol+s-1);
                % but as this matlab operation is very slow, we use IJV/CSR format ...
                %
                % B in IJV format, sorted by rows
                [iB,jB,vB] = find(B);
                [iB,perm] = sort(iB); jB = jB(perm); vB = vB(perm);
                %
                % B in CSR format (row_ptr, jB, vB)
                % row_ptr vector stores the locations in the val vector 'vB' that start a row
                offset = zeros(size(B,1),1); % offset(i) = number of nz in row 'i'
                for i=1:length(iB), offset(iB(i)) = offset(iB(i))+1; end
                row_ptr = cat(1, 1, cumsum(offset)+1);

                %wgtVector = 100*ones(length(ConstraintRhs),1); %FIXME for visualization
                wgtVector = ones(length(ConstraintRhs),1);
                numUnderWeighted=0;
                %global gwgt
                %gwgt = 100*ones(length(ConstraintRhs),1);   % for later plotting purposes
                %gstep = 1;

                % Main loop
                for ii=1:size(Ppattern,1)
                    s=nnzs(ii);           % number of unknowns for the current row in
                    % B ( == nnz(Ppattern(ii,:)) ).

                    if s > 0,
                        % Extract a submatrix M of B
                        rowstart = nulldim*(ii-1)+1;
                        rowend = nulldim*ii;
                        % M = B(rowstart:rowend,icol:icol+s-1);         % -> slow
                        ind    = row_ptr(rowstart):row_ptr(rowend+1)-1; % faster
                        iBpart = iB(ind) - (rowstart-1);                %
                        jBpart = jB(ind) - (icol-1);                    %
                        vBpart = vB(ind);                               %
                        M = sparse(iBpart, jBpart, vBpart,rowend-rowstart+1,s);

                        [U,S,dummy]=svd(full(M));    % SVD on the submatrix
                        z = abs(S);
                        % Not sure what is right. Do we want to retain largest singular
                        % values even if all are small (meaning a little submatrix in
                        % BBt is small).  Consider, for example, scaling the null space
                        % by 10^-5. Then, the constraint matrix is small, but we want
                        % to keep it. On the other hand, consider a null space with SOME
                        % small entries within a region. Interpolation in this region
                        % should be unconstrained (because small coarse null space values
                        % interpolate to small values without constraint).

                        k = nnz(z > tol*max(max(abs(z))));       % k : rank of M

                        if size(M,1) == k
                            % M is full rank: the deflation is not needed.
                            %  It seems that using the svd results when it is not mandatory reduces the quality of the constraint matrix.
                            %  For instance, with R2011b and 10 emin steps, ElasticityTest(2, 200, 100, {'emin'}) does not converge
                            %  when newB and newConstraintRhs are computed using U(:,1:k)'.
                            %  The problem seems to (1) depends on the MATLAB version, (2) appears only with a lot of emin steps.
                            newConstraintRhs(irow:irow+k-1) = ConstraintRhs(rowstart:rowend);
                        else
                            % Deduce newB and newConstraintRhs
                            M = U(:,1:k)' * M;
                            newConstraintRhs(irow:irow+k-1) = U(:,1:k)' * ConstraintRhs(rowstart:rowend);
                        end

                        % newB(irow:irow+k-1,icol:icol+s-1) = M;         % -> slow
                        [iM,jM,vM] = find(M);                            % faster
                        iMleng = length(iM);                             %
                        I(newnz:newnz+iMleng-1) = iM(1:iMleng) + irow-1; %
                        J(newnz:newnz+iMleng-1) = jM(1:iMleng) + icol-1; %
                        V(newnz:newnz+iMleng-1) = vM(1:iMleng);          %
                        newnz = newnz + iMleng;                          %

                        if(toWeight(ii)),
                            wgtVector(irow:irow+k-1) = wgtVector(irow:irow+k-1) * wgtFactor;
                            numUnderWeighted = numUnderWeighted + k;
                            %gwgt(gstep  : gstep+k-1)          = wgtVector(irow:irow+k-1);
                        end
                        %if k < size(S,1) gwgt(gstep+k: gstep+size(S,1)-1 ) = -100; end
                        %gstep = gstep + size(S,1);

                        irow = irow + k;
                        icol=icol+s;
                    end
                end %for ii=1:size(Ppattern,1)

                clear M U;

                fprintf('removing %d rows from constraints due to redundancy\n',size(B,1)-(irow-1));

                B = sparse(I,J,V);
                clear I J V
                ConstraintRhs = newConstraintRhs(1:irow-1);

                if wgtFactor ~= 1
                    fprintf('Underweighting %d constraint equations\n',numUnderWeighted);
                    wgtVector = diag( sparse( wgtVector(1:irow-1) ) );
                    B = wgtVector*B;
                end
                clear wgtVector;

                % Now squeeze out any rows in B that may have been zeroed
                % due to weighting. %FIXME: only useful if wgtFactor ~=1 ?
                [iB,jB,vB] = find(B);
                [iB,perm] = sort(iB); jB = jB(perm); vB = vB(perm);
                uiB = unique(iB);  % used to pick out that part of the RHS corresponding to nonzero constraints
                pat = ones( size(B,1),1 );
                newrow = 1;
                last = iB(1);
                for ii = 1:length(iB)
                    if iB(ii) == last
                        iB(ii) = newrow;
                    else
                        last = iB(ii);
                        newrow = newrow+1;
                        iB(ii) = newrow;
                    end
                end

                [m,n] = size(B);
                B = sparse(iB,jB,vB,newrow,n);
                clear iB jB vB
                ConstraintRhs = ConstraintRhs(uiB);
                clear uiB

                clear newConstraintRhs;
                clear icol irow l k;
            end %if (opt_local_svd)

            % Method 2 : Global SVD on B : simpler than 'local SVD' method but slower
            if (opt_global_svd)
                [U,S,dummy]=svd(full(B));
                clear dummy
                k = nnz(abs(S) > tol);
                clear S V;

                fprintf('removing %d rows from constraints due to redundancy\n',size(B,1)-k);

                B = U(:,1:k)'*B;
                ConstraintRhs = U(:,1:k)'*ConstraintRhs;

                clear k U;
            end

            if (opt_debug)
                disp('B (final) ='); disp(''); disp(full(B));
            end

        end %BuildConstraints()

    end %static methods

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
