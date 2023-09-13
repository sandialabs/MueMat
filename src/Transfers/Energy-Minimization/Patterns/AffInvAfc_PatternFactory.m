classdef AffInvAfc_PatternFactory < PatternFactory
    %AffInvAfc_PATTERN Factory
    % factory to create prolongator sparsity patterns based on optimal
    % prolongation operator

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties                                                                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (Access = private)
        FCSplitting_ = [];
        AffInverseApproxMethod_ = 'neumann';
        AffInverseIterations_ = 25; % number of iterations for approximating Aff^{-1}
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Public methods                                                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        function [this] = AffInvAfc_PatternFactory(filterType, FCSplitting, InitPFact, options)
            %AffInvAfc_PATTERN Constructor
            %
            %   SYNTAX   obj = PatternFactory(filterType, options);
            %
            %     filterType  - filter object for filtering pattern (e.g.
            %                   Thresholding), default = []
            %     InitPFact   - factory for initial P
            %     options     -
            %     obj         -

            % Copy constructor
            if nargin == 1 && isa(filterType, class(this)), this.Copy_(filterType,[]); return; end

            if varexist('FCSplitting'), this.FCSplitting_ = FCSplitting;
            else error('AffInvAfc_PatternFactory: no FCSplitting\n'); end;

            if varexist('InitPFact'), this.InitPFact_ = InitPFact; end;
            % default behaviour: this.InitPFact should be empty
            % then EminPFactory automatically syncronizes InitPFact_ for
            % sparsity pattern and EminP transfer operators

            %
            if varexist('filterType') this.filter_ = filterType;
            else this.filter_ = []; end;

            %       if varexist('options'), this.options_ = options; end
            this.type_ = 'AffInvAfc';
        end %ctor


        function SetNeeds(this, FineLevel, CoarseLevel)
            % Obtain any cross factory specifications
            SetNeeds@PatternFactory(this,FineLevel,CoarseLevel);

            CoarseLevel.Request('NullSpace');
            FineLevel.Request('NullSpace');
            if ~FineLevel.IsRequested('FCSplitting',this.FCSplitting_),
                this.FCSplitting_.SetNeeds(FineLevel);
            end;
            FineLevel.Request('FCSplitting',this.FCSplitting_);
%             if ~FineLevel.IsRequested('P',this.InitPFact_),
%                 this.InitPFact_.SetNeeds(FineLevel,CoarseLevel);
%             end;
            CoarseLevel.Request('P', this.InitPFact_); % request InitPFact_ second time?
        end

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
        function SetAffInverseIterations (this, numIters)
            this.AffInverseIterations_ = numIters;
        end

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
        function SetAffInverseApproxMethod (this, meth)
            this.AffInverseApproxMethod_ = meth;
        end

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

        function [Pattern] = BuildPattern(this, FineLevel, CoarseLevel, Specs)
            %BUILD Build a pattern based upon smoothed aggregation.
            %
            %   SYNTAX   Pattern = obj.Build(FineLevel, CoarseLevel, Specs)
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
            %     Specs            - specs (input)
            %     Pattern          - binary sparsity pattern (output)

            if nargout() ~= 1, error('Build method returns 1 argument'); end

            % check if all needed input is provided
            %if ~FineLevel.IsAvailable('FCSplitting',this.FCSplitting_) error('no FCSplitting? You need FCSplitting as parameter in Build function'); end;
            if ~FineLevel.IsAvailable('NullSpace')   error('no Bzero? You need Bzero as parameter in Build function'); end;
            if ~CoarseLevel.IsAvailable('NullSpace') error('no Bone? You need Bone as parameter in Build function'); end;

            % provide variables
            Amatrix  = FineLevel.Get('A');
            Bzero    = FineLevel.Get('NullSpace'); FineLevel.Release('NullSpace');
            Bone     = CoarseLevel.Get('NullSpace'); CoarseLevel.Release('NullSpace');
            if ~FineLevel.IsAvailable('FCSplitting', this.FCSplitting_),
                this.FCSplitting_.Build(FineLevel);
            end
            AggInfo  = FineLevel.Get('FCSplitting', this.FCSplitting_); FineLevel.Release('FCSplitting',this.FCSplitting_);
            LevelID  = FineLevel.GetLevelId();

            initialP = this.GetInitialP(FineLevel,CoarseLevel);

            if isempty(this.type_) == true
                error('PatternFactory: type of pattern not set.');
            end


            nFine = size(initialP.GetMatrixData(),1);   % use maps for this? %TODO DOF <-> NODE
            nCoarse = size(initialP.GetMatrixData(),2);

            roots = Node2DOF(AggInfo.cpoints,initialP.GetRowMap()); % determine DOFs of root nodes
            fpoints = ones(nFine,1); fpoints(roots) = 0; fpoints = find(fpoints); % DOFs of fine level nodes

            % extract Aff
            AA = Amatrix.GetMatrixData();
            Aff = AA(fpoints,fpoints);

            rightnull = Bzero;    % default: same approx for left and right nullspace (maybe not the best idea for nonsymmetric problems!)
            rightnullout = Bone;
            leftnull = Bzero;
            leftnullout = Bone;

            % TODO: improve nullspaces

            % setup sparse pattern matrix (n x n_c)
            %Pattern = sparse(nFine,length(roots));
            %Pattern(roots,:) = speye(length(roots),length(roots));
            Pattern = sparse(nFine,nCoarse);
            if(length(roots) == nCoarse)
                Pattern(roots,:) = speye(nCoarse,nCoarse);
            else
                Ptent = initialP.GetMatrixData();
                Pattern(roots,:) = spones(Ptent(roots,:));
                clear Ptent;
            end
            invAffAfc = this.AffSolve(AA,Aff,fpoints,roots,Pattern,sprintf('%s',this.AffInverseApproxMethod_),this.AffInverseIterations_,[]); % params (used for symgausseidel, but not needed)

            Pattern(fpoints,:) = invAffAfc(fpoints,:); % return Pattern matrix data

        end %Build()


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
                this.InitPFact_.Build(FineLevel,CoarseLevel);   % bad: this overwrites CoarseLevel.'P'
            end
            initialP = CoarseLevel.Get('P', this.InitPFact_);
            CoarseLevel.Release('P', this.InitPFact_);

            if ~isempty(this.filter_)
                % provide information for filter
                CoarseLevel.Set('PtentForFilter',initialP);
            end
        end %GetInitialP()
    end %public methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Private methods                                                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access=private)

        %% AffSolve
        % approximately solve the equation
        %
        %             Aff Z + Afc = 0
        %
        % for Z, where Z may be a matrix or a vector
        % -> Z = -Aff^{-1} Afc
        %
        %  SYNTAX z = AffSolve(A,Aff,fpoints,roots,z0,strategy,iters,params)
        %
        %  A                - input matrix (matlab format)
        %  Aff              - only fpoint part of A
        %  fpoints          - fine grid points
        %  roots            - root poits (=coarse grid points)
        %  z0               - initial guess
        %  strategy         - can be 'direct', 'neumann', 'iterative', method
        %                     for solving above system of equations
        %  iters            - number of iterations
        %  params           - some parameters
        %
        % some comments:
        % if strategy = iterative -> use iters sym gauss seidel iterations
        % if strategy = neumann -> use truncated neumann series (iters =
        %                          truncation index)
        % note: z(roots,:) = z0(roots,:)!
        function [z] = AffSolve(this, A, Aff, fpoints, roots, z0, strategy, iters, params)

            if(strcmp(strategy,'direct'))
                %% direct solve
                z=z0;
                z(fpoints,:) = -Aff\(A(fpoints,roots)*z0(roots,:));
            elseif (strcmp(strategy,'iterative'))
                %% iterative solution with symmetric gauss seidel
                % initial vector
                z = z0;

                % symmetric gauss seidel (iters iterations)
                temprhs = -A(fpoints,roots)*z(roots,:);
                tempsol = z0(fpoints,:);

                for i=1:iters
                    % AffSolve with Nits too big -> numerical breakdown
                    res = (temprhs-Aff*tempsol); %
                    if this.GetOutputLevel() > 7, fprintf('norm(res) = %e\n', norm(res,'fro')); end;
                    % so we check res for convergence (just for safety)
                    if abs(norm(res,'fro')) < 1e-12
                        break;
                    end

                    tempsol = this.symgausseidel(Aff,tempsol,temprhs, params);
                end;
                z(fpoints,:) = tempsol;
            elseif (strcmp(strategy,'jacobi'))
                %% iterative solution with jacobi iteration
                % z0(roots,:) contains the identity block
                z = z0;
                b = -A(fpoints,roots) * z0(roots,:); % rhs
                AA = A(fpoints,fpoints);
                ddd = spdiags(1./diag(AA),[0],size(AA,1),size(AA,2));

                dddb = ddd*b;

                %omega = 0.3;
                maxnormAff = norm(ddd*AA,inf);
                omega = 1/maxnormAff;
                fprintf('CHECK ME: omega for pattern Jacobi: %f\n',omega);
%                 [V,D] = eig(full(speye(size(AA,1)) - omega*ddd*AA));
%                 omega = (min(abs(diag(D)))+max(abs(diag(D))))/2;
%                 fprintf('optimal omega: %f\n',omega);
%                 omega = eigs(speye(size(AA,1)) - omega*ddd*AA,1)/2;

                tempsol = z0(fpoints,:); % fine level part of solution
                for i=1:iters
                    tempsol = (speye(size(AA,1)) - omega*ddd*AA) * tempsol + omega*ddd*dddb;
                end
                z(fpoints,:) = tempsol;
            elseif (strcmp(strategy,'richardson'))
 %% iterative solution with jacobi iteration
                % z0(roots,:) contains the identity block
                z = z0;
                b = -A(fpoints,roots) * z0(roots,:); % rhs
                AA = A(fpoints,fpoints);

                maxnormAff = norm(AA,inf);

                omega = 1/maxnormAff; % TODO: check me!!! EW zu teuer!
%                 omega = 0.3;
%                 fprintf('CHECK ME: omega for pattern Jacobi: %f\n',omega);
%                 [V,D] = eig(full(speye(size(AA,1)) - omega*ddd*AA));
%                 omega = (min(abs(diag(D)))+max(abs(diag(D))))/2;
%                 fprintf('optimal omega: %f\n',omega);
%                 omega = eigs(speye(size(AA,1)) - omega*ddd*AA,1)/2;

                tempsol = z0(fpoints,:); % fine level part of solution
                for i=1:iters
                    tempsol = (speye(size(AA,1)) - omega * AA) * tempsol + omega*b;
                end
                z(fpoints,:) = tempsol;
            elseif (strcmp(strategy,'neumann'))
                %% approximate inverse of Aff with neumann series

                maxnormAff = norm(Aff,inf);

                gamma = 1/(maxnormAff+0); % TODO: check me!!! EW zu teuer!

                IminusGammaA = speye(size(A)) - gamma * A;
                if this.GetOutputLevel() > 7, fprintf('AffInvAfc_PatternFactory(neumann): norm(I-gamma * A) = %e\n',norm(IminusGammaA,inf)); end;
                if norm(IminusGammaA,inf) > 1.0, warning('AffInvAfc_PatternFactory(neumann): norm(I-gamma * A) = %e > 1\n',norm(IminusGammaA,inf)); end;
                clear IminusGammaA;

                naff = size(Aff,1);
                IminusgammaAff = speye(naff)-gamma*Aff;
                IminusgammaAffpot = speye(naff);

                %Affinv = IminusgammaAff; % error in computations from July
                %1st
                Affinv = speye(naff);

                for i=1:iters
                    IminusgammaAffpot = IminusgammaAffpot * IminusgammaAff;  % matrix-matrix multiplication
                    Affinv = Affinv + IminusgammaAffpot;
                end

                Affinv = gamma * Affinv;

                z = z0;
                z(fpoints,:) = - Affinv * A(fpoints,roots) * z0(roots,:);   % matrix-matrix-matrix?
            else
                error('ERROR: strategy for AffInv not defined! Use direct, iterative or neumann.');
            end
        end % function AffSolve()

        %%
        function [v]= symgausseidel(this,A,v,f, parameters)
            % SYMGS Gauss-Seidel relaxation
            %
            % SYNTAX [v] = gs(f,A,v,its)
            %
            %   f  - right hand side
            %   A  - matrix
            %   v  - initial guess

            L = tril(A,-1);
            U = triu(A,1);
            %fprintf('sizeA in SymGaussSeidel %i, %i\n',size(A,1),size(A,2));
            D = spdiags(diag(A),[0],size(A,1),size(A,1));

            %%% forward sweep
            v = (D + L) \ (f - U * v);
            %%% backward sweep
            v = (D + U) \ (f - L *v);
        end % function symgausseidel
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

end %class PatternFactory
