%% PgPFactory Build a prolongator and restrictor via the PG-AMG algorithm
% builds a prolongator using the PG-AMG smoothing algorithm
% implements the basic variant given in
% "Sala, Tuminaro; A new Petrov-Galerkin smoothed aggregation preconditioner for
% nonsymmetric linear systems; SIAM J. Sci. Comput., 31, 2008"
%% method
% the following steps occur
% # produce graph from discretization matrix with amalgamation and dropping
% # aggregate vertices of graph
% # make tentative prolongator from the aggregates and the fine level
% nullspace
% # calculate local damping factors omega_pi
% # produce final prolongator
% #   P_final = (I - omega_pi Dinv A) P_tent
% # do the same for the restrictor (in an analog way)
%%
classdef PgPFactory < PFactory

    properties (Access = private)
        InitPFact_
        diagonalView_ = 'current' % point view necessary?
%         forceSmoothing_  = false; % if true, make sure that Ptent is always smoothed (no wiggles!)
        AForSmoothingName_;
    end

    methods   % public functions
        %% constructor
        function [this] = PgPFactory(InitPFact,diagonalView)
            % copy constructor
            if nargin == 1 && isa(InitPFact, class(this)), this.Copy_(InitPFact,[]); return; end

            % constructor sets options
            if varexist('InitPFact'), this.InitPFact_ = InitPFact;
            else this.InitPFact_ = 'default'; end; %TentativePFactory(); end;

            % set diagonal view ?
            if varexist('diagonalView'), this.diagonalView_ = diagonalView; end

            if this.GetOutputLevel() > 5
                fprintf('PgPFactory constructor\n');
            end

            this.AForSmoothingName_ = 'A';
        end

        %% set functions
        function SetDiagonalView(this, diagonalView)
            this.diagonalView_ = diagonalView;
        end

        function SetAForSmoothing(this, name)
            % indicates use of filtered version of A (small entries dropped) within prolognator smoothing
            this.UseAfilteredName_ = name;
        end

        function [ToF] = SupportsRestrictionMode(this)
            ToF = true;
        end

        %% "management" of needs
        function SetNeeds(this, FineLevel, CoarseLevel)
            % check if prolongator is already built
            if this.CheckForReUsableP(FineLevel,CoarseLevel,this) == true
                if this.GetOutputLevel()>5, fprintf('PgPFactory: reUse P\n'); end;
                return;
            end

            TwoLevel = TwoLevels(FineLevel,CoarseLevel);
            TwoLevel.Request('P',this.InitPFact_,CoarseLevel);
        end

        %% build prolongator
        function flag = Build(this,FineLevel,CoarseLevel)
            % construct PG-AMG prolongator
            flag = true;

            %% 0) check if prolongator is already built -> reuse it if possible
            if this.CheckForReUsableP(FineLevel,CoarseLevel,this) == true
                %if CoarseLevel.IsAvailable('P',this.InitPFact_),
                CoarseLevel.Release('P', this.InitPFact_);
                %end
                return;
            end

            TwoLevel = TwoLevels(FineLevel,CoarseLevel);

            %% 1) check for system matrix and prepare fine level nullspace
            Amat = this.GetA(FineLevel);  % get system matrix from level storage (depending on prolongation_mode_)

            %% 2) build initial prolongation operator (based on aggregates)
            P = TwoLevel.Get('P',this.InitPFact_,CoarseLevel);

            %% 3) get diagonal of Amat
            Amat.SwitchToPointView();
            BlkDiag = Amat.GetDiagonal([], this.diagonalView_);   % diagonalView should be "point"-based?
            Amat.SwitchToDefaultView();
            if isempty(BlkDiag.GetApplyInverse())
                BlkDiag.FactorBlkDiag();
            end

            %% 4) calculcate D^{-1} * A * P0
            AP0 = Amat * P;
            DinvAP0 = BlkDiag \ AP0;


            %% 5) calculate diag matrizes with local damping factors
            PP = P.GetMatrixData();

            AA = Amat.GetMatrixData();
            XX = DinvAP0.GetMatrixData();


            %% 6) calculate column based omegas
            nf = size(PP,1);   % number of rows (point view!)
            nc = size(PP,2);   % number of aggregates
            omegajp = zeros(nc,1);    % omegas are supposed to be scalar!
            omegajp = full(diag(PP'*AA'*AA*XX) ./ diag(XX'*AA'*AA*XX));

            % check if there are negative omegas
            if ~isempty(find(omegajp < 0.0))
                if this.prolongation_mode_ == true, fprintf('PG-AMG (prolongation): %i negative colbased omegas found (in omegajp)\n',length(find(omegajp < 0.0)));
                else  fprintf('PG-AMG (restriction ): %i negative colbased omegas found (in omegajp)\n',length(find(omegajp < 0.0))); end;
                % TODO: ML style would be to set these omegas to zero ->
                % what for indefinite problems??
            end

            % debug:
%             if this.prolongation_mode_==true
%                 CoarseLevel.Request('omegajp',this);
%                 CoarseLevel.Set('omegajp',omegajp,this);
%             end;
            % end debug

            if this.GetOutputLevel() > 5
                if this.prolongation_mode_ == true
                    fprintf('PG-AMG (prolongation): minOmegaP: %f    maxOmegaP: %f\n', min(abs(omegajp)), max(abs(omegajp)));
                else
                    fprintf('PG-AMG (restriction ): minOmegaR: %f    maxOmegaR: %f\n', min(abs(omegajp)), max(abs(omegajp)));
                end
            end


            %% 7) convert omegas from column-based to row-based
            % omega_int * P0 = P0 * omegajp
            % -> omega_int (make use of the fact, that nnz(P0) in the i-th
            % row of P0 is exactly 1 for all rows i

            omegapfine = -666*ones(nf,1);
            AtLeastOneDefined = 0;
            DinvAP0data = DinvAP0.GetMatrixData();
            for i=1:nf  % loop over all rows (=DOFs) -> use maps??
                bindx = find(DinvAP0data(i,:));

                for j=1:length(bindx);
                    tmpomega = omegajp(bindx(j));   % get corresponding colbased omega
                    if(omegapfine(i) == -666) omegapfine(i) = tmpomega;
                    elseif (tmpomega < omegapfine(i)) omegapfine(i) = tmpomega; end
                end
                if (omegapfine(i) < 0.0) omegapfine(i) = 0.0; end % make sure all omegas > 0! TODO-> indefinite problems??
            end

            omegapfine = spdiags(omegapfine,0,nf,nf);
            omegapfine = Operator(omegapfine,P.GetRowMap(),P.GetColMap(),@MatlabApply);

            % debug:
%             if this.prolongation_mode_==true
%                 CoarseLevel.Request('omegajfine',this);
%                 CoarseLevel.Set('omegajfine',omegapfine,this);
%             end;
            % end debug

            %           % omega_int * P0 = P0 * omegajp
            %           % -> omega_int (make use of the fact, that nnz(P0) in the i-th
            %           % row of P0 is exactly 1 for all rows i
            %
            %           omegajp = spdiags(omegajp,0,nc,nc);   % diagonal matrix (MATLAB format)
            %           omegajp = Operator(omegajp,Map(nc,1),Map(nc,1),@MatlabApply);
            %
            %
            %           Prhs = P * omegajp;
            %
            %           omegapint = zeros(nf,1);
            %
            %
            %           for i=1:nf % better: use maps
            %               % extract row information of P
            %               rowp = PP(i,:);
            %               Pind = find(rowp);    % find nnz indices within current row
            %               if length(Pind)==0
            %                   error('empty row in P? cannot be!'); % see Axels email
            %                   %2/5/111
            %               elseif length(Pind)>1
            %                   error('more than one entry per row in Ptent? should not be! block matrices not supported!');
            %               else
            %                   PPrhs = Prhs.GetMatrixData();
            %                   Prhsind = find(PPrhs(i,:));
            %                   if(length(Prhsind)~=1)
            %                       error('scaled Ptent has more than one entry per row? hows that possible?');
            %                   end
            %                   omegapint(i) = 1/rowp(Pind) .* PPrhs(i,Prhsind);
            %               end
            %
            %           end
            %
            %           omegapfine = zeros(nf,1);
            %
            %           for i=1:nf
            %               Pind = find(AA(i,:));
            %               Arow = omegapint(Pind);
            %               Arow = min(Arow);
            %               omegapfine(i) = max(0,Arow);
            %
            %           end
            %
            %           clear omegapint;
            %           clear omegajp;
            %
            %           omegapfine = spdiags(omegapfine,0,nf,nf);
            %           omegapfine = Operator(omegapfine,P.GetRowMap(),P.GetColMap(),@MatlabApply);

            %% 8) smooth prolongator/restrictor
            DinvAP0 = DinvAP0.GetMatrixData();

            P = P - omegapfine * DinvAP0;%abs(DinvAP0.GetMatrixData());

            %% 9) some more checks on quality of smoothed prolongator
            % some "mean" checks (for detecting strange wiggling behaviour)
            wigglecols = [];
%             PP = P.GetMatrixData();
%             for j=1:size(PP,2)
%                 if abs(min(PP(:,j))+max(PP(:,j))) < 0.2
%                     wigglecols = [wigglecols;j];
%                 end
%             end

            if ~isempty(wigglecols) %(full(max(abs(mean(P.GetMatrixData())))) < 1e-4) || ...
               %(~isempty(find(abs(mean(P.GetMatrixData())) < 1e-5)) || ..

               warning('prolongation operator basis functions have strange wiggles and/or are nearly zero\n');

                %                wigglecols = find(abs(mean(P.GetMatrixData())) < 1e-5);

               fprintf('repair %i/%i wiggly basis functions (columns) in P \n', length(wigglecols),size(P.GetMatrixData(),2));
               PP = P.GetMatrixData();

               if CoarseLevel.IsAvailable('P',this.InitPFact_)
                    PPqr = CoarseLevel.Get('P',this.InitPFact_).GetMatrixData();
                    fprintf('use P (generated by InitPFact_) for repairing wiggly cols\n');
               else
                    error('oops, what to use for PPqr?');
               end

               % repair zero rows
               PP(:,wigglecols) = PPqr(:,wigglecols);

               P = Operator(PP,P.GetRowMap(),P.GetColMap(),@MatlabApply);
               wigglecols = find(abs(mean(P.GetMatrixData())) < 1e-5);
               if ~isempty(wigglecols)
                    warning('repair of wiggly columns in P failed?');
               end
            end

            %% 10) catch spurious zero rows in P
            % (robustness!)
            % fill these rows and columns with the values of Pinitial
            % This is the same as setting the omega values to zero for these
            % rows and columns
            zerorows = find(max(abs(P.GetMatrixData()'))==0);
            if ~isempty (zerorows)
                fprintf('repair %i zero rows in P \n', length(zerorows));
                PP = P.GetMatrixData();

                if CoarseLevel.IsAvailable('P',this.InitPFact_)
                    PPqr = CoarseLevel.Get('P',this.InitPFact_).GetMatrixData();
                    fprintf('use P (generated by InitPFact_) for repairing zero rows/cols\n');
                else
                    error('oops, what to use for PPqr?');
                end

                % repair zero rows
                PP(zerorows,:) = PPqr(zerorows,:);

                P = Operator(PP,P.GetRowMap(),P.GetColMap(),@MatlabApply);
                zerorows = find(max(abs(PP'))==0);
                if ~isempty(zerorows)
                    warning('repair of zero rows in P failed?'); % see Axels email
                    %2/5/11
                end
            end



            %% 11) store PG-AMG transfer operator in CoarseLevel level class
            this.SetTransferOperator(CoarseLevel,P);

            % this is done in GetInitialP
            %% 12) release requested variables
            CoarseLevel.Release('P',this.InitPFact_);

        end

    end   %  end public methods

    methods (Access = protected)
        %%
        % TODO: this functions should be moved to PFactory
        function [A] = GetA(this, FineLevel);
            %GETA
            %
            %   SYNTAX   [A] = GetA(FineLevel);
            %
            %     FineLevel - Level object for fine level
            %     A         - level system matrix for prolongation smoothing
            %
            % default implementation for PFactory derived prolongation operators
            % with support of "restriction" mode
            % The prolongation_mode_ flag controls what is used as system matrix
            % for smoothing the prolongation operator: A in "prolongation" mode
            % and the transposed of A in the "restriction" mode
            if this.prolongation_mode_ == true
                % PgPFactory is in prolongation mode
                % use system matrix A
                A = FineLevel.Get(this.AForSmoothingName_);
            else
                % PgPFactory is in restriction mode
                % use the transposed of A (downwinding)
                A = FineLevel.Get(this.AForSmoothingName_)';
            end
            FineLevel.Release(this.AForSmoothingName_);
        end

        function SetTransferOperator(this,CoarseLevel,TransferOperator)
            %SETTRANSFEROPERATOR
            %
            %   SYNTAX   SetTransferOperator(CoarseLevel, TransferOperator;
            %
            %     CoarseLevel      - Level object for coarse level
            %     TransferOperator - prolongator (or transposed of restrictor in
            %                        "restriction" mode)
            %
            % stores the result of transfer operator smoothing in level data
            % structure as prolongator (if prolongation_mode_==true) or as
            % restrictor (if prolongation_mode_==false).
            % This is the default implementation for prolongation operator with
            % support of "restriction mode".
            if this.prolongation_mode_ == true
                % PgPFactory is in prolongation mode
                % set prolongation operator
                CoarseLevel.Set('P', TransferOperator, this);
            else
                % PgPFactory is in restriction mode
                % set restriction operator
                CoarseLevel.Set('R', TransferOperator', this);
            end
        end
        %%

        % for copy constructor
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
