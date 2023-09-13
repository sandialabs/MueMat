%% AP_PatternFactory
% concrete implementation of sparsity pattern factory for prolongation and
% restricton operator

classdef AP_PatternFactory < PatternFactory
    properties (Access = private)
        degree_ = [];
        options_ = [];     % options for fattening
        coarseNSDim_ = []; %FIXME should default to something better
        GetCPtsFromFactory_ = []; % used to put injection at Cpts
    end

    methods
        function [this] = AP_PatternFactory(filterType, InitPFact, options)
            % AP_PatternFactory constructor
            %
            %   SYNTAX obj = AP_PatternFactory(filterType,options)
            %
            %     filterType - optional filtering (default = [])
            %     InitPFact  - PFactory object for initial guess 'Ptent'
            %     options    - options stored in options_ (used for
            %                  fattening)
            %
            % AP_Pattern Constructor: default pattern A*P
            % allows for a fine grid filtered matrix A

            % copy constructor
            if nargin == 1 && isa(filterType, class(this)), this.Copy_(filterType,[]); return; end

            if varexist('options'), this.options_ = options; end

            if varexist('InitPFact'), this.InitPFact_ = InitPFact; end;
            % default behaviour: this.InitPFact should be empty
            % then EminPFactory automatically syncronizes InitPFact_ for
            % sparsity pattern and EminP transfer operators

            if varexist('filterType') this.filter_ = filterType;
            else this.filter_ = []; end  % TODO: fixme (call PatternFactory constructor?)

            % set default values
            this.type_ = 'AP';
            this.degree_ = 1;
        end

        function SetNeeds(this, FineLevel, CoarseLevel)
        % Obtain any cross factory specifications

            SetNeeds@PatternFactory(this,FineLevel,CoarseLevel);

            % A-filtering TODO: replace by a Aname_ option
            if this.UseAfiltered_,
               FineLevel.Request('Afiltered');  % released by GetAForPattern
            end
            if ~isempty(this.GetCPtsFromFactory_)
               FineLevel.Request('Aggregates',this.GetCPtsFromFactory_);
            end
            CoarseLevel.Request('P', this.InitPFact_);
        end

        function SetDegree(this,degree)
            % set degree of "pattern polynomial" (A*P)^{degree}
            %
            % SYNTAX SetDegree(degree)
            %
            %  degree   - input integer, degree of pattern polynomial
            this.degree_ = degree;
        end

        function SetFactoryForCPoints(this, Factory)
            % set Factory which will supply Cpts that will be used to
            % put an identity in the pattern correspond to injection
            % at the Cpts
            %
            % SYNTAX SetFactoryForCPoints(Factory)
            %
            %  degree   - input Factory, already invoked Factory which has
            %             left Cpts in the hierarchy.
            this.GetCPtsFromFactory_ = Factory;
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

            if isempty(this.type_) == true
                error('AP_Pattern: type of pattern not set.');
            end
            if ~(strcmp(this.type_,'AP') || strcmp(this.type_,'fatten'))
                error('AP_Pattern: pattern type must be "AP" or "fatten"\n');
            end

            initialP = this.GetInitialP(FineLevel,CoarseLevel);

            Amatrix = this.GetAForPattern(FineLevel,CoarseLevel);

            switch this.type_
                case 'AP'
                      % provide variables from levels
                      Pattern = spones(Amatrix.GetMatrixData()) * spones(initialP.GetMatrixData());
                      for ii=2:this.degree_
                          Pattern = spones(Amatrix.GetMatrixData()) * spones(Pattern);
                      end
%                 case 'fatten'
%                     % provide variables from levels
%                     AggInfo  = FineLevel.Get('Aggregates');
%
%                     % Boonen-Inspired Aggregate Fattening.  Think of it as
%                     % multigrid in the foie gras style.  Do AP first, though.
%                     Pattern = spones(Amatrix.GetMatrixData()) * spones(initialP.GetMatrixData());
%                     for ii=2:this.degree_
%                         Pattern = spones(Amatrix.GetMatrixData()) * spones(Pattern);
%                     end
%                     Pattern = spones(Pattern);
%
%                     Pattern = this.FattenP(AggInfo,Pattern,Amatrix,this.options_);

                otherwise
                    error('PatternFactory.Build():  invalid type of pattern.');
            end %switch
            if ~isempty(this.GetCPtsFromFactory_)
               roots=FineLevel.Get('Aggregates',this.GetCPtsFromFactory_).Roots;
               FineLevel.Release('Aggregates',this.GetCPtsFromFactory_);
               if size(Pattern,2) ~= length(roots),
                  fprintf('Cannot put I@Cpts in pattern when size(pattern,2) neq to the number of roots\n');
                  keyboard;
               end
               Pattern(roots,:) = speye(length(roots),length(roots));
            end

            %this.RestoreOutputLevel();
        end %Build()

        function [Amatrix] = GetAForPattern(this,FineLevel,CoarseLevel)
              %GETAFORPATTERN communication layer function between levels and
              %this factory
              %
              %   SYNTAX   Amatrix = obj.GetAForPattern(FineLevel, CoarseLevel)
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
              if this.UseAfiltered_,
                  fprintf('AP_Pattern: Using filtering for A\n');
                  if ~FineLevel.IsAvailable('Afiltered')
                      fprintf('Filtered matrix was not stored\n'); keyboard;
                  end
                  Amatrix  = FineLevel.Get('Afiltered');
                  FineLevel.Release('Afiltered');
              else
                  Amatrix  = FineLevel.Get('A');
              end
        end %GetAForPattern()

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

              if isempty(this.InitPFact_)
                 error('no initial P factory! error!');
              end

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

    end % end public methods

    methods (Access=private)
        function APpattern=FattenP(this,AggInfo,APpattern,Amat,options)
            %FATTENP Idea: Each row of P must be fat enough s.t.
            % nnz(P(:,i)) >= ncoarse.
            %
            %   SYNTAX   APpattern = obj.FattenP(AggInfo, APpattern, Amat, options);
            %
            %     AggInfo   -
            %     APpattern -
            %     Amat      - an Operator object
            %     options   -

            %[N,NullDim] = size(NS);
            %TODO check w/ CMS that the following two lines are correct
            N = size(Amat,2);
            NullDim = this.coarseNSDim_;
            if ~isfield(options,'NCoarseDofPerNode'), Ncpn = NullDim;
            else Ncpn =options.NCoarseDofPerNode;end

            NAggregates  = length(AggInfo.Roots);
            NCoarse      = NAggregates*Ncpn;

            % This only needs to run if NCoarseDofPerNode<NullDim,
            % otherwise any sane aggregation strategy will get this right.
            if(Ncpn<NullDim),
                nfat=0;
                afat=0;

                for ii=1:NAggregates,
                    nodes = find(AggInfo.NodesInAgg(ii,:));
                    rows  = Node2DOF(nodes,Amat.GetRowMap());
                    nRows = length(rows);
                    CDofs = ((ii-1)*Ncpn+1:ii* Ncpn);

                    Alocal= Amat.matrixData_(rows,rows);
                    Plocal= APpattern(rows,:);

                    % Local agglomerated matrix
                    Anode=sparse(length(nodes),length(nodes));
                    Pnode=sparse(length(nodes),size(Plocal,2)/Ncpn);
                    for jj=1:Ncpn,
                        Anode=Anode + (abs(Alocal(jj:Ncpn:end,jj:Ncpn:end))>0);
                        Pnode=Pnode + (abs(Plocal(jj:Ncpn:end,jj:Ncpn:end))>0);
                    end

                    % Calculate *node* fatness
                    for jj=1:length(nodes),
                        LEN(jj)=length(find(sum(Plocal((jj-1)*Ncpn+[1:Ncpn],:))));
                    end

                    for jj=1:length(nodes),
                        if(LEN(jj)< NullDim && AggInfo.Roots(ii)~=nodes(jj)),
                            % Does my node have a neighbor that isn't skinny?
                            AIDX=find(Anode(jj,:)); IDX1=find(Pnode(jj,:));
                            Nrows=Node2DOF(nodes(jj),Amat.GetRowMap());

                            flag=0;
                            for kk=1:length(AIDX),
                                if(LEN(AIDX(kk))>=NullDim),
                                    IDX2=find(Pnode(AIDX(kk),:));
                                    NIDX=setdiff(IDX2,IDX1);
                                    APpattern(Nrows,(NIDX(1)-1)*Ncpn+[1:Ncpn])=1;
                                    nfat=nfat+1;
                                    flag=1;
                                    break;
                                end
                            end
                            if(flag==0),
                                % Since all of my neighbors are skinny, find
                                % someone in my aggregate that isn't.
                                for kk=1:length(rows),
                                    if(LEN(kk)>=NullDim),
                                        IDX2=find(Pnode(kk,:));
                                        NIDX=setdiff(IDX2,IDX1);
                                        APpattern(Nrows,(NIDX(1)-1)*Ncpn+[1:Ncpn])=1;
                                        afat=afat+1;
                                        flag=1;
                                        break;
                                    end
                                end
                            end
                            if(flag==0),
                                fprintf('My aggregate is very, very lonely\n');
                                keyboard;
                            end
                        end
                    end
                end
                fprintf('Fattening: %d(%3.1f%%) w/ neighbors %d(%3.1f%%) w/ aggs out of %d\n',nfat*Ncpn,100*Ncpn*nfat/N,afat*Ncpn,100*afat*Ncpn/N,N);

                % Sanity check the fattening to make sure it worked
                Pnodeg=sparse(N/Ncpn,size(APpattern,2)/Ncpn);
                for jj=1:Ncpn, Pnodeg=Pnodeg + (abs(APpattern(jj:Ncpn:end,jj:Ncpn:end))>0);end
                LENg=full(sum(Pnodeg,2));
                isRoot = zeros(size(Pnodeg,1),1);
                for jj=1:length(AggInfo.Roots)
                    RIDX = AggInfo.Roots(jj);
                    isRoot(RIDX) = 1;
                end
                IDX=find(LENg < NullDim & ~isRoot);
                if(length(IDX)>0), fprintf('ERROR: Fattening missed something\n');keyboard;end
            end
        end %function FattenP()
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
