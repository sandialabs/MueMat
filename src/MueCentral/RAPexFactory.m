classdef RAPexFactory < RAPFactory

    properties (Access = private)
      AggDropTol_ = 0; % for 'AfilteredDecoalesced'
      DropFunc_     % User supplied function which takes the RAP matrix and drops entries in it
      DropData_     % User supplied data passed to DropFunc_()
    end

    % Factory which builds coarse discretizations via R*A*P
    methods
        function [this] = RAPexFactory()
         this.DropData_  = [];
         this.DropFunc_  = [];
        %TODO: Cpy constructor
        end

        function SetNeeds(this, FineLevel, CoarseLevel)
            % Obtain any cross factory specifications
        end
        function SetDropSpecs(this,DropFunc,DropData)
      % Assign predrop function used to decide matrix elements ignored before amalgamation occurs
         this.DropFunc_ = DropFunc;
         if varexist('DropData'), this.DropData_ = DropData; end
      end


        function SetAggDropTol(this, AggDropTol)
            this.AggDropTol_ = AggDropTol; % for 'AfilteredDecoalesced'
        end

        function flag = Build(this,FineLevel, CoarseLevel)
            % build CoarseLevel A via (CoarseLevel R)*(FineLevel A)*(CoarseLevel P)

            % call simple standard RAPFactory.Build
            flag = Build@RAPFactory(this,FineLevel,CoarseLevel);

            % add projection of additional variables
            Pmat = CoarseLevel.Get('P');
            Rmat = CoarseLevel.Get('R');

            % NonGalerkin coarse grids

            if ~isempty(this.DropData_) && ~isempty(this.DropFunc_),
               Amat = CoarseLevel.Get('A');
               ActualA  = Amat.GetMatrixData;
               Afiltered = this.DropFunc_(ActualA,CoarseLevel,this.DropData_);
               if isfield(this.DropData_,'FixRowSums'),
                  nnn     = size(ActualA,1);
                  GoodRowSums = ActualA*ones(nnn,1);
                  BadRowSums  = Afiltered*ones(nnn,1);
                  % do it like this for rounding errors
                  Afiltered = Afiltered - spdiags(BadRowSums,0,nnn,nnn);
                  Afiltered = Afiltered + spdiags(GoodRowSums,0,nnn,nnn);
               end
               NewA=Operator(Afiltered,Amat.GetRowMap,Amat.GetColMap,@MatlabApply);
               CoarseLevel.Set('A', NewA);
            end

            if FineLevel.IsAvailable('xcoords')

                Paux = Pmat;
                Raux = Rmat;
                Type='AverageVector';
                Data = FineLevel.Get('xcoords');
                if CoarseLevel.IsAvailable('coordP') Paux = CoarseLevel.Get('coordP'); end
                if CoarseLevel.IsAvailable('coordR') Raux = CoarseLevel.Get('coordR'); end

                CoarseLevel.Set('xcoords', this.ProjectIt(Type,Raux,Paux,Data));
            end
            if FineLevel.IsAvailable('ycoords')

                Paux = Pmat;
                Raux = Rmat;
                Type='AverageVector';
                Data = FineLevel.Get('ycoords');
                if CoarseLevel.IsAvailable('coordP') Paux = CoarseLevel.Get('coordP'); end
                if CoarseLevel.IsAvailable('coordR') Raux = CoarseLevel.Get('coordR'); end

                CoarseLevel.Set('ycoords', this.ProjectIt(Type,Raux,Paux,Data));
            end
            if FineLevel.IsAvailable('zcoords')

                Paux = Pmat;
                Raux = Rmat;
                Type='AverageVector';
                Data = FineLevel.Get('zcoords');
                if CoarseLevel.IsAvailable('coordP') Paux = CoarseLevel.Get('coordP'); end
                if CoarseLevel.IsAvailable('coordR') Raux = CoarseLevel.Get('coordR'); end

                CoarseLevel.Set('zcoords', this.ProjectIt(Type,Raux,Paux,Data));
            end

            % Project AuxMatrix using 'ProjectIt' when there is no
            % AuxMatrixFunc to define how to do so.
            if FineLevel.IsAvailable('AuxMatrix') && ~FineLevel.IsAvailable('AuxMatrixFunc')
                NRdofs  = CoarseLevel.Get('R').GetColMap().NDOFs();
                Paux = Pmat;
                Raux = Rmat;
                Type='StandardMatrix';
                Data = FineLevel.Get('AuxMatrix');
                if CoarseLevel.IsAvailable('AuxMatP') Paux = CoarseLevel.Get('AuxMatP'); end
                if CoarseLevel.IsAvailable('AuxMatR')
                    Raux = CoarseLevel.Get('AuxMatR');
                    NRdofs  = CoarseLevel.Get('AuxMatR').GetColMap().NDOFs();
                end
                if FineLevel.Get('AuxMatrix').GetRowMap().NDOFs() ~= NRdofs,
                    Type='AverageMatrix';
                end

                CoarseLevel.Set('AuxMatrix', this.ProjectIt(Type,Raux,Paux,Data));
            end

            if FineLevel.IsAvailable('AuxMatrixFunc')
                AMF=FineLevel.Get('AuxMatrixFunc');
                fprintf('RAP: Regenerating AuxMatrix using the Aux Matrix Func (AMF)\n');
                ExtractFact = CoalesceDropFactory();
                Anode=ExtractFact.Build_(CoarseLevel.Get('A'),CoarseLevel);
                CoarseLevel.Set('AuxMatrixFunc', AMF);

                xcoords = []; if CoarseLevel.IsAvailable('xcoords'), xcoords = CoarseLevel.Get('xcoords'); end
                ycoords = []; if CoarseLevel.IsAvailable('ycoords'), ycoords = CoarseLevel.Get('ycoords'); end
                zcoords = []; if CoarseLevel.IsAvailable('zcoords'), zcoords = CoarseLevel.Get('zcoords'); end

                CoarseLevel.Set('AuxMatrix', AMF(Anode,[xcoords,ycoords,zcoords]));
            end

            if FineLevel.IsAvailable('AfilteredDecoalesced') && CoarseLevel.IsAvailable('AuxMatrix')
                CAmat = CoarseLevel.Get('A');
                AuxMat = CoarseLevel.Get('AuxMatrix');

                FilterFact = CoalesceDropFactory();
                FilterFact.SetPreDropSpecifications(@SAdrop, this.AggDropTol_);
                filteredAuxMat = FilterFact.Build_(AuxMat);

                filteredA = BuildDecoalescedFilteredMatrix(CAmat, filteredAuxMat);
                CoarseLevel.Set('AfilteredDecoalesced', filteredA);
            end

        end
    end
    methods (Access = private)

        function Result = ProjectIt(this, Type, Rmat, Pmat, Data)

            if strcmp(Type,'StandardMatrix'),
                Result = Rmat*(Data*Pmat);
            elseif strcmp(Type,'AverageMatrix'),
                ExtractFact = CoalesceDropFactory();
                ExtractFact.SetPostDropSpecifications([],[]);
                ExtractFact.SetPreDropSpecifications([],[]);
                ExtractFact.SetBinary();
                R = ExtractFact.Build_(Rmat,[]);
                P = ExtractFact.Build_(Pmat,[]);
                % should probably do some kind of scaling
                Result = R*(Data*P);
            elseif strcmp(Type,'AverageVector'),
                NRdofs = Rmat.GetColMap().NDOFs();
                if length(Data) == NRdofs,
                    oo = ones(NRdofs,1);
                    temp = Rmat.Copy();   temp.SetMatrixData(Rmat.GetMatrixData() ~= 0);
                    %TODO  replace the following with temp * oo
                    rowsums = temp.Apply(oo, []);
                    %TODO  replace the following with temp ./ rowsums
                    Result =  temp.Apply(Data, []) ./ rowsums;
                else
                    ExtractFact = CoalesceDropFactory();
                    ExtractFact.SetPostDropSpecifications([],[]);
                    ExtractFact.SetPreDropSpecifications([],[]);
                    ExtractFact.SetBinary();
                    R = ExtractFact.Build_(Rmat, []);
                    NRdofs = R.GetColMap().NDOFs();
                    oo = ones(NRdofs,1);
                    applyFcn = R.GetApply();
                    rowsums = applyFcn(R, oo, []);
                    Result =  applyFcn(R, Data, []) ./ rowsums;
                end
            end
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
