classdef TestHybrid2x2PFactory < PFactory
    
    properties (Access = private)
        Factory11_
        Factory22_
    end
    
    methods
        
        function [this] = TestHybrid2x2PFactory(Factory11, Factory22)
            this.Factory11_ = Factory11;
            this.Factory22_ = Factory22;
        end
        
        function [ToF] = SupportsRestrictionMode(this)
            ToF = false;
        end
               
        function [ToF] = ReUseP(this, ToF)
            % not supported!
            ToF = false;
        end
        
        function SetNeeds(this, FineLevel, CoarseLevel)
            % Obtain any cross factory specifications      
            if ~isempty(this.Factory11_),
                if ~CoarseLevel.IsRequested('P', this.Factory11_)
                    this.Factory11_.SetNeeds(FineLevel,CoarseLevel);
                end
                CoarseLevel.Request('P', this.Factory11_);
            end
            if ~isempty(this.Factory22_),
                if ~CoarseLevel.IsRequested('P', this.Factory22_)
                    this.Factory22_.SetNeeds(FineLevel,CoarseLevel);
                end
                CoarseLevel.Request('P', this.Factory22_);
            end
            
            FineLevel.Request('Arr');
            FineLevel.Request('dof2node');            
        end
        
        function flag = Build(this,FineLevel,CoarseLevel)
            
            flag     = true;
                   
            oldA = FineLevel.Get('A');
            FineLevel.Set('A',FineLevel.Get('Arr'));
            flag11 = this.Factory11_.Build(FineLevel,CoarseLevel); %internally needs dof2node
            if flag11 == false, flag = false; return; end;
            P11  = CoarseLevel.Get('P',this.Factory11_).GetMatrixData();
            n11  = size(P11,1);
            n22  = oldA.GetRowMap().NDOFs() - n11;
            BigP = [P11 sparse(n11,n22) ; sparse(n22,size(P11,2)) speye(n22,n22)];
            CoarseLevel.Set('P', Operator(BigP,1,1),this);
            FineLevel.Set('A',oldA);

            if FineLevel.IsAvailable('AuxMatrix')
                % Set AuxMatTransferts
                CoarseLevel.Set('AuxMatP',CoarseLevel.Get('P',this.Factory11_));
                CoarseLevel.Set('AuxMatR',Transpose(CoarseLevel.Get('P',this.Factory11_)));
                
                % Set Coord Transfers
                CoarseLevel.Set('coordP',CoarseLevel.Get('P',this.Factory11_));
                CoarseLevel.Set('coordR', Transpose(CoarseLevel.Get('P',this.Factory11_)));
            end

            % release variables
            if CoarseLevel.IsAvailable('P', this.Factory11_), CoarseLevel.Release('P', this.Factory11_); end;
            if CoarseLevel.IsAvailable('P', this.Factory22_), CoarseLevel.Release('P', this.Factory22_); end;
            FineLevel.Release('Arr');
            FineLevel.Release('dof2node');
        end
    end
end

