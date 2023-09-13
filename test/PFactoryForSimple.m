classdef PFactoryForSimple < PFactory
    
    properties (Access = private)
        Factory11_
        Factory22_
        readme_ = []
    end
    
    methods
        
        function [this] = PFactoryForSimple(Factory11, Factory22)
            this.Factory11_ = Factory11;
            this.Factory22_ = Factory22;
        end
        function [this] = SetReadFile(this,str)
            this.readme_ = str;
        end;
               
        function [ToF] = SupportsRestrictionMode(this)
            ToF = false;
        end
        
        function SetNeeds(this, FineLevel, CoarseLevel)
            % check if prolongator is already built
            if this.CheckForReUsableP(FineLevel,CoarseLevel,this) == true
                if this.GetOutputLevel()>5, fprintf('PFactoryForSimple: reUse P\n'); end;
                return;
            end
            
            % Obtain any cross factory specifications
            if ~isempty(this.Factory11_),
                CoarseLevel.Request('P', this.Factory11_);
                this.Factory11_.SetNeeds(FineLevel,CoarseLevel);
            end
            if ~isempty(this.Factory22_),
                CoarseLevel.Request('P', this.Factory22_);
                this.Factory22_.SetNeeds(FineLevel,CoarseLevel);
            end
        end
        
        function flag = Build(this,FineLevel,CoarseLevel)
            
            flag     = true;
            
            % check if prolongator is already built
            if this.CheckForReUsableP(FineLevel,CoarseLevel,this) == true
                return;
            end
            
            NS    = FineLevel.Get('NullSpace');
            
            
            oldA  = FineLevel.Get('A');
            oldNS = FineLevel.Get('NullSpace');
            FineN11  = FineLevel.Get('N11');
            FineLevel.Set('A', FineLevel.Get('A11'));
            FineLevel.Set('NullSpace', NS.top);
            flag11 = this.Factory11_.Build(FineLevel,CoarseLevel);
            if flag11 == false, flag = false; return; end;
            P11  = CoarseLevel.Get('P',this.Factory11_).GetMatrixData();
            C11Nsp = CoarseLevel.Get('NullSpace');
            
            FineN22  = FineLevel.Get('N22');
            FineLevel.Set('A', FineLevel.Get('A22'));
            FineLevel.Set('NullSpace', NS.bot);

            flag22 = this.Factory22_.Build(FineLevel,CoarseLevel);
            if flag22 == false, flag = false; return; end;
            P22  = CoarseLevel.Get('P',this.Factory22_).GetMatrixData();
            C22Nsp = CoarseLevel.Get('NullSpace');
            
            CoarseN11 = size(P11,2);
            CoarseN22 = size(P22,2);
            CoarseLevel.Set('N11', CoarseN11);
            CoarseLevel.Set('N22', CoarseN22);
            
            BigP = [P11 sparse(FineN11,CoarseN22) ; sparse(FineN22,CoarseN11) P22 ];
            
            CoarseLevel.Set('P', Operator(BigP,1,1),this);
            NewNull.top = C11Nsp;
            NewNull.bot = C22Nsp;
            CoarseLevel.Set('NullSpace', NewNull);

            FineLevel.Set('A', oldA);
            FineLevel.Set('NullSpace', oldNS);
        end
    end
end

