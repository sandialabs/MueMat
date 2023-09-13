classdef RAPXfemFactory < RAPFactory
    % Factory which builds coarse discretizations via R*A*P
    % this RAPFactory extends the basic RAP factory by
    % some Xfem specific projection stuff
    properties (Access = private)
        FCSplit_
    end

    methods
        function [this] = RAPXfemFactory(FCfact)
            % copy constructor
            if nargin == 1 && isa(FCfact,class(this)), this.Copy_(FCfact,[]); return; end

            if varexist('FCfact'), this.FCSplit_ = FCfact; end;
        end

        function SetNeeds(this, FineLevel, CoarseLevel)
            % Obtain any cross factory specifications
            SetNeeds@RAPFactory(this,FineLevel,CoarseLevel);

            FineLevel.Request('NullSpace');
            if FineLevel.IsAvailable('Arr'),        FineLevel.Request('Arr');      end;
            if FineLevel.IsAvailable('dof2node'),   FineLevel.Request('dof2node'); end;
            if FineLevel.IsAvailable('A11'),        FineLevel.Request('A11'); end;
            if FineLevel.IsAvailable('A22'),        FineLevel.Request('A22'); end;
            if CoarseLevel.IsAvailable('Arr'),      CoarseLevel.Request('Arr'); end;

            if FineLevel.IsAvailable('dof2node'),
                if isempty(this.FCSplit_), error('RAPXfemFactory needs an FCSplitting object\n'); end;
                if ~FineLevel.IsRequested('FCSplitting',this.FCSplit_),
                    this.FCSplit_.SetNeeds(FineLevel);
                end
                FineLevel.Request('FCSplitting', this.FCSplit_);
            end

        end

        function flag = Build(this,FineLevel, CoarseLevel)
            % build CoarseLevel A via (CoarseLevel R)*(FineLevel A)*(CoarseLevel P)

            % call simple standard RAPFactory.Build
            flag = Build@RAPFactory(this,FineLevel,CoarseLevel);

            % there has to be a better way to find out what methods are in here
            % only used by TestXfemLevel!
            this.ProjectMyData(FineLevel,CoarseLevel);
        end
    end
    methods  (Access = private)

        function [CoarseLevel] = ProjectMyData(this,FineLevel,CoarseLevel)
            % Projects Arr to the coarse level
            if FineLevel.IsAvailable('Arr')
                NdofsPerCnode = size(FineLevel.Get('NullSpace'),2);
                NdofsPerFnode = FineLevel.Get('Arr').GetRowMap().ConstBlkSize;
                n11 = size(FineLevel.Get('Arr').GetMatrixData(),1);
                Rdata = CoarseLevel.Get('R').GetMatrixData();
                Nrows = size(Rdata,1);
                Ncols = size(Rdata,2);
                if Ncols ~= n11,
                    Rdata = Rdata(1:Nrows-Ncols+n11,1:n11);
                    Pdata = CoarseLevel.Get('P').GetMatrixData();
                    Pdata = Pdata(1:n11,1:Nrows-Ncols+n11);
                    CoarseLevel.Set('Arr',Operator(Rdata,NdofsPerCnode,NdofsPerFnode)*...
                        FineLevel.Get('Arr')*Operator(Pdata,NdofsPerFnode,NdofsPerCnode));
                else
                    CoarseLevel.Set('Arr',CoarseLevel.Get('R')*FineLevel.Get('Arr')*CoarseLevel.Get('P'));
                end
            end
            if FineLevel.IsAvailable('dof2node'),
                if ~FineLevel.IsAvailable('FCSplitting', this.FCSplit_),
                    this.FCSplit_.Build(FineLevel);
                end
                Roots = FineLevel.Get('FCSplitting',this.FCSplit_).cpoints;
                NdofsPerCnode = size(FineLevel.Get('NullSpace'),2);
                NdofsPerFnode = FineLevel.Get('Arr').GetRowMap().ConstBlkSize;
                temp      = FineLevel.Get('dof2node');
                temp      = temp(Roots*NdofsPerFnode);
                NCdof     = size(CoarseLevel.Get('Arr').GetMatrixData,1);
                dof2node = zeros(NCdof,1);
                for i=1:NdofsPerCnode
                    dof2node(i:NdofsPerCnode:end) = temp;
                end
                CoarseLevel.Set('dof2node',dof2node);
            end
            % Projects A11 to coarse level
            FNull = FineLevel.Get('NullSpace');
            Rmat = CoarseLevel.Get('R');
            Pmat = CoarseLevel.Get('P');
            if FineLevel.IsAvailable('A11'),
                NdofPerCnode = size(FNull.top,2);
                NdofPerFnode = FineLevel.Get('A11').GetRowMap().ConstBlkSize;
                cstart = 1; fstart = 1;
                cend   = CoarseLevel.Get('N11'); fend   = FineLevel.Get('N11');
                data   = Rmat.GetMatrixData();
                Rdata  = data(cstart:cend, fstart:fend);
                data   = Pmat.GetMatrixData();
                Pdata  = data(fstart:fend, cstart:cend);
                CoarseLevel.Set('A11',Operator(Rdata,NdofPerCnode,NdofPerFnode)*...
                    FineLevel.Get('A11')*Operator(Pdata,NdofPerFnode, NdofPerCnode));
                CoarseLevel.Set('N11', size(CoarseLevel.Get('A11'),1));
            end
            % Projects A22 to coarse level
            if FineLevel.IsAvailable('A22'),
                NdofPerCnode = size(FNull.bot,2);
                NdofPerFnode = FineLevel.Get('A22').GetRowMap().ConstBlkSize;
                cstart = CoarseLevel.Get('N11')+1;        fstart= FineLevel.Get('N11')+1;
                cend  = cstart+CoarseLevel.Get('N22')-1; fend= fstart+FineLevel.Get('N22')-1;
                data  = Rmat.GetMatrixData();
                Rdata = data(cstart:cend, fstart:fend);
                data  = Pmat.GetMatrixData();
                Pdata = data(fstart:fend, cstart:cend);
                CoarseLevel.Set('A22',Operator(Rdata,NdofPerCnode,NdofPerFnode)*...
                    FineLevel.Get('A22')*Operator(Pdata,NdofPerFnode,NdofPerCnode));
                CoarseLevel.Set('N22', size(CoarseLevel.Get('A22'),1));
            end

            % release variables
            if FineLevel.IsAvailable('NullSpace'),  FineLevel.Release('NullSpace'); end;
            if FineLevel.IsAvailable('Arr'),        FineLevel.Release('Arr'); end;
            if FineLevel.IsAvailable('dof2node'),   FineLevel.Release('dof2node'); end;
            if FineLevel.IsAvailable('A11'),        FineLevel.Release('A11'); end;
            if FineLevel.IsAvailable('A22'),        FineLevel.Release('A22'); end;
            if CoarseLevel.IsAvailable('Arr'),      CoarseLevel.Release('Arr'); end;
            if FineLevel.IsAvailable('FCSplitting',this.FCSplit_), FineLevel.Release('FCSplitting',this.FCSplit_); end;
        end

    end % methods

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
