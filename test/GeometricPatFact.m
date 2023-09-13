%% GeometricPatFact
% concrete implementation of sparsity pattern factory for prolongation and
% restricton operator

classdef GeometricPatFact < PatternFactory
    properties (Access = private)
    end

    methods
        function [this] = GeometricPatFact(filterType, InitPFact, options)
            % GeometricPatFact constructor
            %
            %   SYNTAX obj = GeometricPatFact(filterType,options)
            %
            %     filterType - optional filtering (default = [])
            %     InitPFact  - PFactory object for initial guess 'Ptent'
            %     options    - options stored in options_ (used for
            %                  fattening)
            %
            % AP_Pattern Constructor: default pattern A*P

            % copy constructor
            if nargin == 1 && isa(filterType, class(this)), this.Copy_(filterType,[]); return; end
             this.type_ = 'default';

        end
        function SetNeeds(this, FineLevel, CoarseLevel)
        % Obtain any cross factory specifications
           SetNeeds@PatternFactory(this,FineLevel,CoarseLevel);

           FineLevel.Request('xcoords');
           FineLevel.Request('ycoords');
           FineLevel.Request('DomainList');
           FineLevel.Request('CrossPointList');
           FineLevel.Request('DroppedCrossPoints');
        end
%         function Build(this, FineLevel, CoarseLevel)
%           Pattern = BuildPattern(this, FineLevel, CoarseLevel);
%           CoarseLevel.Set('Ppattern', Pattern);
%           Ncoarse = size(Pattern,2);
%           CoarseLevel.Set('NullSpace', ones(Ncoarse,1));
%         end

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

          Amatrix       = FineLevel.Get('A');
          xcoords       = FineLevel.Get('xcoords');
          ycoords       = FineLevel.Get('ycoords');
          DomainList    = FineLevel.Get('DomainList');
          CrossPointList= FineLevel.Get('CrossPointList');
          dcp = FineLevel.Get('DroppedCrossPoints');
%          Aggs          = FineLevel.Get('Aggregates');
          n      = Amatrix.GetRowMap.NDOFs;
          %Ncoarse= CrossPointList.NCrossPoints;
          Ncoarse= CrossPointList.NCrossPoints - nnz(dcp);
          Pattern= spalloc(n,Ncoarse,10*n);
          IJV = ones(10*n,3);

          % could be slow
          myAgg = 1;
          offset = 1;
          for i=1:CrossPointList.NCrossPoints
             if dcp(i) == 0
               inside   = ones(n,1);
               mydomains= CrossPointList.Domains(i,:);
               minvalue = min(DomainList.LowerLeftCorner(mydomains,1));
               maxvalue = max(DomainList.UpperRightCorner(mydomains,1));
               inside(find(xcoords <  minvalue)) = 0;
               inside(find(xcoords >= maxvalue)) = 0;
               if DomainList.dim > 1,
                 minvalue= min(DomainList.LowerLeftCorner(mydomains,2));
                 maxvalue= max(DomainList.UpperRightCorner(mydomains,2));
                 inside(find(ycoords <  minvalue)) = 0;
                 inside(find(ycoords >= maxvalue)) = 0;
               end
               indices = find(inside);
               %Pattern(indices,myAgg) = 1;
               IJV(offset:offset+length(indices)-1,1) = indices;
               IJV(offset:offset+length(indices)-1,2) = myAgg;
               IJV(offset:offset+length(indices)-1,3) = 1;
               offset = offset + length(indices);
               myAgg = myAgg+1;
             end
          end
%          Pattern(Aggs.Roots,:) = speye(Ncoarse,Ncoarse);
          IJV = IJV(1:offset-1,:);
          Pattern = spconvert(IJV);

           FineLevel.Release('xcoords');
           FineLevel.Release('ycoords');
           FineLevel.Release('DomainList');
           FineLevel.Release('CrossPointList');

        end %Build()
        
        
    end % end public methods

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
