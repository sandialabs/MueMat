%% GenericPRFactory
% This class combines a prolongation operator object (derived from
% PFactory) and a restriction operator object (derived from RFactory) and
% creates a PRFactory, that handles both the prolongator and the
% restrictor.
%%

classdef GenericPRFactory < PRFactory
    % This class combines a prolongation operator object (derived from
    % PFactory) and a restriction operator object (derived from RFactory) and
    % creates a PRFactory, that handles both the prolongator and the
    % restrictor.

    properties (Access = private)
        PFact_ = [];                % internal PFactory derived object
        RFact_ = [];                % internal RFactory derived object
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% public functions
    methods
        function [this] = GenericPRFactory(PFact, RFact)
            %GENERICPRFACTORY Constructor
            %
            %   SYNTAX obj = GenericPRFactory(PFact, RFact)
            %
            %    default parameters for GenericPRFactory constructor:
            %    PFact = TentativePFactory()
            %    RFact = TransPFactory()
            %       -> plain aggregation, no smoothing

            % copy constructor
            if nargin == 1 && isa(PFact,class(this))
                this.Copy_(arg,[]); fprint('Copy constructor\n'); return;
            end;

            if varexist('PFact') this.PFact_ = PFact;
            else this.PFact_ = TentativePFactory(); end;

            if varexist('RFact') this.RFact_ = RFact;
            else this.RFact_ = TransPFactory(); end;

        end

        function SetNeeds(this, FineLevel, CoarseLevel)
          % Obtain any cross factory specifications
          if ~isempty(this.PFact_),
            this.PFact_.SetNeeds(FineLevel, CoarseLevel);
          end
          if ~isempty(this.RFact_),
            this.RFact_.SetNeeds(FineLevel, CoarseLevel);
          end

          % requests for prolongation and restriction operator (results)
          CoarseLevel.Request('P',this.PFact_);
          CoarseLevel.Request('R',this.RFact_);
        end

        function flag = Build(this,FineLevel,CoarseLevel)
            % Build transfer operators (default implementation)
            %
            %   SYNTAX flag = Build(FineLevel, CoarseLevel, Specs)
            %
            %      FineLevel    - FineLevel (Level)
            %      CoarseLevel  - CoarseLevel (Level)

            flag = true;

            % check for break condition (max coarse size)
            if FineLevel.Get('A').GetRowMap().NNodes() <= this.MaxCoarseSize_, flag = false; return; end

            % compute and set prolongator (calls CoarseLevel.SetP)
            this.PFact_.Build(FineLevel,CoarseLevel);

            % store result as general P for the multigrid method
            CoarseLevel.Set('P', CoarseLevel.Get('P', this.PFact_));
            CoarseLevel.Release('P', this.PFact_);


            % compute and set restrictor (calls CoarseLevel.SetR)
            this.RFact_.Build(FineLevel,CoarseLevel);

            % store the result as general R for the multigrid method
            CoarseLevel.Set('R', CoarseLevel.Get('R', this.RFact_));
            CoarseLevel.Release('R', this.RFact_);


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
            [cmd, dummy, mc] = this.CopyCmd_(src,mc);
            eval(cmd);
        end

    end % methods
end
