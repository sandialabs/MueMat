% A factory class that translates aggregates from an aggregation factory to
% a FC splitting object
classdef FCSplittingFactory < VerboseObject
    properties (Access = private)
        AggFact_
    end

    methods
        function [this] = FCSplittingFactory(AggFact)
            % Copy constructor
            if nargin == 1 && isa(AggFact, class(this)), this.Copy_(AggFact,[]); return; end;

            % create an FCSplitting factory
            if varexist('AggFact'), this.AggFact_ = AggFact;
            else this.AggFact_ = 'default'; end;
        end

        function SetNeeds(this, CurrentLevel)
            if CurrentLevel.IsAvailable('FCsplitting',this), return; end;

            if ~CurrentLevel.IsRequested('Aggregates',this.AggFact_),
                this.AggFact_ = CurrentLevel.InterpretHandle('Aggregates',this.AggFact_);
                this.AggFact_.SetNeeds(CurrentLevel);
            end;
            CurrentLevel.Request('Aggregates',this.AggFact_);
        end

        function flag = Build(this, CurrentLevel, Specs)
            flag = true;

            if ~CurrentLevel.IsAvailable('Aggregates',this.AggFact_)
                fprintf('FCSplittingFactory: build aggregates\n');
                this.AggFact_.Build(CurrentLevel);
            end

            AggInfo = CurrentLevel.Get('Aggregates',this.AggFact_);
            CurrentLevel.Release('Aggregates',this.AggFact_);

            fcsplitting.cpoints = AggInfo.Roots;
            temp = ones(size(AggInfo.AggId,1),1);
            temp(AggInfo.Roots) = 0;
            temp = find(temp);
            fcsplitting.fpoints = temp;
            clear temp;

            CurrentLevel.Set('FCSplitting', fcsplitting,this);
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