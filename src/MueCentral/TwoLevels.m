classdef TwoLevels < VerboseObject
    properties (Access = private)
        FineLevel_
        CoarseLevel_
    end

    methods
        function this = TwoLevels(FineLevel, CoarseLevel)
            if nargin == 1 && isa(FineLevel, class(this)), this.Copy_(FineLevel,[]); return; end
            this.FineLevel_ = FineLevel;
            this.CoarseLevel_ = CoarseLevel;
        end

        % Indicate that an object is needed. This increments the storage counter.
        % This method is used during the pre-setup phase.
        function Request(this, ename, ehandle, Level)
            if Level ~= this.FineLevel_ && Level ~= this.CoarseLevel_
                error('TwoLevels.Request: Level parameter must be either FineLevel or CoarseLevel.\n');
            end

            % translate ehandle (e.g. 'default' to default factory)
            ehandle = Level.InterpretHandle(ename,ehandle);

            % call SetNeeds
            if ~Level.IsRequested(ename, ehandle)
                if ~ismethod(ehandle,'SetNeeds'),
                    error('TwoLevels.Request: ehandle must be a factory\n');
                end
                ehandle.SetNeeds(this.FineLevel_,this.CoarseLevel_);
            end
            Level.Request(ename,ehandle);

        end % Request()

        % Get data. This does not decrement the storage counter.
        function data = Get(this, ename, ehandle, Level, forcebuild)
            if Level ~= this.FineLevel_ && Level ~= this.CoarseLevel_
                error('TwoLevels.Get: Level parameter must be either FineLevel or CoarseLevel.\n');
            end
            if ~varexist('forcebuild'), forcebuild = 0; end;

            % translate ehandle (e.g. 'default' to default factory)
            ehandle = Level.InterpretHandle(ename,ehandle);

            if ~Level.IsAvailable(ename,ehandle) || forcebuild == 1
                if ~ismethod(ehandle,'Build'),
                    error('TwoLevels.Get: ehandle must be a factory\n');
                end
                ehandle.Build(this.FineLevel_, this.CoarseLevel_);
            end

            data = Level.Get(ename,ehandle);
        end % Get()
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