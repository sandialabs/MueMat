% TODO: should be a singleton

classdef NoFactory < VerboseObject
    properties (Access = private)
        label_ = 'no label'
    end

    methods
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

        function [this] = NoFactory(label)
            % Copy constructor
            if nargin == 1 && isa(label, class(this)), this.Copy_(label,[]); return; end
            %

            if varexist('label'), label_ = label; end
        end

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

        function SetNeeds(this, FineLevel, CoarseLevel)

        end

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

        function flag = Build(this, FineLevel, CoarseLevel)
            fprintf('NoFactory::Build() error');
            FineLevel.Print();
            CoarseLevel.Print();
            % label_
            keyboard;
        end

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
        % FIXME: Remove. Needed now for using NoFactory as a PatternFactory
        % in EminPFactory
        function SetInitialPFactory(this, InitPFact)
        end
        function [InitPFact] = GetInitialPFactory(this)
            InitPFact = [];
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