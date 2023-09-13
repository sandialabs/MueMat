classdef RAPFactory < TwoLevelFactoryBase
    % Factory which builds coarse discretizations via R*A*P
    % this is the basic RAP factory, that just computes A_c = R*A*P
    % no projection of any other data
    % if you have some variables to be projected to coarser levels, use and
    % extend RAPexFactory or implement your own problem-specific RAPFactory
    % derived class (see RAPXfemFactory as example)
    methods
        function [this] = RAPFactory()
        end

        function SetNeeds(this, FineLevel, CoarseLevel)
            % Obtain any cross factory specifications
        end

        function flag = Build(this,FineLevel, CoarseLevel)
            % build CoarseLevel A via (CoarseLevel R)*(FineLevel A)*(CoarseLevel P)
            flag = true;

            if ~CoarseLevel.IsAvailable('P'), fprintf('No Pmat??\n'); keyboard; end
            if ~CoarseLevel.IsAvailable('R'), fprintf('No Rmat??\n'); keyboard; end
            if ~FineLevel.IsAvailable('A'), fprintf('No Amat??\n'); keyboard; end
            Pmat = CoarseLevel.Get('P');
            Rmat = CoarseLevel.Get('R');
            Amat = FineLevel.Get('A');
            CoarseLevel.Set('A', Rmat*(Amat*Pmat));

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
