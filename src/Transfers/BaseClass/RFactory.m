%% RFactory
% This class provides an interface for a concrete implementation of a
% restriction operator
%%

classdef RFactory < VerboseObject
    % This factory provides an interface for a concrete implementation of a
    % restriction operator
    % For a concrete implementation the user has to overwrite the virtual
    % Build method

    properties (Access = private)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% public functions
    methods
        function [this] = RFactory(arg)
        %RFACTORY Constructor

            % copy constructor
            if nargin == 1 && isa(arg,class(this)), this.Copy_(arg,[]); return; end;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% abstract functions
    methods (Abstract)
        flag = Build(this, FineLevel, CoarseLevel, Specs)

        SetNeeds(this, FineLevel, CoarseLevel);
        % To place request on the data of the level (increment counter
        % by calling Level.Request())
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% protected functions
    methods (Access = protected)

        function Copy_(this,src,mc)
            % COPY_
            % syntax obj.Copy_(src,mc);
            % src: object to copy
            % mc: MATLAB Metaclass
            [cmd,data,mc] = this.CopyCmd_(src,mc);
            eval(cmd);
        end
    end

end
