%SINGLELEVELFACTORYBASE Abstract factory for building objects that require only a single level.
%
% Each factory must provide a function build() which takes one
% multigrid levels and returns that level
% which is now populated with additional operators (e.g., P).
%
classdef SingleLevelFactoryBase < VerboseObject
  methods (Abstract)
     flag = Build(this, CurrentLevel, Specs);
     % Build method can use CurrentLevel information to populate an operator (or operators) on the CurrentLevel

     % TODO SetNeeds(this, Level);
     % To place request on the data of the level (increment counter
     % by calling Level.Request())
  end

  methods
  function [this] = SingleLevelFactoryBase(arg)
    % Copy constructor
    if nargin == 1 && isa(arg, class(this)), this.Copy_(arg,[]); return; end
    %
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
