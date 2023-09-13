%TWOLEVELFACTORYBASE Abstract factory for building P, R, and A on coarse levels
%
% Each factory must provide a function build() which takes 2
% multigrid levels and returns the coarse level (and the fine level)
% which is now populated with additional operators (e.g., P). It is
% generally assumed that the fine level is sufficiently populated
% so that the factory is able to build the new operator on the
% coarse level.
% Note: the fine level is returned if things like aggregates are stored
%       on the finest level or a matrix diagonal is computed on the finest
%       level.
%
classdef TwoLevelFactoryBase < VerboseObject
  methods (Abstract)
     flag = Build(this, FineLevel, CoarseLevel, Specs);
     % Build method can use FineLevel and CoarseLevel information to populate an operator (or operators) on the CoarseLevel or FineLevel
     %
     % These factories are generally used for constructing grid transfers and coarse level discretizations

     SetNeeds(this, FineLevel, CoarseLevel);
     % To place request on the data of the level (increment counter
     % by calling Level.Request())
  end

  methods
  function [this] = TwoLevelFactoryBase(arg)
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
