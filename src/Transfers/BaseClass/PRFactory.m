%% PRFactory
% This class provides an interface for multigrid transfer operators (both
% restriction and prolongation).
% The user has to implement the Build function. The default implementation
% is GenericPRFactory.
%%

classdef PRFactory < TwoLevelFactoryBase
  % This class provides an interface for multigrid transfer operators (both
  % restriction and prolongation).
  % The user has to implement the Build function. The default implementation
  % is GenericPRFactory.

  properties (Access = protected)

    MaxCoarseSize_   = 50;        % max size of coarsest level before stopping aggregation (default: 50)

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% public functions
  methods
    function [this] = PRFactory(arg)
      % copy constructor
      if nargin == 1 && isa(arg,class(this)), this.Copy_(arg,[]); return; end;
    end

    function SetMaxCoarseSize(this, MaxCoarseSize)
      % Determines size of coarse level when further coarsening should stop
      this.MaxCoarseSize_ = MaxCoarseSize;
    end

  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% abstract functions
  methods (Abstract)
    SetNeeds(this, FineLevel, CoarseLevel);
    flag = Build(this, FineLevel, CoarseLevel, Specs);
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