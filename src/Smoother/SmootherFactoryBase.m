classdef SmootherFactoryBase < VerboseObject
  % Base class for smoother factories
  % with only one real capability which is to record (for Hierarchy class)
  % whether a smoother should be built on the coarsest level.
  %
  % See also SmootherFactory, Hierarchy
  %
  methods (Abstract)
    % Build method must be provided by subclasses.
    [PreSmoo,PostSmoo,Amat] = Build(this, Level, Specs);

    %TODO SetNeeds(this, Level);
    % To place request on the data of the level (increment counter
    % by calling Level.Request())
  end

  methods
    function [this] = SmootherFactoryBase(obj)
      %SMOOTHERFACTORYBASE Constructor
      %
      %   SYNTAX   SmootherFactoryBase();
      %

      %Copy constructor
      if nargin == 1 && isa(obj, class(this)), this.Copy_(obj,[]); return; end
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
end % classdef
