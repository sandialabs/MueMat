classdef SaCoarseNSFactory < VerboseObject
  % A factory for generating a coarse grid representation of the
  % nullspace.
  %
  % This factory use the tentative prolongator of Smoothed
  % Aggregation to build the coarse nullspace.

  properties (Access = private)
    QR_ = true
  end % properties

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Public methods                                                            %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  methods

    function [this] = SaCoarseNSFactory(arg)
      %COARSENSFACTORY
      %
      %   SYNTAX   obj = CoarseNSFactory();
      %
      %     obj -

      % Copy constructor
      if nargin == 1 && isa(arg, class(this)), this.Copy_(arg,[]); return; end
      %
    end

    function SetNeeds(this, FineLevel, CoarseLevel)
    end

    function TentativeWithQR(this, value)
      % Specify whether orthogonalization performed when forming tentative prolongator
      this.QR_ = value;
    end

    function  [Ptent,CNull] = Build(this,AggInfo,Amat,FNull, options,OutputLevel,varargin)
      %BUILD Build a coarse grid representation of the near nullspace.
      %
      %   SYNTAX   [Ptent, CNull] = obj.Build(AggInfo, Amat, FNull, Aggs, options, OutputLevel, varargin);
      %
      %     AggInfo     - Aggregation object
      %     Amat        - an Operator object
      %     FNull       - fine nullspace
      %     options     - option structure
      %     OutputLevel - verbosity level
      %     varargin    - optional parameter
      %     Ptent       - optional tentative prolongator
      %     CNull       - coarse grid nullspace

      [Ptent,CNull] = TentativePFactory.MakeTentative(AggInfo, Amat, FNull, this.QR_, this.GetOutputLevel());
    end

  end %public methods

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

end %class CoarseNSFactory
