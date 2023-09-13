classdef SmootherPrototype < CopiableHandle
% Base class for smoother prototypes
% A smoother prototype is a smoother which can be in two states:
%  - ready to be duplicated (parameters defined)
%  - ready to be used (setup phase completed)
%
% 'Smoother prototypes' can be fully copied using the Copy() method
% or the copy constructor. They can also copy the parameters of an
% other smoother object of the same type (CopyParameters()). Both
% capabilities are used by Smoother Factories.
%
% See also: SmootherBase, SmootherFactory
  properties(Access = private)
   % isConfig_   % unused - true if parameters are defined
   isSetup_ = false; % (boolean)
  end

  methods

    function [bool] = isSetup(this)
      %ISSETUP Get the state of a smoother prototype.
      %
      %   SYNTAX   bool = obj.isSetup();
      %
      %     bool - isSetup (boolean)
      bool = this.isSetup_;
    end % function

    function SetNeeds(this, Level)
      % Obtain any cross factory specifications
    end

  end % methods

  methods (Access = protected)
    function SetIsSetup(this, isSetup)
      %SETISSETUP Set the state of a smoother prototype.
      % This method should be called by Setup()
      %
      %   SYNTAX   obj.SetIsSetup(isSetup);
      %
      %     isSetup - (boolean)
      this.isSetup_ = isSetup;
    end

  end % methods
  methods (Abstract)

    CopyParameters(this, src);
      %COPYPARAMETERS Copy the parameter of another smoother prototype.
      % If the setup phase have been done already, this function
      % *must* decides if parameter modifications modify the setup
      % phase. If yes, this method must change the state of the smoother
      % prototype (isSetup=false). The same must be done for each
      % methods which modify parameters of the smoother.
      % See also : ILUSmoother
      %
      %   SYNTAX   obj.CopyParameters(src);
      %
      %     src - Object (SmootherPrototype of same type)

    Setup(this, Level, Specs)
      %SETUP Run the setup phase of the smoother.
      %
      %   SYNTAX   Amat = obj.Setup(Level, Specs);
      %
      %     Level - level of the MG hierachy (Level)
      %     Specs - specifications (CrossFactory)

  end % methods

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

end % class

% TODO: 'Specific' setup