classdef TestSmoother < SmootherBase & SmootherPrototype
  % Test smoother class
  % This smoother does *nothing* and exists only for test purpose.
  % It allows us to test the 'Setup Reuse' capability of
  % SmootherFactory, Hybrid2x2Smoother and MergedSmoother
  %
  % This smoother can change is type name and have a fake parameter
  % that change the setup phase.
  %
  % See also: SmootherReUseSetupTest TODO
  properties (Access = private)

  Param_ = 0 % fake parameter that change the setup phase (integer)
  
  SetupDone_ = false % true if the setup phase have been done
                     % specifically for this smoother (== no reuse)
  end
  
  methods
    function [this] = TestSmoother(Type, Param)
      %TESTSMOOTHER Constructor
      % 
      %   SYNTAX   TestSmoother(Param);
      % 
      %     type  - Name of the fake smoother
      %     param - Smoother parameter (integer)
      
      % Copy constructor
      if nargin == 1 && isa(Type, class(this)), this.Copy_(Type,[]); return; end
      %
      
      if varexist('Type'), this.SetType(Type); end;
      if varexist('Param'), this.Param_ = Param; end;
    end % function
    
    function SetParameters(this, Param)
      %SETPARAMETERS Set smoother parameter
      % 
      %   SYNTAX   obj.SetParameters(Param);
      % 
      %     Param - Smoother parameter (integer)
      %
      if (this.Param_ ~= Param)
        this.Param_ = Param;
        this.SetIsSetup(false);
      end

    end % function
    
    function [Param] = GetParameters(this)
      %GETPARAMETERS Get smoother parameter
      % 
      %   SYNTAX   Param = obj.GetParameters();
      %
      %     Param - Smoother parameter (integer)
      Param = this.Param_;
    end % function
  
    function [SetupDone] = SetupDone(this)
      %GETPARAMETERS Get smoother parameter
      % 
      %   SYNTAX   SetupDone = obj.SetupDone();
      %
      %     SetupDone - true if the setup phase have been done
      SetupDone = this.SetupDone_;
    end % function
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameters / Config
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function CopyParameters(this, src)
      %COPYPARAMETERS Copy the parameters of another smoother prototype.
      % See also: SmootherPrototype.CopyParameters
      % 
      %   SYNTAX   obj.CopyParameters(src);
      % 
      %     src - Object (SmootherPrototype of same type)
      this.SetParameters(src.GetParameters());
    end % function
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setup phase
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function Setup(this, Level, Specs)
      %SETUP Run the setup phase of the smoother.
      % See also: SmootherPrototype.Setup
      % 
      %   SYNTAX   Amat = obj.Setup(Level, Specs);
      % 
      %     Level - level of the MG hierachy (Level)
      %     Specs - specifications (CrossFactory)
      if this.isSetup(), return; end

      % fprintf('Setup\n');
      this.SetupDone_ = true;
      
      this.SetIsSetup(true);
    end % function
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Apply
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [iu,SolStatus] = Apply(this, iu, irhs, SolStatus)

     error('Cannot apply a TestSmoother because it is just an empty shell for tests.');

    end % function

    function Print(this,prefix)
      %PRINT Print smoother information
      %
      %   SYNTAX   obj.Print()
      %
      %     prefix  - optional string that is prepended to each print line

      if ~varexist('prefix'), prefix = ''; end
      
      % Print
      fprintf('%s%s: name=''%s'', param=%d\n', prefix, 'TestSmoother', this.GetType(), param);
    end %Print()

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
      
      %% IMPORTANT:
      this.SetupDone_ = false;
    end
    
  end % methods
  
end % class
