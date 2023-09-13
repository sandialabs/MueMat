classdef ILUSmoother < SmootherBase & SmootherPrototype
  % ILU smoother class
  % This smoother applies an incomplete factorization to
  % the residual equation. The ILU factors are computed by Matlab's ILU.
  % The default is ILU(0), but an optional parameter struct can be passed in
  % to set any option accepted by Matlab's ILU.
  %
  % See also: ilu

  properties (Access = private)
    % Parameters of the smoother
    Params_ % parameters for Matlab's ILU

    % Data after Setup phase
    A_      % matrix of the level
    L_      % factor L of A=LU (unit lower triangular matrix)
    U_      % factor U of A=LU (upper triangular matrix)
    P_      % permutation P
  end

  methods
    function [this] = ILUSmoother(Params)
      %ILUSMOOTHER Constructor
      %
      %   SYNTAX   ILUSmoother(Params);
      %
      %     Params - optional parameters for MATLAB's ILU (struct,optional,default=Params.type='nofill')

      % Copy constructor
      if nargin == 1 && isa(Params, class(this)), this.Copy_(Params,[]); return; end
      %

      this.Params_.type = 'nofill'; % ILU(0) is default
      if varexist('Params')
        if ~isempty(Params),  this.SetParameters(Params);    end
      end

      this.SetType('ILU');
    end % function

    function SetParameters(this, Params)
      %SETPARAMETERS Set smoother parameters.
      %
      %   SYNTAX   obj.SetParameters(Params);
      %
      %     Params - optional parameters for MATLAB's ILU (struct)
      %
      % Note: if parameters of ILU are changed, Setup() must be done again.

      if ~structeq(this.Params_, Params)
        this.A_ = []; this.L_ = []; this.U_ = []; this.P_ = [];
        this.SetIsSetup(false);
      end

      this.Params_ = Params;
    end % function

    function [Params] = GetParameters(this)
      %GETPARAMETERS Get ILU parameters
      %
      %   SYNTAX   Params = obj.GetParameters();
      %
      %     Params - optional parameters for MATLAB's ILU (struct)
      Params = this.Params_;
    end % function

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameters / Config
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function CopyParameters(this, src)
      %COPYPARAMETERS Copy the parameters of another smoother prototype.
      % See also: SmootherPrototype.CopyParameters, ILUSmoother.SetParameters
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

      this.A_ = Level.Get('A');
      [this.L_,this.U_,this.P_] = ilu(this.A_.GetMatrixData(),this.Params_);

      this.SetIsSetup(true);
    end % function

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Apply
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [iu,SolStatus] = Apply(this, iu, irhs, SolStatus)
      %APPLY Apply the smoother
      % See also: SmootherBase.Apply
      %
      %   SYNTAX   [iu, SolStatus] = obj.Apply(iu, irhs, SolStatus);
      %
      %     iu        - vector to smooth (could be a MultiVector or in Matlab format)
      %     irhs      - right-hand side  (could be a MultiVector or in Matlab format)
      %     SolStatus - when InitGuessStatus==ALLZEROS, iu is assumed to be all zeros and initial matvec operation could be avoid (optional,default=NOTALLZEROS)

      if ~this.isSetup(), error('apply'); end

      if   isa(iu,'MultiVector'), u = iu.GetVectors();
      else                        u = iu;                   end
      if   isa(irhs,'MultiVector'), rhs = irhs.GetVectors();
      else                          rhs = irhs;             end

      mue_include

      % Solve correction equations
      if (SolStatus == ALLZEROS)
        ue = this.U_\(this.L_\(this.P_*rhs));
      else
        ue = this.U_\(this.L_\(this.P_*(rhs-this.A_*u)));
      end
      u = u + ue;

      SolStatus = NOTALLZEROS;

      if   isa(iu,'MultiVector'), iu.SetVectors(u);
      else                        iu = u;                   end
    end % function

    function Print(this,prefix)
      %PRINT Print smoother information
      %
      %   SYNTAX   obj.Print()
      %
      %     prefix  - optional string that is prepended to each print line

      if ~varexist('prefix'), prefix = ''; end

      % Get the description of the struct Params_
      str = evalc('disp(this.Params_)');

      % Format output (, and ':' per '=')
      str = strrep(str, '  ', '');            % - skip some blank space
      str = strrep(str, sprintf('\n'), ', '); % - replace \n per ','
      str = strrep(str, ', ,', '');           %   and handle final \n
      str = strrep(str, ': ', '=');           % - replace ': ' per '='

      % Print
      fprintf('%s%s: %s\n', prefix, this.GetType(), str);
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
    end

  end % methods

end % class

% TODO:
% P depends on When SETUP.type == 'nofill'
% or when SETUP.type == 'crout', P is always an identity matrix as
% neither of these methods performs pivoting.  When SETUP.type =
% 'ilutp', the role of P is determined by the value of SETUP.milu.  When
% SETUP.milu ~= 'row', P is returned such that L and U are incomplete
% factors of P*A.  When SETUP.milu == 'row', P is returned such that L
% and U are incomplete factors of A*P.
