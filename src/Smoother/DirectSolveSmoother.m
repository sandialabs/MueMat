classdef DirectSolveSmoother < SmootherBase & SmootherPrototype
  % Direct solve smoother class
  % This smoother should not be used alone but rather inside of a
  % Hybrid2x2 Smoother.
  %
  % See also: Hybrid2x2Smoother, Hybrid2x2SmootherFactory, lu.

  properties (Access = private)
    % Parameters of the smoother
    % *no parameters*

    % Data after Setup phase: P*(R\A)*Q = L*U
    L_ % unit lower triangular matrix
    U_ % upper triangular matrix
    P_ % permutation matrix
    Q_ % permutation matrix
  end

  methods
    function [this] = DirectSolveSmoother(obj)
      %DIRECTSOLVESMOOTHER Constructor
      %
      %   SYNTAX   DirectSolveSmoother();

      % Copy constructor
      if nargin == 1 && isa(obj, class(this)), this.Copy_(obj,[]); return; end

      this.SetType('directsolve');
     end % function

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameters / Config
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function CopyParameters(this, src)
      %COPYPARAMETERS Copy the parameters of another smoother prototype.
      % See also: SmootherPrototype.Copy
      %
      %   SYNTAX   obj.CopyParameters(src);
      %
      %     src - Object (SmootherPrototype of same type)

      % *nothing to do*
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

      Amat = Level.Get('A');
      [this.L_,this.U_,this.P_,this.Q_] = lu(Amat.GetMatrixData());

      this.SetIsSetup(true);
    end % function

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Apply
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [iu,SolStatus] = Apply(this, iu, irhs, SolStatus)
      %APPLY  Apply the smoother
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

        u = this.Q_*(this.U_\(this.L_\(this.P_*rhs)));
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
      fprintf('%s%s: *no parameters*\n', prefix, this.GetType());
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
