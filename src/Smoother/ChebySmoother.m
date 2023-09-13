classdef ChebySmoother < SmootherBase & SmootherPrototype
  % Chebyshev smoother class

  properties (Access = private)
    % Parameters of the smoother
    nIts_         = 1         % number of iterations
    LambdaRatio_  = 1/27      % ratio between maximum and minimum eigenvalues
    diagonalView_ = 'current' % diagonal view label (default == current view)

    % Data after Setup phase
    Amat_          % matrix of the level
    lambda_1_      % maximum eigenvalue
    % lambda_n_    % unused
  end

  methods
    function [this] = ChebySmoother(nIts, LambdaRatio, diagonalView)
      %CHEBYSMOOTHER Constructor
      %
      %   SYNTAX   ChebySmoother(nIts, LambdaRatio);
      %
      %     nIts        - number of iterations
      %     LambdaRatio - ratio between maximum and minimum eigenvalue

       % Copy constructor
       if nargin == 1 && isa(nIts, class(this)), this.Copy_(nIts,[]); return; end
       %

       %if nargout() < 1, error('nargout'); end
       if varexist('nIts'),         this.SetNIts(nIts);                 end
       if varexist('LambdaRatio'),  this.SetLambdaRatio(LambdaRatio);   end
       if varexist('diagonalView'), this.SetDiagonalView(diagonalView); end

       this.SetType('chebychev');
     end % function

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameters / Config
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function SetNIts(this, nIts)
      %SETNITS Set the number of iterations
      %
      %   SYNTAX   obj.SetNIts(nIts);
      %
      %     nIts - number of iterations
      this.nIts_ = nIts;
    end % function

    function [nIts] = GetNIts(this)
      %GETNITS Get the number of iterations
      %
      %   SYNTAX   nIts = obj.GetNIts();
      %
      %     nIts - number of iterations
      nIts = this.nIts_;
    end % function

    function SetLambdaRatio(this, LambdaRatio)
      %SETLAMBDARATIO Set the ratio between maximum and minimum eigenvalue
      %
      %   SYNTAX   obj.SetLambdaRatio(LambdaRatio);
      %
      %     LambdaRatio - ratio between maximum and minimum eigenvalue
      this.LambdaRatio_ = LambdaRatio;
    end % function

    function [LambdaRatio] = GetLambdaRatio(this)
      %GETLAMBDARATIO Get the ratio between maximum and minimum eigenvalue
      %
      %   SYNTAX   LambdaRatio_ = obj.GetLambdaRatio();
      %
      %     LambdaRatio_ - ratio between maximum and minimum eigenvalue
      LambdaRatio = this.LambdaRatio_;
    end % function

    function SetDiagonalView(this, diagonalView)

      if ~strcmp(this.diagonalView_, diagonalView)
        this.Amat_ = []; this.lambda_1_ = [];
        this.SetIsSetup(false);
      end

      this.diagonalView_ = diagonalView;
    end % function

    function [diagonalView] = GetDiagonalView(this)
      diagonalView = this.diagonalView_;
    end % function

    function CopyParameters(this, src)
      %COPYPARAMETERS Copy the parameters of another smoother prototype.
      % See also: SmootherPrototype.CopyParameters
      %
      %   SYNTAX   obj.CopyParameters(src);
      %
      %     src -  Object (SmootherPrototype of same type)
      this.SetNIts(src.GetNIts());
      this.SetLambdaRatio(src.GetLambdaRatio());
      this.SetDiagonalView(src.GetDiagonalView()); %TODO isSetup()
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
      % unneeded: previousView = Amat.SwitchToView(this.diagonalView_);

      %
      % diagonal    = Amat.GetDiagonal();
      % ApplyInvFcn = diagonal.GetApplyInverse();

      %
      this.Amat_     = Amat;
      this.lambda_1_ = Amat.GetDinvALambda(this.diagonalView_); % Estimate of maximum eigenvalue
      % this.lambda_n_ = this.LambdaRatio_ * this.lambda_1_;    % Estimate of minimum eigenvalue
                                                                                    % = lambda_1/27 or lambda_1/30

      % scalar diag (for tests)
      if (0)
       DinvA = Amat.GetMatrixData() / diag(diag(Amat.GetMatrixData()));
       this.lambda_1_ = eigs(DinvA,1);               % Estimate of maximum eigenvalue
       this.lambda_n_ = LambdaRatio * this.lambda_1_;% Estimate of minimum eigenvalue
      end

      % unneeded: Amat.SwitchToView(previousView);

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

       Amat = this.Amat_;
       previousView = Amat.SwitchToView(this.diagonalView_);

       if   isa(iu,'MultiVector'), u = iu.GetVectors();
       else                        u = iu;                   end
       if   isa(irhs,'MultiVector'), rhs = irhs.GetVectors();
       else                          rhs = irhs;             end

       lambda_n = this.LambdaRatio_ * this.lambda_1_;

       [u,SolStatus] = ChebySmoother.cheby(Amat, u, SolStatus, rhs, this.nIts_, this.lambda_1_, lambda_n);
       % note: x0 = u; b  = rhs;

       if   isa(iu,'MultiVector'), iu.SetVectors(u);
       else                        iu = u;                   end

       Amat.SwitchToView(previousView);
    end% function

    function Print(this,prefix)
      %PRINT Print smoother information
      %
      %   SYNTAX   obj.Print()
      %
      %     prefix  - optional string that is prepended to each print line

      if ~varexist('prefix'), prefix = ''; end
      fprintf('%sChebyshev: degree=%d, lambda_max=%g\n',...
               prefix,this.nIts_,lambda_1_);

      % TODO: complete this method
    end %Print()

  end % methods

  methods (Access = private, Static = true)
    function [x, SolStatus]=cheby(A,x,SolStatus,b,nits,lambda_1,lambda_n)
      %CHEBY Runs a Chebychev polynomial method
      % on the matrix A to solve the equation Ax=b.Uses initial guess x0.
      % This method is an implemenation of Algorithm 12.1 in
      % Saad's "Iterative Methods for Sparse Linear Systems".
      %
      %   SYNTAX   [x, SolStatus] = obj.cheby(A, x, SolStatus, b, nits, lambda_1, lambda_n);
      %
      %     A         - Matrix of the system (MueMat Operator)
      %     x         - Initial guess / Solution vector from Chebyshev
      %     SolStatus - NOTALLZEROS if x != 0. When InitGuessStatus==ALLZEROS, iu is assumed to be all zeros and initial matvec operation could be avoid (optional,default=NOTALLZEROS,output=NOTALLZEROS)
      %     b         - RHS of system
      %     nits      - Number of iterations to perform (defaults to 1). nits is the degree of the chebyshev polynomial.
      %     lambda_1  - Estimate of maximum eigenvalue
      %     lambda_n  - Estimate of minimum eigenvalue. For a multigrid smoother, using lambda_1/30 is probably fine.
      %
      % Version History
      % 0.3 : 06/28/10 - Diagonal scaling <jngaida@sandia.gov>
      % 0.2 : 05/05/10 - Modified for MueMat <jngaida@sandia.gov>
      % 0.1 : 08/02/06 - Initial Version <csiefer@sandia.gov>

      mue_include

       % deg == 0
       if (nits == 0), return; end

       % Get diagonal
       diagonal = A.GetDiagonal();                 % blk diag
       ApplyInvFcn = diagonal.GetApplyInverse();   %
       if (0)
        diagonal = diag(A.GetMatrixData());       % scalar diag (for tests)
        zerosEntries = find( diagonal == 0.);     %
        diagonal(zerosEntries) = ones(size(zerosEntries,1),1);
        ApplyInvFcn = @ScalarDiagApplyInvFcn;
       end
       % debug diag view: size(diagonal.diagData_)

       % Constants
       theta = (lambda_1 + lambda_n) / 2;
       delta = (lambda_1 - lambda_n) / 2;
       sigma = theta / delta;

       % Initialize
       rho_prev = 1 / sigma;

       % deg == 1
       % This function has an ability to avoid a matvec if the initial
       % guess is zero
       if (SolStatus == NOTALLZEROS)
        % Compute initial residual norm
        r        = (b - A*x);
        scaled_r = ApplyInvFcn(diagonal, r); % = r ./ diagonal;
        d        = scaled_r / theta;
        x        = x + d; % Update solution
       else
        % Avoid matvec when initial guess is all zeros
        r        = b;
        scaled_r = ApplyInvFcn(diagonal, r); % = r ./ diagonal;
        d        = scaled_r / theta;
        x        = d;     % Update solution
       end

       SolStatus = NOTALLZEROS;

       % deg > 1
       for i=1:nits-1,
        % Matvec, update residual, constants, next step
        r        = (r - A*d);
        scaled_r = ApplyInvFcn(diagonal, r); % = r ./ diagonal;
        rho      = 1 / (2*sigma - rho_prev);
        d        = rho * rho_prev * d   +   ( (2*rho/delta) * scaled_r );
        rho_prev = rho;

        % Update solution
        % at the end of the loop (!= Saad's book)
        x = x + d;
       end

    end % function
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

% function [vec] = ScalarDiagApplyInvFcn(diagonal, vec)
%       %SCALARDIAGAPPLYINVFCN Scalar diag inverse apply function (for tests)
%       %
%       %   SYNTAX   vec = obj.ScalarDiagApplyInvFcn(diagonal, vec);
%       %
%       %     diagonal - diagonal (array)
%       %     vec      - vector (array)
%   vec = vec ./ diagonal;
% end

% TODO: diagonalView = label or View ?
