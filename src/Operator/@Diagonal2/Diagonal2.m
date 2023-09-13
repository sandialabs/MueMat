classdef Diagonal2 < OperatorView
%DIAGONAL2 Class that holds data for diagonals (point or block)
% Diagonal2 stores the block diagonal as a MATLAB sparse matrix
% (the original Diagonal class stores the block diagonal as vectors)

  properties
    diagData_
    invDiagData_

    DinvALambda_ = []

    Collection_ = []
  end %properties

  methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Constructor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [this] = Diagonal2(Amat,Collection)
      this@OperatorView([],[],[]);

      % Copy constructor
      if nargin == 1 && isa(Amat, class(this)), this.Copy_(Amat,[]); return; end
      %

      if ~varexist('Collection') || isempty(Collection)
        this.ExtractBlkDiag(Amat);
      else
        fprintf('This case is not handle (yet) by Diagonal2. Use the class Diagonal instead\n');
        keyboard;
      end

    end

    function SetDinvALambda(this,lambda)
      % Set eigenvalue estimate for diagonally-scaled matrix.
      this.DinvALambda_ = lambda;
    end

    function lambda = GetDinvALambda(this)
      % Get eigenvalue estimate for diagonally-scaled matrix.
      lambda = this.DinvALambda_;
    end

    function [diagData] = GetMatrixData(this)
      % Return underlying matrix data.
      diagData = this.diagData_;
    end

    function [invDiagData] = GetDiagonalInvData(this)
      % Return apply-inverse data.
      invDiagData = this.invDiagData_;
    end

    function [Collection] = GetCollection(this)
      % Return apply-inverse data.
      Collection = this.Collection_;
    end

    function  [Cmat] = mldivide(Amat,Bmat)
      %MLDIVIDE Operator solve.
      %
      %   SYNTAX   Cmat = Amat \ Bmat;
      %
      %     Amat - a Diagonal
      %     Bmat - an Operator
      %     Cmat - an Operator
      if (~strcmp(func2str(Amat.GetApplyInverse()),func2str(@DinvBlkApply))),
        fprintf('MatMatSolve:: Amat must be a Dinv-style matrix.\n');
        keyboard;
      end

      ApplyInverse = Amat.GetApplyInverse();

      if isa(Bmat,'Operator')
        if ~strcmp(func2str(Bmat.GetApply()),func2str(@MatlabApply)),
          fprintf('MatMatSolve:: Bmat must be a MatlabApply matrix.\n');
          keyboard;
        end

        Cdata = sparse(ApplyInverse(Amat, Bmat.GetMatrixData(), []));
        Cmat = Operator(Cdata, Amat.GetRowMap(), Bmat.GetColMap(), @MatlabApply);

      else
        Cmat = ApplyInverse(Amat, Bmat, []);
        if issparse(Bmat)
          Cmat = sparse(Cmat);
        end
      end

    end % Operator solve.

  end %methods

  methods (Access = private)

    function SetMatrixData(this,diagData)
      % Assign apply-inverse data.
      this.diagData_ = diagData;
    end

    function SetDiagonalInvData(this,invDiagData)
     % Assign apply-inverse data.
     this.invDiagData_ = invDiagData;
    end

    function SetCollection(this,Collection)
     this.Collection_ = Collection;
    end
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
end %classdef OperatorView
