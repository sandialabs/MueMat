classdef CGEminSolver < VerboseObject

  properties (Access=private)
    nIts_ = 2
    SatisfyConstraintsFunc_ = @EminPFactory.SatisfyConstraints
  end

  methods

    function this = CGEminSolver(nIts)
      if nargin == 1 && isa(nIts, class(this)), this.Copy_(nIts,[]); return; end

      if varexist('nIts'), this.nIts_ = nIts; end
    end

    function SetNumIterations(this, nIts)
      this.nIts_ = nIts;
    end

    function SetNeeds(this,FineLevel, CoarseLevel)
    % no needs
    end

    function P = Iterate(this, A, Pattern, B, P0, FineLevel, CoarseLevel, prolongation_mode)
      P = this.Iterate_(A.GetMatrixData(), Pattern, B, P0.GetMatrixData(), CoarseLevel);
      P = Operator(P, P0.GetRowMap(), P0.GetColMap(), P0.GetApply());
    end

    function P = Iterate_(this, A, Pattern, B, P0, CoarseLevel)
      %
      % This is CG applied to the problem
      %
      %  min  (1/2) sum (p_i)^T A p_i
      %   P          i
      %
      %  subject to  P(i,j) = 0  if  (i,j) is not in SparsityPattern
      %      and     P cnull =  fnull.
      %

      % Initialize
      %
      %blecko = (A ~= 0);
      %A(options.difrows,:) = options.S(options.difrows,:) .* blecko(options.difrows,:);
      %A(:,options.difrows) = options.S(:,options.difrows) .* blecko(:,options.difrows);
      %A = options.S;
      n       = size(A,1);
      ddd     = spdiags(1 ./ diag(A),0,n,n);
      %ddd(options.difrows,options.difrows) = sparse(length(options.difrows),length(options.difrows));
      %scaleme = ones(n,1); scaleme(options.difrows) =  -0.00004; ddd = spdiags(scaleme,0,n,n)*ddd;
      %ddd(options.difrows,options.difrows) = sparse(length(options.difrows),length(options.difrows));

      P=P0;

      % CG iterations

      rk = this.SatisfyConstraintsFunc_(-A*P, [], CoarseLevel, Pattern, B);

      %solve iteratively for the correction
      % Note: we need to have P, rk, pk, and ap

      fprintf('initial emin residual = %4.3e\n', full(max(max(abs(rk)))));
      i = 1;
      while ( (i<=this.nIts_) && (max(max(abs(rk)))>1e-12))
        zk = ddd * rk;    % precondition by inverse diagonal
        newsum = full(sum(diag(rk'*zk)));
        if (i == 1)
          pk = zk;
        else
          betak = newsum/oldsum;
          pk = zk + pk * betak;
        end
        oldsum = newsum;

        % multiply A pk

        ap = this.SatisfyConstraintsFunc_(A*pk, [], CoarseLevel, Pattern, B);

        alphak = newsum/sum(diag(pk'*ap));

        P   = P  + alphak*pk;
        rk  = rk - alphak*ap;
        %fprintf('%d: r = %20.13e\n', i, full(max(max(abs(rk)))));
        i = i+1;
      end
      fprintf('final emin residual   = %4.3e\n', full(max(max(abs(rk)))));

    end

  end
end