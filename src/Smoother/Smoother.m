classdef Smoother < SmootherBase & SmootherPrototype
  % Jacobi, Gauss-Seidel smoother class (and more)
  %
  % Fairly general class for defining smoothers and block smoothers such as
  % Gauss-Seidel, block Gauss-Seidel, domain decomposition and Jacobi
  %
  % The basic idea is that one uses BDiag and/or Collection to define a set of
  % small square block matrices. These matrices are inverted with residual
  % calculations appearing between inversions. In this way, several different
  % general methods are incorporated.
  %
  % Conceptually, one iteration of the smoother looks like
  %
  %     for m=1: Nblks
  %        r^(m)   = b - A u^(m)
  %        u^(m+1) = u^(m) + omega*Pinv^(m)* r^(m)
  %     end
  %
  % where Pinv^(m) is a matrix of all zeros except for a square diagonal block
  % that corresponds to the inverse of a dense square diagonal block.
  % Note: Pinv's nonzero rows/cols might be scattered (i.e. Pinv has one dense
  %       block matrix corresponding to the diagonal after a suitable permutation).
  % Note: Since Pinv^(m) is mostly zeros, it is not necessary to compute
  %       all elements of r^(m). That is, we only compute r(row_inds,:)
  %       corresponding to nonzero locations in Pinv^(m). This uses the
  %       selective multiply option.
  %
  % Here are examples:
  %
  %    point Gauss-Seidel           |    Nblks = size(A,1)
  %                                 |    Pinv^(m)(m,m) = 1/a(m,m)
  %    -------------------------------------------------------------------
  %    blk Gauss-Seidel             |    Nblks = size(A,1)/2
  %    with 2x2 blocks              |    Pinv^(m)(ii,ii) = inv(a(ii,ii))
  %                                 |    ii = (2*m-1:2*m).
  %    -------------------------------------------------------------------
  %    overlapping or nonoverlapping|    Nblks = # domains
  %    multiplicative domain        |    Pinv^(m)(ii,ii) = inv(a(ii,ii))
  %    decomposition                |    ii = unknowns residing in mth domain
  %    -------------------------------------------------------------------
  %
  % Additive versions (Jacobi, block Jacobi, additive domain decomposition)
  % are available and correspond to replacing r^(m)   = b - A u^(m) with a
  % single residual calculation before the loop above.
  %
  % Implementation of Domain Decomposition
  % =======================================
  % For domain decomposition, Pinv is implemented as a list of Nblks matrices.
  % Components (corresponding to the mth ii) from r are gathered before the solve
  % associated with Pinv^(m) and results are scattered back to elements in u^(m)
  %
  %
  % See also SmoothingTest.m, SmootherFactory.m
  %
  %
  properties (Access = private)
    % Parameters of the smoother
    nIts_  = 2                 % number of iteration (integer)
    Omega_ = 1.0               % damping factor
    ForwardSweep_ = true
    BackwardSweep_
    JacobiStyle_               % todo:remove
    diagonalView_ = 'current'  % diagonal view label ([] == current view)

    % Data after Setup phase
    Amat_           % matrix of the level
    Atrans_         % transpose of the Matlab matrix, for fast access to subset of rows
    GSFastD_        % internal data for fast version of algorithms
    GSFastL_        % internal data for fast version of algorithms
    GSFastU_        % internal data for fast version of algorithms
    DL_Lfactor_
    DL_Ufactor_
    DU_Lfactor_
    DU_Ufactor_
  end

  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [this] = Smoother(string, nIts, Omega, diagonalView)
      %SMOOTHER Constructor
      %
      %   SYNTAX   Smoother(string, nIts, Omega, diagonalView);
      %
      %     string         - smoother name (string)
      %     nIts           - number of iteration (integer, optional, default=2)
      %     Omega          - damping factor (double, optional, default=1)
      %     diagonalView   - Block or point smoother (optional)
      %
      %   EXAMPLE
      %
      %     Smoo = Smoother('BlkGaussSeidel', 2, 1);

      % Copy Constructor
      if nargin == 1 && isa(string, class(this)), this.Copy_(string,[]); return; end
      %

      if varexist('nIts'),  this.SetNIts(nIts); end
      if varexist('Omega'), this.SetOmega(Omega); end
      if varexist('diagonalView'), this.SetDiagonalView(diagonalView); end;

      % TODO : Smoother(string,...) vs JacobiStyle= and BackwardSweep=
      if strcmp(string,'Jacobi')
        this.JacobiStyle_ = true;
        this.SetBackwardSweep(false);
      elseif strcmp(string,'GaussSeidel')
        this.JacobiStyle_ = false;
        this.SetBackwardSweep(true);
      elseif strcmp(string,'SymmGaussSeidel')
        this.JacobiStyle_ = false;
        this.SetForwardSweep(true);
        this.SetBackwardSweep(true);
      else
        fprintf('Smoother: Unknown smoother %s\n',string)
        keyboard;
      end

      this.SetType(string);
    end% Smoother()

    function SetNeeds(this, Level)
       % Obtain any cross factory specifications
    end

    function SetNIts(this, nIts)
      %SETNITS Set the number of iterations performed each time Apply() is invoked.
      %
      %   SYNTAX   obj.SetNIts(nIts);
      %
      %     nIts - number of iteration (integer)
      this.nIts_ = nIts;
    end
    function SetOmega(this, Omega)
      %SETOMEGA Set the damping parameter used each time Apply() is invoked on constructed smoother objects.
      %
      %   SYNTAX   obj.SetOmega(Omega);
      %
      %     Omega - damping factor (double)
      this.Omega_ = Omega;
    end
    function SetForwardSweep(this, ForwardSweep)
      %SETFORWARDSWEEP Indicates whether forwards sweeps are used each time Apply() is invoked on constructed smoother objects.
      %
      %   SYNTAX   obj.SetForwardSweep(ForwardSweep);
      %
      %     ForwardSweep - forward sweep (boolean)
      this.ForwardSweep_ = ForwardSweep;
    end
    function SetBackwardSweep(this, BackwardSweep)
      %SETBACKWARDSWEEP Indicates whether backward sweeps are used each time Apply() is invoked on constructed smoother objects.
      %
      %   SYNTAX   obj.SetBackwardSweep(BackwardSweep);
      %
      %     BackwardSweep - backward sweep (boolean)
      this.BackwardSweep_ = BackwardSweep;
    end
    function SetJacobiStyle(this, JacobiStyle)
      %SETJACOBISTYLE
      %
      %   SYNTAX   obj.SetJacobiStyle(JacobiStyle);
      %
      %     JacobiStyle - (boolean)
      this.JacobiStyle_  = JacobiStyle;
    end

    function SetDiagonalView(this, diagonalView)
      %SETBLOCKINGSCHEME Indicates a strategy for how blocks are define each time Apply() is invoked on constructed smoother objects.
      %
      %   SYNTAX   obj.SetDiagonalView(diagonalView);
      %

      % TODO: if default==point or default==block or current ...
      if ~strcmp(this.diagonalView_, diagonalView)
        this.Amat_ = []; this.GSFastD_ = []; this.GSFastL_ = []; this.GSFastU_ = [];
        this.SetIsSetup(false);
      end

      this.diagonalView_ = diagonalView;
    end

    function [nIts] = GetNIts(this)
      %GETNITS Get the number of iterations
      %
      %   SYNTAX   nIts = obj.GetNIts();
      %
      %     nIts - number of iteration (integer)
      nIts  = this.nIts_;
    end
    function [Amat] = GetA(this)
      %GETNITS Get the matrix A associated with smoother
      %
      %   SYNTAX   Amat = obj.Get('A');
      %
      %     Amat - Matrix associated with smoother
      Amat  = this.Amat_;
    end
    function [Omega] = GetOmega(this)
      %GETOMEGA Get the damping parameter
      %
      %   SYNTAX   Omega = obj.GetOmega();
      %
      %     Omega - damping factor (double)
      Omega  = this.Omega_;
    end
    function [ForwardSweep] = GetForwardSweep(this)
      %GETFORWARDSWEEP Get if forward sweeps are used
      %
      %   SYNTAX   ForwardSweep = obj.GetForwardSweep();
      %
      %     ForwardSweep - forward sweep (boolean)
      ForwardSweep  = this.ForwardSweep_;
    end
    function [BackwardSweep] = GetBackwardSweep(this)
      %GETBACKWARDSWEEP Get if backward sweeps are used
      %
      %   SYNTAX   BackwardSweep = obj.GetBackwardSweep();
      %
      %     BackwardSweep - backward sweep (boolean)
      BackwardSweep  = this.BackwardSweep_;
    end
    function [JacobiStyle] = GetJacobiStyle(this)
      %GETJACOBISTYLE
      %
      %   SYNTAX   JacobiStyle = obj.GetJacobiStyle();
      %
      %     JacobiStyle - (boolean)
      JacobiStyle  = this.JacobiStyle_;
    end
    function [diagonalView] = GetDiagonalView(this)
      %GETBLOCKINGSCHEME Get the block strategy
      %
      %   SYNTAX   diagonalView = obj.GetDiagonalView();
      diagonalView  = this.diagonalView_;
    end

    function CopyParameters(this, src)
      %COPYPARAMETERS Copy the parameters of another smoother prototype.
      % See also: SmootherPrototype.CopyParameters
      %
      %   SYNTAX   obj.CopyParameters(src);
      %
      %     src - Object (SmootherPrototype)
      this.SetNIts(src.GetNIts());
      this.SetOmega(src.GetOmega());
      this.SetForwardSweep(src.GetForwardSweep());
      this.SetBackwardSweep(src.GetBackwardSweep());
      this.SetJacobiStyle(src.GetJacobiStyle());
      this.SetDiagonalView(src.GetDiagonalView());
    end % function

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setup phase
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function Setup(this, Level, Specs)
      %SETUP Run the setup phase of the smoother.
      % See also: SmootherPrototype.Setup, Smoother.SetupDiagonal, Smoother.Setup
      %
      %   SYNTAX   Amat = obj.Setup(Level, Specs);
      %
      %     Level - level of the MG hierachy (Level)
      %     Specs - specifications (CrossFactory)

      % This Smoother object might have an alternative definition
      % of the block diagonal or an alternative definition of 'block'
      % as well as a different damping parameter and number of iterations
      % associated with Smoother.Apply() depending on factory specifications.

      if this.isSetup(), return; end

      %TODO        if varexist('Specs'), this.TempOutputLevel(Specs); end;

      this.Amat_ = Level.Get('A');
      this.Atrans_ = this.Amat_.GetMatrixData()';

      this.SetupDiagonal(); % 'view' cuisine
      this.SetupCache();    % initialize GSFastL_ U_ D_ for fast version of the algo

      this.SetIsSetup(true);

      %TODO        if varexist('Specs'), this.RestoreOutputLevel(); end;
    end

    function SetupDiagonal(this)
      %OVERLOADSETUP Submethod of Setup()
      % See also: Smoother.Setup
      %
      %   SYNTAX   Amat = obj.SetupDiagonal();
      %

      Amat = this.Amat_;

      % Builds an appropriate Smoother object
      %
      % This Smoother object might have an an alternative definition
      % of the block diagonal or an alternative definition of 'block'
      % as well as a different damping parameter and number of iterations
      % associated with Smoother.Apply() depending on factory specifications.
      % TODO remove BDiag

      if strcmp(this.diagonalView_,'Random NonOverlapping') || strcmp(this.diagonalView_,'Random Overlapping'),
        Smoother.AddDiagTestView(Amat, this.diagonalView_);
      end

      % Setup Diagonal
      BDiag = Amat.GetDiagonal([], this.diagonalView_); % Diag is build only if this have not be done yet
      if isempty(BDiag.GetApplyInverse())
        FactorBlkDiag(BDiag);
      end

    end %SetupDiagonal()

    function SetupCache(this)
      %SUBSETUP Submethod of Setup()
      % See also: Smoother.Setup
      %
      %   SYNTAX   obj.SubSetup();
      %

      Amat = this.Amat_;

      % special code to make the case of contiguous blocks faster
      %
      BDiag = Amat.GetDiagonal([], this.diagonalView_);
      if ~this.JacobiStyle_ && isempty(BDiag.GetCollection()) && BDiag.GetRowMap().HasConstBlkSize(),
        nPDE = BDiag.GetRowMap().ConstBlkSize();
        nn   =  BDiag.GetRowMap().NDOFs();
        [aaa,bbb,ccc] = find(Amat.GetMatrixData());
        DiagInds  = find( floor(aaa/nPDE - .00001) == floor(bbb/nPDE - .00001));
        UpperInds = find( floor(aaa/nPDE - .00001)  < floor(bbb/nPDE - .00001));
        LowerInds = find( floor(aaa/nPDE - .00001)  > floor(bbb/nPDE - .00001));
        this.GSFastL_=sparse(aaa(LowerInds),bbb(LowerInds),ccc(LowerInds),nn,nn);
        this.GSFastU_=sparse(aaa(UpperInds),bbb(UpperInds),ccc(UpperInds),nn,nn);
        this.GSFastD_=sparse(aaa( DiagInds),bbb( DiagInds),ccc( DiagInds),nn,nn);
        [this.DU_Lfactor_,this.DU_Ufactor_] = lu(this.GSFastD_ + this.Omega_*this.GSFastU_);
        [this.DL_Lfactor_,this.DL_Ufactor_] = lu(this.GSFastD_ + this.Omega_*this.GSFastL_);
      end

    end % SetupCache()


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
      %     SolStatus - when InitGuessStatus==ALLZEROS, iu is
      %     assumed to be all zeros and initial matvec operation
      %     could be avoid (optional,default=NOTALLZEROS)
      %
      % Apply block Gauss-Seidel iterations according to the information in
      % Smoother.Setup() should be invoked prior to Smoother.Apply() to
      % set Smoother. This includes setting the total number of iterations
      % and also allows for specify an alternative block diagonal matrix
      % and defining alternative blocking schemes.
      %
      % Note: In the most general case, the inverted block diagonal matrix
      %       may be stored differently than Amat. In particular, if
      %       BlkDiag correspond to some arbitrary blocking (grouping) of
      %       matrix rows in Amat, these rows have been gathered together
      %       and are stored consecutively in BlkDiag. Here is an example
      %       to clarify things.
      %             MatrixData = [a1; a2; a3; a4]; % ai is the ith row of A.
      %             newMap = Map(4, 1);
      %             A.MatrixData     = MatrixData;
      %             A.RowMap         = Map;
      %             A.ColMap         = Map;
      %             A.Apply              = @MatlabApply;
      %
      %       We want to do block smoothing using a diagonal matrix with
      %       two blocks defined by grouping together a1,a3,a4 and grouping
      %       together a2,a4. To do this, the user sets Collection via
      %             Collection.NSubsets = 2;
      %             Collection.Subsets(1) = CreateDOFSubset(A.RowMap,'Scattered',-1,-1,[1 3 4] );
      %             Collection.Subsets(2) = CreateDOFSubset(A.RowMap,'Scattered',-1,-1,[2 4] );
      %       and then invokes
      %             [Smoo,A] = Smoo.Setup(A,1,1.,[],Collection);
      %
      %       This builds a block diagonal matrix in Smoother.BlkDiag that has a
      %       total of blocks. The first is 3x3 and the second is 2x2 which are
      %       stored contiguously. When we do Gauss-Seidel, the first matrix-vector
      %       product works only on the subset corresponding to a1,a3,a4 which
      %       is Collection.Subsets(1). The result is placed in a vector of length 3.
      %       When we apply Dinv, we use a Subset corresponding to the first block
      %       which is rows d1,d2,d3. The second matrix vector product within
      %       a GS iteration corresponds to a subset a2,a4. The result of this
      %       is applied to the 2nd block (d4,d5) of Dinv.
      if ~this.isSetup(), error('apply'); end

      if   isa(iu,'MultiVector'), u = iu.GetVectors();
      else                        u = iu;                   end
      if   isa(irhs,'MultiVector'), rhs = irhs.GetVectors();
      else                          rhs = irhs;             end

      Amat = this.Amat_;
      previousView = Amat.SwitchToView(this.diagonalView_);

      mue_include

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % fast version which can only be used for
      % constant block size and contiguously defined blocks (with omega = 1).
      %
      if ~isempty(this.GSFastD_)
        omega = this.Omega_;
        OneMinusomega = 1 - this.Omega_;
        for k=1:this.nIts_
          if this.ForwardSweep_,
            if SolStatus == ALLZEROS,
              u = omega*(this.DL_Ufactor_ \ (this.DL_Lfactor_\rhs) );
            else u = this.DL_Ufactor_ \ (this.DL_Lfactor_\ ...
                (omega*(rhs - this.GSFastU_*u) + OneMinusomega*this.GSFastD_*u));
            end
            SolStatus = NOTALLZEROS;
          end
          if this.BackwardSweep_,
            if SolStatus == ALLZEROS,
              u = omega*(this.DU_Ufactor_ \ (this.DU_Lfactor_\rhs) );
            else u = this.DU_Ufactor_ \ (this.DU_Lfactor_\ ...
                (omega*(rhs - this.GSFastL_*u) + OneMinusomega*this.GSFastD_*u));
            end
            SolStatus = NOTALLZEROS;
          end
          %fprintf('%20.13e\n',norm(u));
        end
      else
        AmatRowMap = this.Amat_.GetRowMap();
        BlkDiag = this.Amat_.GetDiagonal();
        BlkDiagRowMap = BlkDiag.GetRowMap();
        Collection = BlkDiag.GetCollection();
        AmatApplyFcn = this.Amat_.GetApply();
        for k=1:this.nIts_
          if this.ForwardSweep_,
            Nsteps = BlkDiagRowMap.NNodes();

            % if JacobiStyle compute the residual at the start
            if this.JacobiStyle_,
              res = rhs;
              if SolStatus ~= ALLZEROS,
                res = res - AmatApplyFcn(this.Amat_, u, ...
                  CreateDOFSubset(AmatRowMap,'All'));

              end
              % fast version of Jacobi
              if isempty(Collection),
                Nsteps = 1;
                BDiagSubset=CreateDOFSubset(BlkDiagRowMap,'Contiguous',1,...
                  BlkDiagRowMap.NNodes(),[]);
              end
            end
            for i=1:Nsteps
              % fast version of Jacobi
              if ~this.JacobiStyle_ || (this.JacobiStyle_ && ~isempty(Collection))
                %TODO: simplified test :-)
                BDiagSubset=CreateDOFSubset(BlkDiagRowMap,'Contiguous',i,i,[]);
              end
              if ~isempty(Collection), % simplify also this one
                AmatSubset= Collection.Subsets(i);
              else  AmatSubset = BDiagSubset;
              end
              indices = Subset2DOF(AmatSubset,AmatRowMap);

              % compute sub residual. If this is a Jacobi-style algorithm, then
              % we just need to take a subset of the already computed residual.
              % If not Jacobi-style, we must perform a local A*x.
              if this.JacobiStyle_,  rsub = res(indices);
              else
                if SolStatus == ALLZEROS, rsub = rhs(indices);
                else                      rsub=rhs(indices) - (u'*this.Atrans_(:,indices))';
                end
              end

              applyInvFcn = BlkDiag.GetApplyInverse();
              mue_include; % #define LOCAL
              u(indices)=u(indices)+this.Omega_*applyInvFcn(BlkDiag, rsub, BDiagSubset, LOCAL);
              SolStatus = NOTALLZEROS;
              % code to check local residual after individual Dinv's
              %t = norm(rhs(indices)-this.Amat.Apply(this.Amat,u,AmatSubset));
              %if ( t > 1e-6), keyboard; end
            end
          end
          if this.BackwardSweep_,
            if this.JacobiStyle_,
              res = rhs;
              if SolStatus ~= ALLZEROS,
                res = res - this.Amat_.Apply( u, ...
                  CreateDOFSubset(AmatRowMap,'All'));
              end
            end
            for i=BlkDiagRowMap.NNodes():-1:1
              BDiagSubset=CreateDOFSubset(BlkDiagRowMap,'Contiguous',i,i,[]);
              if  ~isempty(Collection),
                AmatSubset= Collection.Subsets(i);
              else AmatSubset = BDiagSubset;
              end
              indices = Subset2DOF(AmatSubset,AmatRowMap);

              % compute sub residual. If this is a Jacobi-style algorithm, then
              % we just need to take a subset of the already computed residual.
              % If not Jacobi-style, we must perform a local A*x.
              if this.JacobiStyle_,  rsub = res(indices);
              else
                if SolStatus == ALLZEROS, rsub = rhs(indices);
                else                      rsub=rhs(indices) - (u'*this.Atrans_(:,indices))';
                end

                applyInvFcn = BlkDiag.GetApplyInverse();
                LOCAL=2; % input type for the apply function
                u(indices)=u(indices)+this.Omega_*applyInvFcn(BlkDiag, rsub, BDiagSubset, LOCAL);
                SolStatus = NOTALLZEROS;
                % code to check local residual after individual Dinv's
                %t = norm(rhs(indices)-this.Amat.Apply(this.Amat,u,AmatSubset));
                %if ( t > 1e-6), keyboard; end
              end % for i=BlkDiagRowMap.NNodes():-1:1
            end
          end % if  this.sym_,
        end %for k=1:this.nIts_
      end % if strcmp(this.type_,'GSFast'),

      Amat.SwitchToView(previousView);

      if   isa(iu,'MultiVector'), iu.SetVectors(u);
      else                        iu = u;                   end
    end% function

    function Print(this,prefix)
      %PRINT Print smoother information
      %
      %   SYNTAX   obj.Print()
      %
      %     prefix  - optional string that is prepended to each print line

      if ~varexist('prefix'), prefix = ''; end
      sweepMode = 'unknown';
      if this.ForwardSweep_  && this.BackwardSweep_
        sweepMode = 'symmetric';
      elseif this.ForwardSweep_
        sweepMode = 'forward';
      elseif this.BackwardSweep_
        sweepMode = 'backward';
      end
      fprintf('%s%s: sweeps=%d, omega=%g, mode=''%s''\n',...
              prefix,this.GetType(),this.nIts_,this.Omega_,sweepMode);

    end %Print()

  end % methods


  % Auxiliary methods for SETUP
  methods (Access = private, Static = true)

    function AddDiagTestView(Amat, diagonalView)
      %GETBLKDIAGSPECS Sets up the diagonal, block diagonal, and the collection according
      % to specifications.
      %
      %   SYNTAX   [Amat, Collection] = obj.GetBlkDiagSpecs(Amat, BlockingScheme, Collection);
      %
      %     Amat           - matrix of the current level (Operator)
      %     BlockingScheme - Block or point smoother (string=BlkDiagonal or PointDiagonal)
      %     Collection     - how to group blocks together for the smoothing process (Collection)
      %
      %
      % Note: a collection is a set of indices (rows and columns) which is to
      %       be treated as an entire block when performing something like
      %       block Jacobi or Gauss-Seidel. This lets us do some form of
      %       additive or multiplicative Schwarz method.
      %

      % Note: the view (point or block) must be set
      % before calling Diagonal(Amat,Collection).
      % Here, this is done by the calling function (OverloadSetup).
      % It's important to set the view because method
      % ExtractNonContigBlkDiag use it
      %

      if strcmp(diagonalView,'Random NonOverlapping'),
        Collection = MakeUpRandomBlks(Amat,'NonOverlapping');
      elseif strcmp(diagonalView,'Random Overlapping'),
        Collection = MakeUpRandomBlks(Amat,'Overlapping');
      else
        fprintf('Smoother: Unknown diagonal view\n');
        keyboard;
      end

      BDiag = Diagonal(Amat,Collection);         % build the special diag
      Amat.CreateView(diagonalView,Amat.GetRowMap(),Amat.GetColMap(),Amat.GetApply());
      Amat.SetDiagonal(BDiag, diagonalView);
      % Note:    DBDiag = Diagonal(Amat,Collection)
      %       +  Amat.SetDiagonal(BDiag, diagonalView)
      %       == FineA.GetDiagonal(Collection, 'ovblock');

    end % GetBlkDiagSpecs()

  end % method

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

% TODO: cleanup comments
