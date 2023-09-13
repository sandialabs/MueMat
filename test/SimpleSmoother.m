classdef SimpleSmoother < SmootherBase & SmootherPrototype
  % SimpleSmoother class
  % This smoother applies one smoother followed by a second smoother 
  %
  
  properties (Access = public)
    % Parameters of the smoother
    F_
    Bt_
    B_
    C_
    SmootherF1_           
    SmootherF2_           
    SmootherS_            
    nv_
    np_
    alpha_
    UseDinv_
  end
  
  methods (Static = true)
    function [OutVec] = SchurApply(Op,Vec, Subset)
      % 
      %   SYNTAX   OutVec = obj.SchurApply(Op, Vec, Subset);
      % 
      mue_include
      if nargout() ~= 1, error('SchurApply method returns vector'); end

      widget = Op.GetMatrixData();
      t1=widget.SmootherF1_.Apply(zeros(widget.nv_,1), widget.Bt_*Vec, NOTALLZEROS);
      OutVec = widget.C_*Vec - widget.B_*t1;
    end %SchurApply
  end
  methods
    function [Smoo] = GetSmootherF1(this)
        Smoo = this.SmootherF1_;
    end;
    function [Smoo] = GetSmootherF2(this)
        Smoo = this.SmootherF2_;
    end;
    function [Smoo] = GetSmootherS(this)
        Smoo = this.SmootherS_;
    end;
    function [alpha] = GetAlpha(this)
        alpha = this.alpha_;
    end;
    function [UseDinv] = GetUseDinv(this)
        UseDinv = this.UseDinv_;
    end;
    function [UseDinv] = SetUseDinv(this,UseDinv)
        this.UseDinv_ = UseDinv;
    end;
    function [this] = SimpleSmoother(SmootherF1,SmootherF2,SmootherS,alpha)
      
      % Copy constructor
      if nargin == 1 && isa(SmootherF1, class(this)), this.Copy_(SmootherF1,[]); return; end

      if nargin ~= 4, fprintf('SimpleSmoother needs 4 arguments\n'); keyboard; end;
      this.SmootherF1_ = SmootherF1;
      this.SmootherF2_ = SmootherF2;
      this.SmootherS_  = SmootherS;
      this.alpha_      = alpha;
      this.UseDinv_    = false;
      this.SetType('SimpleSmoother');
    end % function

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameters / Config
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function CopyParameters(this, src)
      %COPYPARAMETERS Copy the parameters of another smoother prototype.
      % See also: SmootherPrototype.CopyParameters, Hybrid2x2Smoother.CopyParameters
      % 
      %   SYNTAX   obj.CopyParameters(src);
      % 
      %     src - Object (SmootherPrototype of same type)

      srcSmootherF1 = src.GetSmootherF1();
      if ~isempty(this.SmootherF1_) && ~isempty(srcSmootherF1) && ...
            (strcmp(this.SmootherF1_.GetType(), srcSmootherF1.GetType()))
        this.SmootherF1_.CopyParameters(srcSmootherF1);
      else
        this.SmootherF1_ = srcSmootherF1.Copy(); % Copy because we don't want to modify 'src'
      end
      srcSmootherF2 = src.GetSmootherF2();
      if ~isempty(this.SmootherF2_) && ~isempty(srcSmootherF2) && ...
            (strcmp(this.SmootherF2_.GetType(), srcSmootherF2.GetType()))
        this.SmootherF2_.CopyParameters(srcSmootherF2);
      else
        this.SmootherF2_ = srcSmootherF2.Copy(); % Copy because we don't want to modify 'src'
      end
      srcSmootherS = src.GetSmootherS();
      if ~isempty(this.SmootherS_) && ~isempty(srcSmootherS) && ...
            (strcmp(this.SmootherS_.GetType(), srcSmootherS.GetType()))
        this.SmootherS_.CopyParameters(srcSmootherS);
      else
        this.SmootherS_ = srcSmootherS.Copy(); % Copy because we don't want to modify 'src'
      end
      this.alpha_   = src.GetAlpha();
      this.UseDinv_ = src.GetUseDinv();
      
      this.SetIsSetup(this.SmootherF1_.isSetup() && ...
                 this.SmootherF2_.isSetup() && this.SmootherS_.isSetup() );
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

      OrigA    = Level.Get('A');
      Adata    = OrigA.GetMatrixData(); 
      nv       = Level.Get('N11');
      np       = Level.Get('N22');
      this.nv_ = nv;  
      this.np_ = np;
      this.F_  = Adata(1:nv,1:nv);
      this.B_  = Adata(nv+1:nv+np,1:nv);
      this.Bt_ = Adata(1:nv,nv+1:nv+np);
      this.C_  = Adata(nv+1:nv+np,nv+1:nv+np);

      Level.Set('A', Operator(this.F_,Map(nv,1),Map(nv,1),@MatlabApply));
      if ~this.SmootherF1_.isSetup()
        this.SmootherF1_.Setup(Level);
      end
      if ~this.SmootherF2_.isSetup()
        this.SmootherF2_.Setup(Level);
      end

      % Build an implicit approximate Schur complement operator with a matvec
      % function corresponding to 
      %
      %            C*v - B * SmootherF1_.Apply(Bt*v)
      %
      % and an artificial diagonal (stored in the diagonal view 'mydiag')
      % corresponding to 
      %
      %          diag(C - B* inv( sum(abs(F)) )  * Bt)
      %
     
      if this.UseDinv_,
         %ddd = .3*sum(abs(this.F_)')';
         ddd = sum(abs(this.F_)')';
         ddd = this.C_ - this.B_*spdiags(1./ddd,0,nv,nv)*this.Bt_;
         Sapprox = Operator(ddd, 1, 1);
         Sapprox.CreateView('mydiag',Sapprox.GetRowMap(),Sapprox.GetColMap(),...
                      Sapprox.GetApply());
      else
         Sapprox = Operator(sparse(np,np),Map(np,1), Map(np,1), ...
                            @SimpleSmoother.SchurApply,'pfo');

         widget.SmootherF1_= this.SmootherF1_;
         widget.nv_        = this.nv_;
         widget.Bt_        = this.Bt_;
         widget.C_         = this.C_;
         widget.B_         = this.B_;
         Sapprox.SetMatrixData(widget);
         ddd = sum(abs(this.F_)')';
         ddd = diag(this.C_ - this.B_*spdiags(1./ddd,0,nv,nv)*this.Bt_);
         ddd = Diagonal(Operator(spdiags(ddd,0,np,np), Sapprox.GetRowMap(),...
                        Sapprox.GetColMap(), @SimpleSmoother.SchurApply,'poo'));
         Sapprox.CreateView('mydiag',Sapprox.GetRowMap(),Sapprox.GetColMap(),...
                      @SimpleSmoother.SchurApply);
         Sapprox.SetDiagonal(ddd,'mydiag');
      end
   

      Level.Set('A', Sapprox);
      if ~this.SmootherS_.isSetup()
        this.SmootherS_.Setup(Level);
      end
      Level.Set('A', OrigA);
      
      this.SetIsSetup(this.SmootherF1_.isSetup() && this.SmootherF2_.isSetup() && this.SmootherS_.isSetup());
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
      mue_include;      

      if ~this.isSetup(), error('apply'); end
      
      if   isa(iu,'MultiVector'), u = iu.GetVectors();
      else                        u = iu;                   end
      if   isa(irhs,'MultiVector'), rhs = irhs.GetVectors();
      else                          rhs = irhs;             end
      
      nv = this.nv_;
      np = this.np_;

      if SolStatus == NOTALLZEROS,
         rhs = rhs - [this.F_*u(1:nv)+this.Bt_*u(nv+1:nv+np);...
                                         this.B_*u(1:nv)+ this.C_*u(nv+1:nv+np)];
      end;

      t1 = this.SmootherF1_.Apply(zeros(nv,1), rhs(1:nv), SolStatus);
      t2 = rhs(nv+1:nv+np)  - this.B_*t1;
      t3 = this.SmootherS_.Apply(zeros(np,1), t2, SolStatus);
      u(nv+1:nv+np) = u(nv+1:nv+np) + this.alpha_*t3;
      u(1:nv) = u(1:nv) + t1 - spdiags(1 ./sum(abs(this.F_)')',0,nv,nv)*(this.Bt_*t3);
%      u(1:nv) = u(1:nv) + t1 - this.SmootherF1_.Apply(zeros(nv,1),this.Bt_*t3,ALLZEROS);
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
      fprintf('subsmoother for F1\n');
      this.SmootherF1_.Print(prefix);
      fprintf('subsmoother for F2\n');
      this.SmootherF2_.Print(prefix);
      fprintf('subsmoother for S\n');
      this.SmootherS_.Print(prefix);
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
