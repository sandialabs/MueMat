classdef Hybrid2x2Smoother < SmootherBase & SmootherPrototype
% Hybrid smoother: merges two basic smoothers for use on a 2x2 system
%
%       ( A_rr    A_rx  )  ( u_r )  ( u_r )     ( b_r )
%       (               )  (     )  (     )  =  (     )
%       ( A_xr    A_xx  )  ( u_x )  ( u_x )     ( b_x )
%
%  via
%        a) StartTwo its of Smoother2(Axx, u_x, b_x - A_xr u_r)
%        b) MainIts of
%                - MiddleOne its of Smoother1(Arr, u_r, b_r - A_rx u_x)
%                - MiddleTwo its of Smoother2(Axx, u_x, b_x - A_xr u_r)
%        c) EndTwo its of Smoother2(Axx, u_x, b_x - A_xr u_r)
%

  properties (Access = private)
    % Parameters of the smoother
    nMainIts_   % number of iterations (integer)
    MiddleOne_  % number of iterations (integer)
    StartTwo_   % number of iterations (integer)
    MiddleTwo_  % number of iterations (integer)
    EndTwo_     % number of iterations (integer)

    % Parameters (prototypes) and data after setup
    SmootherOne_ % smoother of the first system block (A_rr)
    SmootherTwo_ % smoother of the second system block (A_xx)

    fakeLevelOne_ = [] % Level data for Smoother One (it's just a wrapper around the matrix AOne)
    fakeLevelTwo_ = [] % Level data for Smoother Two (it's just a wrapper around the matrix ATwo)

    % Data after Setup phase
    Amat_         % matrix of the level
  end

  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [this] = Hybrid2x2Smoother(nMainIts, MiddleOne, ...
                                        StartTwo, MiddleTwo, EndTwo, ...
                                        SmootherOne, SmootherTwo)
      %HYBRID2X2SMOOTHER Constructor
      %
      %   SYNTAX   Hybrid2x2Smoother(nMainIts, MiddleOne, StartTwo, MiddleTwo, EndTwo, SmootherOne, SmootherTwo);
      %
      %     nMainIts    - number of iterations (integer)
      %     MiddleOne   - number of iterations (integer)
      %     StartTwo    - number of iterations (integer)
      %     MiddleTwo   - number of iterations (integer)
      %     EndTwo      - number of iterations (integer)
      %     SmootherOne - smoother of the first system block (A_rr)
      %     SmootherTwo - smoother of the second system block (A_xx)

      % Copy constructor
      if nargin == 1 && isa(nMainIts, class(this)), this.Copy_(nMainIts,[]); return; end
      %

      %if nargout() < 1, error('nargout'); end
      % TODO : manage optional parameters

      % TODO Do we need to check if Amat has a diagonal here ?
      this.SetIts(nMainIts, MiddleOne, StartTwo, MiddleTwo, EndTwo);
      this.SetSmoother(SmootherOne, SmootherTwo);

      this.SetType('hybrid2x2');
    end% function

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameters / Config
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TODO: exhaustive Set/Get
    function SetIts(this, nMainIts, MiddleOne, StartTwo, MiddleTwo, EndTwo)
      %SETITS
      %
      %   SYNTAX   obj.SetIts(nMainIts, MiddleOne, StartTwo, MiddleTwo, EndTwo);
      %
      %     nMainIts  - number of iterations (integer)
      %     MiddleOne - number of iterations (integer)
      %     StartTwo  - number of iterations (integer)
      %     MiddleTwo - number of iterations (integer)
      %     EndTwo    - number of iterations (integer)
      this.nMainIts_    = nMainIts;
      this.MiddleOne_   = MiddleOne;
      this.StartTwo_    = StartTwo;
      this.MiddleTwo_   = MiddleTwo;
      this.EndTwo_      = EndTwo;
    end

    function [nMainIts, MiddleOne, StartTwo, MiddleTwo, EndTwo] = GetIts(this)
      %GETITS
      %
      %   SYNTAX   [nMainIts, MiddleOne, StartTwo, MiddleTwo, EndTwo] = obj.GetIts();
      %
      %     nMainIts  - number of iterations (integer)
      %     MiddleOne - number of iterations (integer)
      %     StartTwo  - number of iterations (integer)
      %     MiddleTwo - number of iterations (integer)
      %     EndTwo    - number of iterations (integer)
      nMainIts    = this.nMainIts_;
      MiddleOne   = this.MiddleOne_;
      StartTwo    = this.StartTwo_;
      MiddleTwo   = this.MiddleTwo_;
      EndTwo      = this.EndTwo_;
    end % function

    % TODO: private (following functions too)
    function SetSmoother(this, SmootherOne, SmootherTwo)
      %SETSMOOTHER Set smoothers for each part of the 2x2 system.
      %
      %   SYNTAX   obj.SetSmoother(SmootherOne, SmootherTwo);
      %
      %     SmootherOne - smoother of the first system block (A_rr)
      %     SmootherTwo - smoother of the second system block (A_xx)
      this.SmootherOne_ = SmootherOne;
      this.SmootherTwo_ = SmootherTwo;
      this.UpdateIsSetup();
    end % function

    function SetSmootherOne(this, SmootherOne)
      %SETSMOOTHERONE Get the smoothers  of the first system block (A_rr)
      %
      %   SYNTAX   obj.SetSmootherOne(SmootherOne);
      %
      %     SmootherOne - smoother of the first system block (A_rr)
      this.SmootherOne_ = SmootherOne;
      this.UpdateIsSetup();
    end % function

    function [SmootherOne] = GetSmootherOne(this)
      %GETSMOOTHERONE Get the smoother of the first system block (A_rr)
      %
      %   SYNTAX   SmootherOne = obj.GetSmootherOne();
      %
      %     SmootherOne - smoother of the first system block (A_rr)
      SmootherOne = this.SmootherOne_;
    end % function

    function SetSmootherTwo(this, SmootherTwo)
      %SETSMOOTHERTWO Set the smoother of the second system block (A_xx)
      %
      %   SYNTAX   obj.SetSmootherTwo(SmootherTwo);
      %
      %     SmootherTwo - smoother of the second system block (A_xx)
      this.SmootherTwo_ = SmootherTwo;
      this.UpdateIsSetup();
    end % function

    function [SmootherTwo] = GetSmootherTwo(this)
      %GETSMOOTHERTWO Get the smoother of the second system block (A_xx)
      %
      %   SYNTAX   SmootherTwo = obj.GetSmootherTwo();
      %
      %     SmootherTwo - smoother of the second system block (A_xx)
      SmootherTwo = this.SmootherTwo_;
    end % function

    function SetFakeLevels(this, fakeLevelOne, fakeLevelTwo)
      %SETSMOOTHER Set fake levels for each part of the 2x2 system.
      %
      %   SYNTAX   obj.SetFakeLevels(fakeLevelOne, fakeLevelTwo);
      %
      %     fakeLevelOne - level of the first system block (A_rr)
      %     fakeLevelTwo - level of the second system block (A_xx)
      this.fakeLevelOne_ = fakeLevelOne;
      this.fakeLevelTwo_ = fakeLevelTwo;
      this.UpdateIsSetup();
    end % function

    function [fakeLevelOne, fakeLevelTwo] = GetFakeLevels(this)
      %SETSMOOTHER Get fake levels for each part of the 2x2 system.
      %
      %   SYNTAX   obj.SetFakeLevels(fakeLevelOne, fakeLevelTwo);
      %
      %     fakeLevelOne - level of the first system block (A_rr)
      %     fakeLevelTwo - level of the second system block (A_xx)
      fakeLevelOne = this.fakeLevelOne_;
      fakeLevelTwo = this.fakeLevelTwo_;
    end % function

  end % method
  methods (Access = private)

    function UpdateIsSetup(this)
      %UPDATEISSETUP Update the state of a smoother prototype.
      % This smoother is ready to be used iff:
      % - matrices A_rr and A_xx have been created
      % - SmootherOne and SmootherTwo have been setup
      % See also: SmootherPrototype.isSetup
      %
      %   SYNTAX   obj.UpdateIsSetup();
      %
      fakeLevelIsSetup = ~isempty(this.fakeLevelOne_) && ~isempty(this.fakeLevelTwo_);
      this.SetIsSetup(fakeLevelIsSetup && this.SmootherOne_.isSetup() && this.SmootherTwo_.isSetup());
    end % function

   end % method private
   methods

    function CopyParameters(this, src)
      %COPYPARAMETERS Copy the parameters of another smoother prototype.
      % See also: SmootherPrototype.CopyParameters
      %
      %   SYNTAX   obj.CopyParameters(src);
      %
      %     src - Object (SmootherPrototype of same type)

      [nMainIts, MiddleOne, StartTwo, MiddleTwo, EndTwo] = src.GetIts();
      this.SetIts(nMainIts, MiddleOne, StartTwo, MiddleTwo, EndTwo);

      % Note:
      % If this method is called by SmootherFactory::Build, then
      % - 'this', is a copy of the PreSmoother
      % - 'src' is the PostSmootherPrototype
      %
      % If Pre and Post Smoothers have the same underlying
      % SmootherOne and SmootherTwo, then we try to reuse the setup
      % done during the pre-smoother setup. It is the same idea as in
      % SmootherFactory: the parameters of the post-smoother
      % prototype (== src) are copied and this.SmootherOne/Two.isSetup()
      % switch to 'false' only if the difference between Pre and
      % Post parameters affects setup data.
      %
      % If the PostSmootherPrototype doesn't have the same internal
      % SmootherOne/Two, we throw them away and creates new objects
      % by duplicating the SmootherOne/Two of the PostSmootherPrototype.
      % In this case, the only thing we can reuse is the
      % construction of AOne and ATwo.

      srcSmootherOne = src.GetSmootherOne();
      if ~isempty(this.SmootherOne_) && ~isempty(srcSmootherOne) && ...
            (strcmp(this.SmootherOne_.GetType(), srcSmootherOne.GetType()))
        this.SmootherOne_.CopyParameters(srcSmootherOne);
      else
        this.SmootherOne_ = srcSmootherOne.Copy(); % Copy because we don't want to modify 'src'
      end

      % Same steps for SmootherTwo:
      srcSmootherTwo = src.GetSmootherTwo();
      if ~isempty(this.SmootherTwo_) && ~isempty(srcSmootherTwo) && ...
            (strcmp(this.SmootherTwo_.GetType(), srcSmootherTwo.GetType()))
        this.SmootherTwo_.CopyParameters(srcSmootherTwo);
      else
        this.SmootherTwo_ = srcSmootherTwo.Copy();
      end

      % Reuse AOne and ATwo
      [srcFakeLevelOne, srcFakeLevelTwo] = src.GetFakeLevels();
      this.SetFakeLevels(srcFakeLevelOne, srcFakeLevelTwo);

      % Now, isSetup is true if SmootherOne and SmootherTwo don't need a Setup phase
      this.UpdateIsSetup();
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

      this.Amat_ = Level.Get('A');

      Adata  = Level.Get('A').GetMatrixData();
      n      = size(Adata,1);
      NumR   = Level.Get('Arr').GetRowMap().NDOFs;
      Bsize  = Level.Get('Arr').GetRowMap().ConstBlkSize;

      % AOne (and fakeLevelOne_) is created only if needed.
      % I consider here that the 2x2 system never change between Pre and
      % Post Smoother (same AOne ATwo). Indeed, the 2x2 system is
      % fully defined by Level.Get('Arr') and we use the same for Pre
      % and Post. So we always reuse AOne and ATwo.
      if isempty(this.fakeLevelOne_)
        % TODO: do not re-create AOne. AOne == Level.Get('Arr') !
        AOne = Operator(Adata(1:NumR,1:NumR),Bsize,Bsize);
        this.fakeLevelOne_ = Level.BuildMe();
        this.fakeLevelOne_.SetLevelId(-7);
        this.fakeLevelOne_.Set('A', AOne);
      end

      % If SmootherOne is the same for pre/post smoother, avoid the
      % setup phase for the post smoother.
      % See also: SmootherFactory
      if (this.SmootherOne_.isSetup() == false) % Note: it's a useless test, must also be done at the beginning of SmootherOne.Setup() method.
        this.SmootherOne_.Setup(this.fakeLevelOne_, Specs);
      end

      % Same setup steps for SmootherTwo:
      if isempty(this.fakeLevelTwo_)
        ATwo = Operator(Adata(NumR+1:n,NumR+1:n),1,1);
        this.fakeLevelTwo_ = Level.BuildMe();
        this.fakeLevelTwo_.SetLevelId(-7);
        this.fakeLevelTwo_.Set('A', ATwo);
      end

      if (this.SmootherTwo_.isSetup() == false) % Note: useless test
        this.SmootherTwo_.Setup(this.fakeLevelTwo_, Specs);
      end

      this.UpdateIsSetup();

      % Final check
      if ~this.isSetup(), error('Internal error: Hybrid2x2Smoother cannot setup underlying smoothers.'); end
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
      % to avoid call of Setup twice (for Pre and Post smoothing)
      if ~this.isSetup(), error('apply'); end

      % Additional checks. Should NEVER append because if SmootherOne/Two not setup, then 'this' is also not setup.
      if ~this.SmootherOne_.isSetup(), error('Internal error'); end
      if ~this.SmootherTwo_.isSetup(), error('Internal error'); end

       if   isa(iu,'MultiVector'), u = iu.GetVectors();
        else                        u = iu;                   end
        if   isa(irhs,'MultiVector'), rhs = irhs.GetVectors();
        else                        rhs = irhs;               end

        mue_include

        A   = this.Amat_.GetMatrixData();
        NumR = this.GetSmootherOne.Get('A').GetRowMap().NDOFs;

        reg = (1:NumR); special = (NumR+1:length(u));

        b = rhs(special);
        if SolStatus ~= ALLZEROS, b = b - A(special,reg)*u(reg); end       % TODO: reuse the construction of sub-matrices in Apply()
        this.SmootherTwo_.SetNIts(this.StartTwo_);
        [u(special),SolStatus]= this.SmootherTwo_.Apply( u(special),b,SolStatus);

        for i=1:this.nMainIts_
         b = rhs(reg);
         if SolStatus ~= ALLZEROS, b = b - A(reg,special)*u(special);end

         this.SmootherOne_.SetNIts(this.MiddleOne_);
         [u(reg),SolStatus]= this.SmootherOne_.Apply( u(reg),b,SolStatus);

         b = rhs(special) - A(special,reg)*u(reg);
         this.SmootherTwo_.SetNIts(this.MiddleTwo_);
         [u(special),SolStatus]= this.SmootherTwo_.Apply( u(special),b,SolStatus);
        end
        this.SmootherTwo_.SetNIts(this.EndTwo_);
        [u(special),SolStatus]= this.SmootherTwo_.Apply( u(special),b,SolStatus);

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
      fprintf('%s%s: nMainIts=%d, MiddleOne=%d, MiddleTwo=%d, StartTwo=%d, MiddleTwo=%d, EndTwo=%d\n',...
              prefix,this.GetType(),this.nMainIts_, this.MiddleOne_, this.StartTwo_, this.MiddleTwo_, this.EndTwo_);
      fprintf('%s\n%s    SmootherOne=\n', prefix, prefix);
      this.SmootherOne_.Print(prefix);

      fprintf('%s\n%s    SmootherTwo=\n', prefix, prefix);
      this.SmootherTwo_.Print(prefix);
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

end % class

% Note: fakeLevelOne and fakeLevelTwo only depends on NumR.
%       NumR is stored in the level, so we always
%       assume that PreSmoo.fakeLevelOne == PostSmoo.fakeLevelOne
%       (and the same for fakeLevelTwo).

% Note: We do not try to reuse the setup phase of SmootherOne for SmootherTwo because
%       matrices One and Two are different.
