% Generic Smoother Factory for generating the smoothers of the MG hierarchy
%
% This factory is generic and can produce any kind of Smoother.
%
% The design of the Smoother factory is based on the prototype design
% pattern.  The smoother factory uses prototypical instances of
% smoothers and prototypes are cloned to produce the smoothers of the MG
% hierarchy.
%
% See also:
% - http://en.wikipedia.org/wiki/Factory_method_pattern
% - http://en.wikipedia.org/wiki/Prototype_pattern
% - http://www.oodesign.com/prototype-pattern.html
%
% There is one prototype for the pre-smoother and one for the
% post-smoother. Thus, this factory can produce two different smoothers
% for pre and post-smoothing. Prototypes are stored in this factory.
%
% The type of smoother to create is determined by the prototypical
% instances. Prototypes store their own parameters (nits, omega) and the
% factory doesn't have any knowledge about that.
%
% See also: SmootherFactory.SmootherFactory, SmootherFactory.Build, Smoother
%
classdef SmootherFactory < SmootherFactoryBase
   properties (Access = private)
     PreSmootherPrototype_  % Prototype for pre-smoothers
     PostSmootherPrototype_ % Prototype for post-smoothers
   end
   methods
     function [this] = SmootherFactory(PreSmootherPrototype, ...
                                   PostSmootherPrototype)
      %SMOOTHERFACTORY Constructor
      % The constructor takes as arguments pre-configured smoother
      % object(s). They are used as prototypes by the factory.
      %
      %   SYNTAX   SmootherFactory(PreSmootherPrototype, PostSmootherPrototype);
      %
      %     PreSmootherPrototype  - prototype for pre-smoothers  (SmootherPrototype)
      %     PostSmootherPrototype - prototype for post-smoothers (SmootherPrototype,optional,default=PreSmootherPrototype)
      %
      % EXAMPLES:
      %
      %   SmooFactory = SmootherFactory(ChebySmoother(5,1/30))
      %
      % or
      %
      %   nIts        = 5;
      %   lambdaRatio = 1/30;
      %   Smoo        = ChebySmoother()
      %   Smoo        = Smoo.SetIts(nIts);
      %   Smoo        = Smoo.SetLambdaRatio(1/30);
      %   SmooFactory = SmootherFactory(Smoo);
      %
      % To use different smoothers for pre and post smoothing, two prototypes can be passed in as argument:
      %   PreSmoo     = ChebySmoother(2, 1/30)
      %   PostSmoo    = ChebySmoother(10,1/30)
      %   SmooFactory = SmootherFactory(PreSmoo, PostSmoo);
      %
      % [] is also a valid smoother which do nothing:
      %   PostSmoo    = ChebySmoother(10,1/30)
      %   SmooFactory = SmootherFactory([], PostSmoo);
      %
      % See also: SmootherFactory.Build

      % Copy constructor
      if nargin == 1 && isa(PreSmootherPrototype, class(this)), this.Copy_(PreSmootherPrototype,[]); return; end
      %

      if nargin==1
       % If post-smoother is not explicitily specified, it is the same as pre-smoother
       this.PreSmootherPrototype_  = PreSmootherPrototype;
       this.PostSmootherPrototype_ = PreSmootherPrototype; % handle copy
      elseif nargin==2
       this.PreSmootherPrototype_  = PreSmootherPrototype;
       this.PostSmootherPrototype_ = PostSmootherPrototype;
      else
       fprintf('Incorrect constructor arguments for SmootherFactory()\n');
       keyboard;
      end
     end

     function SetSmootherPrototypes(this, PreSmootherPrototype, PostSmootherPrototype)
      %SETSMOOTHERPROTOTYPES Set smoother prototypes
      %
      %   SYNTAX   obj.SetSmootherPrototypes(PreSmootherPrototype, PostSmootherPrototype);
      %
      %     PreSmootherPrototype  - prototype for pre-smoothers  (SmootherPrototype)
      %     PostSmootherPrototype - prototype for post-smoothers (SmootherPrototype)

      %TODO: same behaviour as constructor (post=optional). Factorize with constructor
      this.PreSmootherPrototype_  = PreSmootherPrototype;
      this.PostSmootherPrototype_ = PostSmootherPrototype;
     end

     function [PreSmootherPrototype,PostSmootherPrototype] = GetSmootherPrototypes(this)
      %GETSMOOTHERPROTOTYPES Get smoother prototypes
      %
      %   SYNTAX   [PreSmootherPrototype, PostSmootherPrototype] = obj.GetSmootherPrototypes();
      %
      %     PreSmootherPrototype  - prototype for pre-smoothers  (SmootherPrototype)
      %     PostSmootherPrototype - prototype for post-smoothers (SmootherPrototype)

      PreSmootherPrototype  = this.PreSmootherPrototype_;
      PostSmootherPrototype = this.PostSmootherPrototype_;
     end

     function SetNeeds(this, Level)
       % Obtain any cross factory specifications
       if ~isempty(this.PreSmootherPrototype_ ), this.PreSmootherPrototype_.SetNeeds(Level); end;
       if ~isempty(this.PostSmootherPrototype_), this.PostSmootherPrototype_.SetNeeds(Level); end;
     end

     function [PreSmoo,PostSmoo] = Build(this, Level) %, Specs)
      %BUILD Creates pre and post smoothers.
      %
      %   SYNTAX   [PreSmoo, PostSmoo] = obj.Build(Level, Specs);
      %
      %     Level    - level of the MG hierachy (Level)
      %     Specs    - specifications (CrossFactory)
      %     PreSmoo  - pre-smoother  (SmootherBase)
      %     PostSmoo - post-smoother (SmootherBase)
      %
      % Factory.Build() clones its prototypes and calls Setup() on the
      % new objects to create fully functional smoothers for the current
      % hierarchy level. Note that smoothers do their own setup
      % (ie: ILUSmoother.Setup() computes the ILU factorization).
      %
      % If pre and post smoother are identical, the Setup() phase
      % is done only once: to create the post-smoother, the
      % pre-smoother is duplicated and its parameters are changed
      % according to the parameters of the post-smoother prototype.
      %
      % If the parameters of pre and post smoothers are not the same,
      % the Setup() phase is also done only once when parameters
      % don't change the result of the setup computation.
      %
      if nargout() < 2, error('nargout'); end

      PreSmoo  = [];
      PostSmoo = [];

      % PRE SMOOTHER
      if ~isempty(this.PreSmootherPrototype_)

       %Copy the prototype which stores the parameters of the pre-smoother.
       PreSmoo = this.PreSmootherPrototype_.Copy();

       if this.GetOutputLevel() > 0
         PreSmoo.Print(['(level ' num2str(Level.GetLevelId()) ') ']);
       end

       % Run the setup phase on the copy.
       PreSmoo.Setup(Level); %, Specs);

      end

      % POST SMOOTHER
      if ~isempty(this.PostSmootherPrototype_)

       % Is post-smoother of the same type as pre-smoother ?
       if ~isempty(this.PreSmootherPrototype_) && ...
             (strcmp(this.PreSmootherPrototype_.GetType(), this.PostSmootherPrototype_.GetType()))

         % YES: post-smoother == pre-smoother
         % => copy the pre-smoother to avoid the setup phase of the post-smoother.
         PostSmoo = PreSmoo.Copy();

         % If the post-smoother parameters are different from
         % pre-smoother, the parameters stored in the post-smoother
         % prototype are copied in the new post-smoother object.
         PostSmoo.CopyParameters(this.PostSmootherPrototype_);

         % If parameters don't influence the Setup phase (it is the case
         % for Jacobi, Chebyshev...), PostSmoo is already setup. Nothing
         % more to do. On the case of ILU, parameters of the smoother
         % are in fact the parameters of the Setup phase. The call to
         % CopyParameters reset the smoother (only if parameters are
         % different) and we must call Setup() again.
         PostSmoo.Setup(Level); %, Specs); % In general, do nothing.

         % TODO: if CopyParameters do not exist, do setup twice.

       else

         % NO: post-smoother != pre-smoother
         % Copy the prototype and run the setup phase.

         PostSmoo = this.PostSmootherPrototype_.Copy();
         PostSmoo.Setup(Level); %s, Specs);

       end
      end
     end

   end %public methods

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

end % classdef

