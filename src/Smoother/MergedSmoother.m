classdef MergedSmoother < SmootherBase & SmootherPrototype
  % MergedSmoother class
  % This smoother applies one smoother followed by a second smoother
  %

  properties (Access = private)
    % Parameters of the smoother
    ReverseOrder_ = false    % Indicates whether to reverse the order of
                             % the smoother invocations
    smoothers_ = []          % list of smoothers

  end

  methods
    function [this] = MergedSmoother( smootherlist )
      % MERGEDSMOOTHER Constructor
      %
      %   SYNTAX   MergedSmoother(SmootherOne,SmootherTwo);
      %
      %     Smootherone - First smoother to be applied
      %     SmootherTwo - Second smoother to be applied

      % Copy constructor // needs work
      if nargin == 1 && isa(smootherlist, class(this))
        this.Copy_(smootherlist,[]);
        return;
      end

      this.smoothers_ = smootherlist;
      this.SetType('MergedSmoother');
    end % function

    function ReverseOrder(this)
      this.ReverseOrder_ = true;
    end

    function StandardOrder(this)
      this.ReverseOrder_ = false;
    end

    function [reverseOrder] = GetReverseOrder(this)
       reverseOrder = this.ReverseOrder_;
    end

    function [smoother] = GetSmoother(this, i)
       smoother = this.smoothers_{i};
    end

    function [n] = NumSmoothers(this)
      n = size(this.smoothers_,2);
    end

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

      this.ReverseOrder_ = src.GetReverseOrder();

      %TODO: factorize with Hybrid and Factory
      for i = 1:size(src.smoothers_,2)
        srcSmoother = src.GetSmoother(i);
%         if ~isempty(this.smoothers_{i}) && ~isempty(srcSmoother) && ...
%               (strcmp(this.smoothers_{i}.GetType(), srcSmoother.GetType()))
%           this.smoothers_{i}.CopyParameters(srcSmoother);
%         else
          this.smoothers_{i} = srcSmoother.Copy(); % Copy because we don't want to modify 'src'
%         end
      end

      setup = true;
      for i = 1:this.NumSmoothers()
        if ( ~this.smoothers_{i}.isSetup() )
          setup = false;
          break
        end
      end
      this.SetIsSetup(setup)

    end % function

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setup phase
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function Setup(this, Level)
      %SETUP Run the setup phase of the smoother.
      % See also: SmootherPrototype.Setup
      %
      %   SYNTAX   Amat = obj.Setup(Level);
      %
      %     Level - level of the MG hierachy (Level)
      if this.isSetup(), return; end

      for i = 1:this.NumSmoothers()
        smo = this.smoothers_{i};
        if ~smo.isSetup()
          smo.Setup(Level);
        end
      end

      setup = true;
      for i = 1:this.NumSmoothers()
        if ( ~this.smoothers_{i}.isSetup() )
          setup = false;
          break
        end
      end
      this.SetIsSetup(setup)

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

      if (~this.ReverseOrder_)
        for i = 1:this.NumSmoothers()
          [u, SolStatus] = this.smoothers_{i}.Apply(u, rhs, SolStatus);
        end
      else
        for i = this.NumSmoothers():-1:1
          [u, SolStatus] = this.smoothers_{i}.Apply(u, rhs, SolStatus);
        end
      end;

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
      fprintf('%sMergedSmoother: ReverseOrder=%g\n',...
               prefix,this.ReverseOrder_);
      for i = 1:this.NumSmoothers()
        fprintf('subsmoother %d\n', i);
        this.smoothers_{i}.Print(prefix);
      end

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

%TODO: Set() methods