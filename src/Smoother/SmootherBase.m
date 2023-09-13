classdef SmootherBase < VerboseObject
% Base class for smoothers
% This has the signature for the required Apply function and contains data that is generic across all smoothers.

  properties(Access = private)
   type_ = 'undefined' % identify the type of a smoother object (string)
  end

  methods

    %TODO: remove SetNIts. It's not 'generic'.
    function SetNIts(this, nIts)
      %SETNITS Prototype of the method to setup the number of iteration of the smoother.
      % This method is needed by Hybrid Smoother. TODO: remove (as DirectSolve don't use it).
      % This superclass method do nothing and should be reimplemented it in derived classes.
      %
      %   SYNTAX   obj.SetNIts(nIts);
      %
      %     nIts - number of iterations (integer)

      % DO NOTHING
    end % function

    function [type] = GetType(this)
      %GETTYPE Get the smoother type.
      %
      %   SYNTAX   type = obj.GetType();
      %
      %     type - identify the type of a smoother object (string)
      type  = this.type_;
    end % function

  end % methods

  methods(Access = protected)
    function SetType(this, type)
      %SETTYPE Set the smoother type.
      % This method must be called by constructors of derived classes.
      %
      %   SYNTAX   obj.SetType(type);
      %
      %     type - identify the type of a smoother object (string)
      this.type_ = type;
    end % function

  end %methods

  methods (Abstract)
    [iu,SolStatus] = Apply(this, iu, irhs, SolStatus);
    %APPLY Apply the smoother
    %
    %    SYNTAX    obj.Apply(iu, irhs, InitGuessStatus)
    %
    %      iu        - vector to smooth (could be a MultiVector or in Matlab format)
    %      irhs      - right-hand side  (could be a MultiVector or in Matlab format)
    %      SolStatus - when InitGuessStatus==ALLZEROS, iu is assumed to be all zeros and initial matvec operation could be avoid (optional,default=NOTALLZEROS)

    Print(this,prefix);
    %PRINT Print information about the smoother
    %
    %   SYNTAX   obj.Print()
    %
    %     prefix  - optional string that is prepended to each print line

  end % abstract methods

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
