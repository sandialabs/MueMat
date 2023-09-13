% Example of class using CopiableHandle
%
% See also: CopiableHandle

classdef Foo < CopiableHandle & FooBase
  properties
    FooA_= 0;
  end
  properties (Access = private)
    FooB_= 0;
  end
  
  methods
    function this = Foo(a,b)
      % Copy constructor
      if nargin == 1 && isa(a, class(this)), this.Copy_(a,[]); return; end
      
      % Default constructor
      if nargin > 0
       this.SetFooA(a);
      end
      if nargin > 1
       this.SetFooA(b);
      end
      
    end
    
    function SetFooA(this, FooA)
      this.FooA_ = FooA;
    end
    
    function [FooA] = GetFooA(this)
      FooA = this.FooA_;
    end
    
    function SetFooB(this, FooB)
      this.FooB_ = FooB;
    end
    
    function [FooB] = GetFooB(this)
      FooB = this.FooB_;
    end

  end % methods
  
  methods (Access = protected)  

    function Copy_(this, src, mc)
      [cmd, data, mc] = this.CopyCmd_(src,mc);
      eval(cmd);    
    end
    
  end % methods
  
end
