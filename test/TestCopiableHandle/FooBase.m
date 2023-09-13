classdef FooBase < CopiableHandle
  properties
    BaseA_= 0;
  end
  properties (Access = private)
    BaseB_= 0;
  end
  
  methods
    function this = FooBase(a,b)
      % Copy constructor
      if nargin == 1 && isa(a, class(this)), this.Copy_(a,[]); return; end
      
      % Default constructor
      if nargin > 0
       this.SetBaseA(a);
      end
      if nargin > 1
       this.SetBaseA(b);
      end
      
    end
    
    function SetBaseA(this, BaseA)
      this.BaseA_ = BaseA;
    end
    
    function [BaseA] = GetBaseA(this)
      BaseA = this.BaseA_;
    end
    
    function SetBaseB(this, BaseB)
      this.BaseB_ = BaseB;
    end
    
    function [BaseB] = GetBaseB(this)
      BaseB = this.BaseB_;
    end

  end % methods
  
  methods (Access = protected)  

    function Copy_(this, src, mc)
      [cmd, data, mc] = this.CopyCmd_(src,mc);
      eval(cmd);    
    end
    
  end % methods
  
end
