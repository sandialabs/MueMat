classdef Map < CopiableHandle
% Like epetra BlockMaps in that they represent both point and block size information
  properties(Access = private)
     NDOF_         %    Total Number of DOFs in vector
     NNodes_       %    Number of nodes in vector
     MaxBlkSize_   %    Maximum number of degrees-of-freedom per node.
     type_         %    = 'ConstBlk': all nodes have the same number of DOFs
                   %    = 'VarBlk   : nodes have varying numbers of DOFs
     ConstBlkSize_ %    DOFs per node when type == 'ConstBlk'
     VBlkPtr_      %    VarBlkPtr(i):VarBlkPtr(i+1)-1 gives DOF ids in
                   %    ith node when  type == 'VarBlk'
  end

  methods
    function this = Map(NNodes,Bsize)
       % construct a BlockMap where the block sizes can either be constant or variable.
       % In particular, if the Bsize argument is a scalar, a constant block size is assumed.
       % Otherwise, a variable block sized map is created.

       % Copy constructor
       if nargin == 1 && isa(NNodes, class(this)), this.Copy_(NNodes,[]); return; end
       %

       this.NNodes_          = NNodes;

       if length(Bsize) == 1,
          this.ConstBlkSize_ = Bsize;
          this.MaxBlkSize_   = Bsize;
          this.type_         = 'ConstBlk';
          this.VBlkPtr_      = -1;
          this.NDOF_         = NNodes*Bsize;

       else

          if length(Bsize) ~= NNodes,
             fprintf('Map: the number of nodes does not match the length of Bsize\n');
             keyboard;
          end

          this.ConstBlkSize_ = -1;
          this.MaxBlkSize_   = max(Bsize);
          this.type_         = 'VarBlk';

          VPtr               = zeros(NNodes+1,1);
          VPtr(1)            = 1;
          for i=1:NNodes, VPtr(i+1) = VPtr(i) + Bsize(i); end

          this.VBlkPtr_     = VPtr;
          this.NDOF_        = sum(Bsize);

       end
    end% function
    function value = NNodes(this)
    % Number of nodes (or blocks) associated with this map.
       value = this.NNodes_;
    end
    function value = NDOFs(this)
    % Number of degrees-of-freedom (or points) associated with this map.
       value = this.NDOF_;
    end
    function value = MaxBlkSize(this)
    % Maximum block size associated with this map.
       value = this.MaxBlkSize_;
    end
    function value = Vptr(this)
    % Returns a vector where the ith value indicates the degree-of-freedom associated with the start of the ith block for variable block size maps (otherwise returns -1 for constant block size maps).
       value = this.VBlkPtr_;
    end
    function size = ConstBlkSize(this)
    % Returns the block size of a map with constant block size (otherwise returns -1 for variable block size maps).
       size = this.ConstBlkSize_;
    end
    function bool = HasConstBlkSize(this)
    % Returns true if constant block size map, otherwise returns false
       if this.ConstBlkSize() == -1, bool = false; else bool = true; end
    end
    function bool = HasVariableBlkSize(this)
    % Returns true if variable block size map, otherwise returns false
       if this.VBlkPtr_ == -1, bool = false; else bool = true; end
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
