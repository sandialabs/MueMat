classdef OperatorView < CopiableHandle
%SINGLEVIEWOPERATOR Class that holds data for a specific view of an Operator
%
% This class maintains a specific row map, column map, apply function, diagonal, and eigenvalue estimate for
% a matrix.  The only data it maintains is for the diagonal.  The actual matrix data is stored with the Operator
% class.

  properties
    rowmap_
    colmap_
    apply_
    applyinverse_
  end %properties

  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Constructor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [this] = OperatorView(rowmap,colmap,matvec)

      % Copy constructor
      if nargin == 1 && isa(rowmap, class(this)), this.Copy_(rowmap,[]); return; end
      %

      % OperatorView constructor.
      if ~varexist('rowmap')
        error('OperatorView ctor missing "rowmap"\n');
      elseif ~varexist('colmap')
        error('OperatorView ctor missing "colmap"\n');
      elseif ~varexist('matvec')
        error('OperatorView ctor missing "matvec"\n');
      end
      this.SetApply(matvec);
      this.SetRowMap(rowmap);
      this.SetColMap(colmap);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Set methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function SetRowMap(this,rowmap)
      % Assign row map.
      this.rowmap_ = rowmap;
    end
    %%%%
    function SetColMap(this,colmap)
      % Assign column map.
      this.colmap_ = colmap;
    end
    %%%%
    function SetApply(this,apply)
      % Assign matrix-vector function handle.
      this.apply_ = apply;
    end
    %%%%
    function SetApplyInverse(this,fcn,overwrite)
      % Assign apply-inverse function handle.
      if ~varexist('overwrite') || overwrite == true
        this.applyinverse_ = fcn;
      end
    end
    %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Get methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [rowmap] = GetRowMap(this,rowmap)
      % Return row map.
      rowmap = this.rowmap_;
    end
    %%%%
    function [colmap] = GetColMap(this,colmap)
      % Return column map.
      colmap = this.colmap_;
    end
    %%%%
    function [apply] = GetApply(this,apply)
      % Return matrix-vector function handle.
      apply = this.apply_;
    end
    %%%%
    function [fcn] = GetApplyInverse(this,fcn)
      % Return apply-inverse function handle.
      fcn = this.applyinverse_;
    end

    %%%

  end %methods

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

end %classdef OperatorView
