classdef Operator < CopiableHandle
%OPERATOR Provides single interface for matrix access
%
% This class provides a unified interface for creating, accessing, and
% manipulating matrices.  The main idea is that the user can access the matrix
% as if it were stored in different formats, independent of the actual underlying
% storage format.  Each way of accessing the matrix is called a "view".  A view
% allows for associating with the matrix a particular row map, column map,
% matrix-vector multiply, and notion of matrix diagonal.
%
% This class is fairly lightweight.  It keeps the underlying Matlab matrix, which
% view is the current one, and a ViewTable.
% On construction, an Operator has just one view, referred to as the default view.
%
% Getting and setting attributes of an Operator are through Get/Set methods, e.g.,
%
%     rowmap = A.GetRowMap();
%     dat = sprand(20,20);
%     B.SetMatrixData(dat);
%
% Note:
% - 'current' refers to the current view.
% - 'default refers to the defaultView_ (ie: 'point' or 'block')
% - 'current' and 'default' are reserved keywords. You cannot create views named 'default' or 'current'.
%
% Often, by a slight abuse of the notation, view labels (= string
% which identify the view in the hashtable) are just named "a view". As
% MATLAB is an untyped languages, it may be confusing.
%

% TODO: simplify defaultView_ currentView_ pointView_ ... etc.
% TODO: add label everywhere: pointViewLabel_

  properties (Access = protected)
    label_              % user-specified label, used only in printing
    matrixData_         % matrix data (native Matlab matrix)
    currentView_        % current view (string label, key into viewTable_) - Note: [] == currentView
    defaultView_        % the view associated with inital Operator construction
    pointView_          % the view associated with point diagonal
    operatorViewTable_  % table of available views
    diagonalViewTable_  % table of available views
    applyFcnHandle_     % handle to current view's matvec (to avoid indirection)
    rowmapPtr_          % current view's rowmap (to avoid indirection)
    matrixPattern_      % non-zero pattern (native Matlab matrix with 1 entries)
  end %properties

  methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Methods implemented:
    %
    % Constructors
    %  function [this] = Operator(A,rowmap,colmap,matvec,label)
    %  function [this] = Operator(A,RowBlkSize,ColBlkSize)
    %
    % Numerical operations
    %  function  [OutVec] = MatlabApply(Op,Vec, Subset)
    %  function  [Cmat] = plus(Amat,Bmat)
    %  function  [Cmat] = minus(Amat,Bmat)
    %  function  [this] = Scale(this,scalar)
    %  function  [Cmat] = mtimes(Amat,Bmat)
    %  function  [Atrans] = ctranspose(Amat)
    %  function  [Cmat] = MatMatSolve(Amat,Bmat)
    %
    % Comparison operations
    %  function  [TorF] = eq(Amat,Bmat)
    %
    % Utilities
    %  function  [] = Print(this)
    %  function  [this] = SetLabel(this,label)
    %  function  [label] = GetLabel(this)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [this] = Operator(A,varargin)
      %OPERATOR Constructor
      %
      %   SYNTAX   obj = Operator(A, RowBlkSize, ColBlkSize, label);
      %
      %     A          - a MATLAB matrix
      %     RowBlkSize - size of row block (integer, optional, default=1)
      %     ColBlkSize - size of row block (integer, optional, default=1)
      %     label      - name of Operator (string, optional, default='unamed')
      %
      %   OR
      %
      %   SYNTAX   obj = Operator(A, rowmap, colmap, matvec, label);
      %
      %     A      - a MATLAB matrix
      %     rowmap - a row map (Map)
      %     colmap - a column map (Map)
      %     matvec - matrix-vector multiplication function (handle)
      %     label  - name of Operator (string, optional, default='unamed')

      %
      % See also MatlabAndMapToOperator and MatlabToOperator for
      % optional arguments
       % Copy constructor
       if nargin == 1 && isa(A, class(this)), this.Copy_(A,[]); return; end
       %

      if nargin >= 2 && isa(varargin{1},'Map')
       this.MatlabAndMapToOperator(A,varargin{:});
      else
       this.MatlabToOperator(A,varargin{:});
      end

    end

    end % method

    methods(Access = private) % Overloaded constructors
      function MatlabAndMapToOperator(this,A,rowmap,colmap,matvec,label)
      %MATLABANDMAPTOOPERATOR Constructor sub-function
      %   SYNTAX   obj.MatlabAndMapToOperator(A, rowmap, colmap, matvec, label);
      %
      %     A      - a MATLAB matrix
      %     rowmap - a row map (Map)
      %     colmap - a column map (Map)
      %     matvec - matrix-vector multiplication function (handle)
      %     label  - name of Operator (string, optional, default='unamed')

      % basic idea
      %   determine whether rowmap is point-based
      %      yes -- call the view "point" or something and save it in the private variable defaultView_
      %      no  -- 1) call the view "block" and save it in the private variable defaultView_
      %             2) create a point view and call it point
      %   two alternatives
      %       1) make a new method, SwitchtoDefaultView, that will use defaultView_
      %       2) make a new method, GetDefaultView, and then do this.SwitchToView(this.GetDefaultView());
      if ~varexist('A')
       fprintf('Operator ctor missing "A"\n');
        keyboard;
      end
      if ~varexist('rowmap')
        fprintf('Operator ctor missing "rowmap"\n');
        keyboard;
      elseif ~varexist('colmap')
        fprintf('Operator ctor missing "colmap"\n');
        keyboard;
      elseif ~varexist('matvec')
        fprintf('Operator ctor missing "matvec"\n');
        keyboard;
      end
      this.matrixData_ = A;

      this.operatorViewTable_ = containers.Map(); % table of available views
      this.diagonalViewTable_ = containers.Map(); % table of available views

      if (rowmap.ConstBlkSize() == 1), viewName = 'point';
      else                             viewName = 'block'; end
      this.defaultView_ = viewName;
      this.currentView_ = viewName; % needed because CreateView() use CurrentView() (to update applyFcnHandle_ and rowmapPtr_).
      this.CreateView(viewName,rowmap,colmap,matvec);
      %this.SwitchToView(viewName);

      if strcmp(viewName,'point') == false
        pointrowmap = Map(rowmap.NDOFs(),1);
        pointcolmap = Map(colmap.NDOFs(),1);
        this.CreateView('point',pointrowmap,pointcolmap,matvec);
        % TODO: create also point diagonal ?
      end
      this.pointView_ = 'point';

      if varexist('label')
        this.label_ = label;
      else
        this.label_ = 'unnamed';
      end
    end

    % Operator constructor.
    % Convert a regular matlab file into a MueMat structure.
    %
    % At present, this only works for constant block size.
    % It would be great if we could use overloading so that
    % we could also handle variable block sizes.
    %
    % varargin = label

    function MatlabToOperator(this, A, RowBlkSize, ColBlkSize, label)
      %MATLABTOOPERATOR Constructor sub-function
      %
      %   SYNTAX   obj.MatlabToOperator(A, RowBlkSize, ColBlkSize, label);
      %
      %     A          - a MATLAB matrix
      %     RowBlkSize - size of row block (integer)
      %     ColBlkSize - size of row block (integer)
      %     label      - name of Operator (string, optional, default='unamed')
      if ~varexist('RowBlkSize'), RowBlkSize = 1; end
      if ~varexist('ColBlkSize'), ColBlkSize = 1; end

      %
      %   Setup RowMap
      %
      RowMap = Map(size(A,1)/RowBlkSize, RowBlkSize);

      %
      %   Setup ColMap
      %
      ColMap = Map(size(A,2)/ColBlkSize, ColBlkSize);

      %
      %   Make MuMat matrices
      %
      if varexist('label')
       this.MatlabAndMapToOperator(A,RowMap,ColMap,@MatlabApply,label);
      else
       this.MatlabAndMapToOperator(A,RowMap,ColMap,@MatlabApply);      % update applyFcnHandle_ and rowmapPtr_
      % (if view ==this.CurrentView(), applyFcnHandle_ and rowmapPtr_ must be updated)
      this.applyFcnHandle_ = this.operatorViewTable_(this.CurrentView()).GetApply();
      this.rowmapPtr_ = this.operatorViewTable_(this.CurrentView()).GetRowMap();
      end
    end

    end % methods(Static, Access = private)

    methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [viewLabel] = CurrentView(this)
      %CURRENTVIEW Return the name of the current view.
      %
      %   SYNTAX   view = obj.CurrentView();
      %
      %     viewLabel - current operator view
      viewLabel = this.currentView_;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [viewList] = GetViewList(this)
      %GETVIEWLIST Return a cell array of all the view names for this Operator.
      %
      %   SYNTAX   ViewList = obj.GetViewList();
      %
      %     ViewList - list of views (cell array)
      viewList = keys(this.operatorViewTable_);
    end

    %TODO: remove ???
    function [viewList] = GetDiagonalViewList(this)
      %GETVIEWLIST Return a cell array of all the view names for the operator diagonal.
      %
      %   SYNTAX   ViewList = obj.GetViewList();
      %
      %     ViewList - list of views (cell array)
      viewList = keys(this.diagonalViewTable_);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function oldView = SwitchToView(this,newView)
      %SWITCHTOVIEW Switch to a new view.  Returns the old view name.
      %
      %   SYNTAX   oldView = obj.SwitchToView(newView);
      %
      %     newView - new operator view (string)
      %     oldView - old operator view (string)
      if strcmp(newView,'current'), oldView = 'current'; return; end
      if strcmp(newView,'default'), newView = this.defaultView_; end

      oldView = this.currentView_;
      this.currentView_ = newView;
      this.applyFcnHandle_ = this.operatorViewTable_(newView).GetApply();
      this.rowmapPtr_ = this.operatorViewTable_(newView).GetRowMap();
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function oldView = SwitchToDefaultView(this)
      %SWITCHTODEFAULTVIEW Switch to default view.  Returns the old view name.
      %
      %   SYNTAX   oldView = obj.SwitchToDefaultView();
      %
      %     oldView - old operator view (string)
      oldView = this.SwitchToView(this.defaultView_);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function oldView = SwitchToPointView(this)
      %SWITCHTOPOINTVIEW Switch to point view.  Returns the old view name.
      %
      %   SYNTAX   oldView = obj.SwitchToPointView();
      %
      %     oldView - old operator view (string)
      oldView = this.SwitchToView(this.pointView_);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function CreateView(this,newView,rowmap,colmap,matvec)
      %CREATEVIEW Create a new view of this Operator.  This does *not* switch the view.
      % Use CreateView() instead of SetView() or use SwitchToView()
      % afterwards. /// ??
      %
      %   SYNTAX   obj.CreateView(newView, rowmap, colmap, matvec);
      %
      %     newView - a operator view (string)
      %     rowmap  - a row map (Map)
      %     colmap  - a column map (Map)
      %     matvec  - matrix-vector multiplication function (handle)

      if strcmp(newView,'current') || strcmp(newView,'default')
        error(['''current'' and ''default'' are reserved keywords. ' ...
               'You cannot named a view ''current'' or ''default''']);
      end;

      %TODO: opView-> rename view
      opView   = OperatorView(rowmap,colmap,matvec);

      % TODO: if view already exists -> error();
      this.operatorViewTable_(newView) = opView;

      % update applyFcnHandle_ and rowmapPtr_
      % (if newView ==this.CurrentView(), applyFcnHandle_ and rowmapPtr_ must be updated)
      this.applyFcnHandle_ = this.operatorViewTable_(this.CurrentView()).GetApply();
      this.rowmapPtr_ = this.operatorViewTable_(this.CurrentView()).GetRowMap();
    end

    function SetView(this,viewLabel,rowmap,colmap,matvec)
      %SETVIEW Change the view to the one associated with 'viewLabel'.
      % If this view does not exist, a new one is created using
      % rowmap,colmap,matvec arguments. If the view exists, these additional
      % arguments are ignored.
      %
      %   SYNTAX   obj.SetView(view, rowmap, colmap, matvec);
      %
      %     viewLabel - an operator view
      %     rowmap - a row map (Map)
      %     colmap - a column map (Map)
      %     matvec - matrix-vector multiplication function (handle)

      % TODO: use varargin for rowmap,colmap,matvec

      % TODO: add the following check when the view doesn't exist:
%       if strcmp(newView,'current') || strcmp(view,'default')
%         error(['''current'' and ''default'' are reserved keywords. ' ...
%                'You cannot named a view ''current'' or ''default''')]);
%       end;

      opView   = OperatorView(rowmap,colmap,matvec);

      % Same code as CreateView but no error() if view exists.
      this.operatorViewTable_(viewLabel) = opView;%,@this.GetMatrixData);

      % Switch view
      this.SwitchToView(viewLabel); % update currentView_ applyFcnHandle_ and rowmapPtr_
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [OutVec] = MatlabApply(Op,Vec, Subset)
      %MATLABAPPLY Multiply the matrix with a vector/matrix using Matlab's native multiplication.
      %
      %   SYNTAX   OutVec = obj.MatlabApply(Op, Vec, Subset);
      %
      %     Op     - an Operator
      %     Vec    - a vector (MATLAB vector)
      %     Subset - Subset of indices on which to operate (optional, default=1:end)
      %     OutVec - output vector
      if nargout() ~= 1, error('MatlabApply method returns vector'); end
      indices = Subset2DOF(Subset,Op.GetRowMap());
      OutVec  = Op.matrixData_(indices,:)*Vec;
    end %MatlabApply

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function  [Cmat] = plus(Amat,Bmat)
      %PLUS Overload of the binary addition operator +.
      %
      %   SYNTAX   Cmat = Amat + Bmat;
      %
      %     Amat - an Operator
      %     Bmat - an Operator or a MATLAB matrix
      %     Cmat - an Operator
      if ~isa(Amat,'Operator')
        tmp = Bmat;
        Bmat = Amat;
        Amat = tmp;
      end
      if ~strcmp(func2str(Amat.GetApply()),func2str(@MatlabApply)),
        fprintf('plus:: Only matrices for MatlabApply matrices.\n');
        keyboard;
      end
      if isa(Bmat,'Operator')
       if xor(Amat.GetRowMap().HasConstBlkSize(),Bmat.GetRowMap().HasConstBlkSize()) || ...
          xor(Amat.GetColMap().HasConstBlkSize(),Bmat.GetColMap().HasConstBlkSize()),
         fprintf('plus:: Only matrices with identical blocking implemented.\n');
         keyboard;
       end
        Bdata = Bmat.GetMatrixData();
      else
        Bdata = Bmat;
      end

      Cdata = Amat.GetMatrixData() + Bdata;
      Cmat  = Operator(Cdata, Amat.GetRowMap(), Amat.GetColMap(), @MatlabApply);
    end %operator plus

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function  [Cmat] = uminus(Amat)
      %UMINUS
      %
      %   SYNTAX   Cmat = -Amat;
      %
      %     Amat - an Operator
      %     Cmat - an Operator
      Cmat = Copy(Amat);
      Cmat.matrixData_ = -1.0 * Cmat.matrixData_;
    end %operator uminus

    function [nrow, ncols] = size(Amat,dim)
      %TODO: to be fixed
      %SIZE Overload of the size function.
      %
      %   SYNTAX   [nrow, ncols] = obj.size(Amat, dim);
      %
      %     Amat  - an Operator
      %     dim   - get size of the dimension specified by the scalar dim (integer, optional, see built-in function)
      %     nrow  - number of rows (integer)
      %     ncols - number of columns (integer)
      if ~varexist('dim'), dim = 1; end
      if dim>2
        nrow = 1; % mimic the builtin size behavior.
        return;
      end
      if dim<1
        error('Operator size: dimension must be positive integer.');
      end
      [nrow,ncols] = size(Amat.matrixData_);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function  [Cmat] = minus(Amat,Bmat)
      %MINUS Overload of the binary subtraction operator -.
      %
      %   SYNTAX   Cmat = Amat - Bmat;
      %
      %     Amat - an Operator
      %     Bmat - an Operator
      %     Cmat - an Operator
      Cmat = Amat + (-Bmat);
    end %operator minus

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [Cmat] = mldivide(Amat,Bmat)
      %MLDIVIDE Overload of the binary left division operator \.
      %
      %   SYNTAX   Cmat = Amat \ Bmat;
      %
      %     Amat - an Operator
      %     Bmat - an Operator, a MultiVector or a MATLAB matrix
      %     Cmat - an Operator, a MultiVector or a MATLAB matrix
      ncA = size(Amat,2);
      [nrB,ncB] = size(Bmat);
      if ncA ~= nrB
        error('Operator left divide: Inner dimensions must agree.');
      end

      if isa(Amat,'Operator'), Adata = Amat.GetMatrixData();
      else                     Adata = Amat;           end
      if isa(Bmat,'Operator')
        Bdata = Bmat.GetMatrixData();
        Cdata = Adata \ Bdata;
        %return an Operator
        Cmat = Operator(Cdata,Bmat.GetRowMap(),Bmat.GetColMap(),@MatlabApply);
      elseif isa(Bmat,'MultiVector')
        Bdata = Bmat.GetVectors();
        Cdata = Adata \ Bdata;
        %return a MultiVector
        Cmat = MultiVector([],[],Cdata);
      else
        Bdata = Bmat;
        %return a Matlab matrix
        Cmat = Adata \ Bdata;
      end
    end %operator mldivide

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Scale(this,scalar)
      %SCALE Multiply all entries of Operator by given scalar.
      %
      %   SYNTAX   obj.Scale(scalar);
      %
      %     scalar - value (double)
      this.matrixData_ = scalar*this.matrixData_;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function  [Cmat] = mtimes(Amat,Bmat)
      %MTIMES Overload of the binary multiplication operator *
      % Bmat can be (but does not have to be) an Operator.
      % If Bmat is a MultiVector, so are Cmat.
      %
      %   SYNTAX   Cmat = Scalar * Amat
      %            Cmat = Amat   * Bmat
      %
      %     Amat   - an Operator
      %     Bmat   - an Operator, a MultiVector or a MATLAB matrix
      %     Cmat   - an Operator, a MultiVector or a MATLAB matrix
      %
      % If Bmat is a MultiVector, so are Cmat.
      % If Bmat is a MATLAB, so are Cmat. % TODO: should be changed

      if ~isa(Amat,'Operator')

        % scalar * Mat
        Cmat = Operator(Amat * Bmat.GetMatrixData(), Bmat.GetRowMap(), Bmat.GetColMap(), @MatlabApply);

      else

        % Mat * Mat
        matvec = Amat.GetApply();
        if isa(Bmat,'Operator')
          if ~strcmp(func2str(Bmat.GetApply()),func2str(@MatlabApply)),
            error('mtimes:: Bmat must be a MatlabApply matrix.\n');
          end
          matrixData = sparse(matvec(Amat, Bmat.matrixData_, []));
          Cmat = Operator(matrixData,Amat.GetRowMap(),Bmat.GetColMap(),@MatlabApply);
        elseif isa(Bmat,'MultiVector')
          Cmat = MultiVector( [], [], matvec(Amat, Bmat.GetVectors(), []) );
        else
          % Assuming that Bmat is a plain old Matlab matrix or vector
          matrixData = matvec(Amat, Bmat, []);
          if issparse(Bmat)
            matrixData = sparse(matrixData);
          end
          Cmat = matrixData;
          % TODO Convert Cmat back to an Operator.
          % For the moment, this break MueMat...
          % ColMap = Map(size(Bmat,2), 1);
          % Cmat = Operator(matrixData, Amat.GetRowMap(), ColMap, @MatlabApply);
        end

      end
    end %operator *

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function  [Cmat] = Apply(this,Bmat,Subset)
      %APPLY Matrix-vector multiplication.
      % Allows for multiplication that uses a subset of unknowns.
      %
      %   SYNTAX   Cmat = obj.Apply(Bmat, Subset);
      %
      %     Bmat   - an Operator
      %     Subset - Subset of indices on which to operate (optional, default=1:end)
      %     Cmat   - an Operator
      %
      if ~varexist('Subset') || isempty(Subset)
        Cmat = this * Bmat;
      else
        matvec = this.applyFcnHandle_;  %get the function handle
        if isa(Bmat,'Operator')
          if ~strcmp(func2str(Bmat.GetApply()),func2str(@MatlabApply)),
             error('Operator.Apply(): Bmat must be a MatlabApply matrix.');
          end
          matrixData = sparse(matvec(this, Bmat.matrixData_, Subset));
          Cmat = Operator(matrixData,this.GetRowMap(),Bmat.GetColMap(),@MatlabApply);
        else
          % Assuming that Bmat is a plain old Matlab matrix or vector
          Cmat = matvec(this, Bmat, Subset);
          if issparse(Bmat)
            Cmat = sparse(Cmat);
          end
        end
      end

    end %operator *

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [Atrans] = ctranspose(Amat)
      %CTRANSPOSE Overload of the unary transpose operator '.
      %
      %   SYNTAX   Atrans = obj.ctranspose(Amat);
      %
      %     Amat   - an Operator
      %     Atrans - A' (Operator)

      %FIXME: we only support 1 view. We remove everything else...

      % VERSION 1
      % TODO: deals with diagonals, invdata, other views...
      % Atrans = Amat;
      % Atrans = Atrans.SetColMap(Amat.GetRowMap());
      % Atrans = Atrans.SetRowMap(Amat.GetColMap());
      % Atrans.matrixData_ = Amat.matrixData_';

      % VERSION 2
      % The following version of ctranspose did not keep other
      % views, diagonals etc. but create a cleaner object:
      Atrans = Operator(Amat.GetMatrixData()', Amat.GetColMap(), Amat.GetRowMap(),@MatlabApply, Amat.GetLabel());
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % % NEVER USED IN MUEMAT. SEE SINGLE VIEW OPERATOR INSTEAD
    %     function  [Cmat] = MatMatSolve(Amat,Bmat)
    %       %MATMATSOLVE Operator solve.
    %       %
    %       %   SYNTAX   Cmat = MatMatSolve(Amat, Bmat);
    %       %
    %       %     Amat - an Operator
    %       %     Bmat - an Operator
    %       %     Cmat - an Operator
    %     end % Operator solve.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [ToF] = eq(Amat,Bmat)
      %EQ Overload of the binary logical operator == (NOT FINISHED).
      %
      %   SYNTAX   TorF = (Amat == Bmat);
      %
      %     Amat - an Operator
      %     Bmat - an Operator
      %     TorF - (boolean)
      if Amat.matrixData_ ~= Bmat.matrixData_
        ToF = false;
        % FIXME Map class needs == operator
      %elseif Amat.rowmap_ ~= Bmat.rowmap_
      %  ToF = false;
      %elseif Amat.colmap_ ~= Bmat.colmap_
      %  ToF = false;
      %elseif ~strcmp(func2str(Amat.apply_),func2str(Bmat.apply_)),
      elseif ~strcmp(func2str(Amat.GetApply),func2str(Bmat.GetApply)),
        ToF = false;
      else
        ToF = true;
      end
    end %Operator ==

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Print(this)
      %PRINT Print the Operator.
      %
      %   SYNTAX   obj.Print();
      %
      A = this.matrixData_;
      Nrows = size(A,1);   Ncols = size(A,2);
      fprintf('%s :\n',this.label_);
      for i=1:Nrows
        for j=1:Ncols
            if full(A(i,j)~=0)
                fprintf('%4.1f ',full(A(i,j)));
            else
                fprintf('  .  ');
            end
        end
        fprintf('\n');
      end
    end %Print

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function SetLabel(this,label)
      %SETLABEL Associate a name with Operator.
      %
      %   SYNTAX   obj.SetLabel(label);
      %
      %     label - name of Operator (string,optional,default='unamed')
      this.label_ = label;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [label] = GetLabel(this)
      %GETLABEL Get the Operator's name (string format).
      %
      %   SYNTAX   label = obj.GetLabel();
      %
      label = this.label_;
    end

% DEPRECATED
%     function [svOp] = GetSVOperator(this)
%       %GETSVOPERATOR Get the currect view's SingleView Operator
%       %
%       %   SYNTAX   svOp = obj.GetSVOperator();
%       %
%       %     svOp - a Single-view Operator (OperatorView)
%       svOp = this.viewTable_.GetOperator(this.currentView_);
%     end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %% Methods to manage Operator views
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [rowmap] = GetRowMap(this) %TODO: arg view
      %GETROWMAP Get the currect view's row map
      %
      %   SYNTAX   rowmap = obj.GetRowMap();
      %
      %     rowmap - a row map (Map)

      %FIXME: optional view argument
      if isempty(this.rowmapPtr_)
        error('Empty rowmap.');
      end
      rowmap = this.rowmapPtr_;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [colmap] = GetColMap(this,viewLabel)
      %GETCOLMAP Get the currect view's column map
      %
      %   SYNTAX   colmap = obj.GetColMap(view);
      %
      %     viewLabel - an operator view (optional, default=defaultView)
      %     colmap - a column map (Map)
      if ~varexist('viewLabel'), viewLabel = this.currentView_; end
      if  varexist('viewLabel') && strcmp(viewLabel,'current'), viewLabel = this.currentView_; end
      if varexist('viewLabel') && strcmp(viewLabel,'default'), viewLabel = this.defaultView_; end

      colmap = this.operatorViewTable_(viewLabel).GetColMap();
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [applyFcn] = GetApply(this,viewLabel)
      %GETAPPLY Get the current view's apply (matvec) function handle.
      %
      %   SYNTAX   applyFcn = obj.GetApply(viewLabel);
      %
      %     viewLabel     - an operator view (optional, default=defaultView)
      %     applyFcn - matrix-vector multiplication function (handle)
      if ~varexist('viewLabel'), viewLabel = this.currentView_; end
      if  varexist('viewLabel') && strcmp(viewLabel,'current'), viewLabel = this.currentView_; end
      if varexist('viewLabel') && strcmp(viewLabel,'default'), viewLabel = this.defaultView_; end;

      if (strcmp(viewLabel,this.currentView_) == false) || isempty(this.applyFcnHandle_)
        this.applyFcnHandle_ = this.operatorViewTable_(viewLabel).GetApply();
      end
      applyFcn = this.applyFcnHandle_;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [data] = GetMatrixData(this)
      %GETMATRIXDATA Get the underlying Matlab matrix.
      %
      %   SYNTAX   data = obj.GetMatrixData(viewLabel);
      %
      %     data - matrix data (MATLAB matrix)

      data = this.matrixData_;
    end

    function [data] = GetPattern(this)
      %GETPATTERN Get the graph/pattern of the matrix.
      %
      %   SYNTAX   data = obj.GetPattern(viewLabel);
      %
      %     data - matrix data (MATLAB matrix)
      %
      % See also: SetPattern

      data = this.matrixPattern_;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [ApplyInv] = GetApplyInverse(this,viewLabel)
      %GETAPPLYINVERSE Get the currect view's matvec
      %
      %   SYNTAX   ApplyInv = obj.GetApplyInverse(viewLabel);
      %
      %     viewLabel     - an operator viewLabel (optional, default=defaultView)
      %     ApplyInv - apply inverse function (handle)
      if ~varexist('viewLabel'), viewLabel = this.currentView_; end
      if  varexist('viewLabel') && strcmp(viewLabel,'current'), viewLabel = this.currentView_; end
      if varexist('viewLabel') && strcmp(viewLabel,'default'), viewLabel = this.defaultView_; end

      ApplyInv = this.operatorViewTable_(viewLabel).GetApplyInverse();
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [bool] = HasDiagonal(this,viewLabel)
      if ~varexist('viewLabel'), viewLabel = this.currentView_; end
      if  varexist('viewLabel') && strcmp(viewLabel,'current'), viewLabel = this.currentView_; end
      if varexist('viewLabel') && strcmp(viewLabel,'default'), viewLabel = this.defaultView_; end

      diagonal = [];
      if this.diagonalViewTable_.isKey(viewLabel)
        diagonal = this.diagonalViewTable_(viewLabel);
      end

      bool = ~isempty(diagonal);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [diagonal, flag] = GetDiagonal(this, Collection, viewLabel) %flag
      %GETDIAGONAL Get the view's diagonal
      %
      % If the diagonal does not exist, it is created.
      %
      %     collection - (optional, only use if the diagonal doesn't exist)
      %     viewLabel     - an operator view (optional, default=defaultView)
      %     diagonal - diagonal matrix (single-view Operator)
      %     flag     - true if the diagonal was already stored, false if the diagonal have been built (optional)
      if ~varexist('Collection'), Collection = []; end
      if ~varexist('viewLabel'), viewLabel = this.currentView_; end
      if  varexist('viewLabel') && strcmp(viewLabel,'current'), viewLabel = this.currentView_; end
      if varexist('viewLabel') && strcmp(viewLabel,'default'), viewLabel = this.defaultView_; end

      flag = this.HasDiagonal(viewLabel);
      if flag
        diagonal = this.diagonalViewTable_(viewLabel);
        %TODO: test to check if when ~isempty(Collection),
        %Diag.Collection==Collecton (else, throw error)
      else
        oldView = this.SwitchToView(viewLabel); %TODO: add extra arg to Diagonal constructor?
        diagonal = Diagonal(this,Collection);
        %diagonal = Diagonal2(this,Collection);
        this.SwitchToView(oldView);        %TODO

        this.SetDiagonal(diagonal,viewLabel);
      end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function SetRowMap(this,rowmap,viewLabel)
      %SETROWMAP Associate rowmap with 'viewLabel'.
      % If 'viewLabel' is empty, this is assumed to be for the current view.
      %
      %   SYNTAX   obj.SetRowMap(rowmap, viewLabel);
      %
      %     rowmap - a row map (Map)
      %     viewLabel - an operator view (optional, default=defaultView)
      if ~varexist('viewLabel'), viewLabel = this.currentView_; end
      if  varexist('viewLabel') && strcmp(view,'current'), viewLabel = this.currentView_; end
      if varexist('viewLabel') && strcmp(viewLabel,'default'), viewLabel = this.defaultView_; end

      this.operatorViewTable_(viewLabel).SetRowMap(rowmap);

      % update rowmapPtr_
      % (if viewLabel==this.CurrentView(), rowmapPtr_ must be updated)
      this.rowmapPtr_ = this.operatorViewTable_(this.CurrentView()).GetRowMap();
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function SetColMap(this,colmap,viewLabel)
      %SETCOLMAP Associate column map with 'viewLabel'.
      % If 'viewLabel' is empty, this is assumed to be for the current view.
      %
      %   SYNTAX   obj.SetColMap(colmap, viewLabel);
      %
      %     colmap - a column map (Map)
      %     viewLabel   - an operator view (optional, default=defaultView)
      if ~varexist('viewLabel'), viewLabel = this.currentView_; end
      if  varexist('viewLabel') && strcmp(viewLabel,'current'), viewLabel = this.currentView_; end
      if varexist('viewLabel') && strcmp(viewLabel,'default'), viewLabel = this.defaultView_; end

      this.operatorViewTable_(viewLabel).SetColMap(colmap);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function SetMatVec(this,matvec,viewLabel)
      %SETMATVEC Associate matrix-vector multiply handle with 'viewLabel'.
      % If 'viewLabel' is empty, this is assumed to be for the current view.
      %
      %   SYNTAX   obj.SetMatVec(matvec, viewLabel);
      %
      %     matvec - matrix-vector multiplication function (handle)
      %     viewLabel   - an operator view (optional, default=defaultView)
      if ~varexist('viewLabel'), viewLabel = this.currentView_; end
      if  varexist('viewLabel') && strcmp(viewLabel,'current'), viewLabel = this.currentView_; end
      if varexist('viewLabel') && strcmp(viewLabel,'default'), viewLabel = this.defaultView_; end

      this.operatorViewTable_(viewLabel).SetApply(matvec);

      % update applyFcnHandle_
      % (if viewLabel ==this.CurrentView(), applyFcnHandle_ must be updated)
      this.applyFcnHandle_ = this.operatorViewTable_(this.CurrentView()).GetApply();
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function SetApplyInverse(this,ApplyInv,viewLabel)
      %SETAPPLYINVERSE Associate diagonal matrix solve function handle with 'viewLabel'.
      % If 'viewLabel' is empty, this is assumed to be for the current view.
      %
      %   SYNTAX   obj.SetApplyInverse(ApplyInv, viewLabel);
      %
      %     ApplyInv - apply inverse function (handle)
      %     viewLabel     - an operator view (optional, default=defaultView)
      if ~varexist('viewLabel'), viewLabel = this.currentView_; end
      if  varexist('viewLabel') && strcmp(viewLabel,'current'), viewLabel = this.currentView_; end
      if varexist('viewLabel') && strcmp(viewLabel,'default'), viewLabel = this.defaultView_; end

      this.operatorViewTable_(viewLabel).SetApplyInverse(ApplyInv);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function SetMatrixInvData(this,MatrixInvData,viewLabel)
      %SETMATRIXINVDATA Associate diagonal matrix inverse data with 'viewLabel'.
      % If 'viewLabel' is empty, this is assumed to be for the current view.
      %
      %   SYNTAX   obj.SetMatrixInvData(MatrixInvData, viewLabel);
      %
      %     MatrixInvData - inverse matrix data (MATLAB matrix)
      %     viewLabel          - an operator view (optional, default=defaultView)
      if ~varexist('viewLabel'), viewLabel = this.currentView_; end
      if  varexist('viewLabel') && strcmp(viewLabel,'current'), viewLabel = this.currentView_; end
      if varexist('viewLabel') && strcmp(viewLabel,'default'), viewLabel = this.defaultView_; end

      this.viewlTable_ = this.operatorViewTable_(viewLabel).SetMatrixInvData(MatrixInvData);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function SetDiagonal(this,diagonal,viewLabel)
      %SETDIAGONAL Associate diagonal with 'viewLabel'.
      % If 'viewLabel' is empty, this is assumed to be for the current viewLabel.
      %
      %   SYNTAX   obj.SetDiagonal(diagonal, viewLabel);
      %
      %     diagonal - diagonal matrix (single-view Operator)
      %     viewLabel     - an operator view (optional, default=defaultView)
      if ~varexist('viewLabel'), viewLabel = this.currentView_; end
      if  varexist('viewLabel') && strcmp(viewLabel,'current'), viewLabel = this.currentView_; end
      if varexist('viewLabel') && strcmp(viewLabel,'default'), viewLabel = this.defaultView_; end

      this.diagonalViewTable_(viewLabel) = diagonal;
    end

    function SetMatrixData(this,data)
      %SETMATRIXDATA Set the underlying Matlab matrix data.
      %
      %   SYNTAX   obj.SetMatrixData(data);
      %
      %     data - matrix data (MATLAB matrix)
      this.matrixData_ = data;
    end

    function SetPattern(this,data)
      %SETPATTERN Define the non-zero pattern of the matrix.
      %
      %   SYNTAX   obj.SetPattern(data);
      %
      %     data - matrix data (MATLAB matrix)
      %
      % The non-zero pattern of the matrix can be stored independently of the actual data.
      % It allows to distinguish between the zeros of the pattern and true
      % zeros inside of the pattern.
      % Ex: Ptent have true zero entries inside of his pattern when the nullspace has zero entries.
      % This is problematic when Ptent is used as the pattern in Emin.
      %
      % Note: the pattern is never set/updated by any Operator operations.

      this.matrixPattern_ = data;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function SetDinvALambda(this,lambda,viewLabel)
      %SETDINVALAMBDA Set eigenvalue of scaled matrix associated with current viewLabel.
      %
      %   SYNTAX   obj.SetDinvALambda(lambda, viewLabel);
      %
      %     lambda - eigenvalue of D^{-1}.A (double)
      %     viewLabel   - an operator view (optional, default=defaultView)
      if ~varexist('viewLabel'), viewLabel = this.currentView_; end
      if  varexist('viewLabel') && strcmp(viewLabel,'current'), viewLabel = this.currentView_; end
      if varexist('viewLabel') && strcmp(viewLabel,'default'), viewLabel = this.defaultView_; end

      diagonal = this.diagonalViewTable_(viewLabel);
      diagonal.SetDinvALambda(lambda);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [lambda] = GetDinvALambda(this,viewLabel)
      %GETDINVALAMBDA Get eigenvalue of scaled matrix associated with current view.
      %
      %   SYNTAX   lambda = obj.GetDinvALambda(viewLabel);
      %
      %     viewLabel   - an operator view (optional, default=defaultView)
      %     lambda - eigenvalue of D^{-1}.A (double)
      if ~varexist('viewLabel'), viewLabel = this.currentView_; end
      if  varexist('viewLabel') && strcmp(viewLabel,'current'), viewLabel = this.currentView_; end
      if  varexist('viewLabel') && strcmp(viewLabel,'default'), viewLabel = this.defaultView_; end

      diagonal = [];
      if this.diagonalViewTable_.isKey(viewLabel)
        diagonal = this.diagonalViewTable_(viewLabel);
      end

      if isempty(diagonal), lambda = [];
      else lambda = diagonal.GetDinvALambda(); end
      if isempty(lambda)
        % Compute lambda (note: diag is created if it doesn't exist

        n = this.GetRowMap().NDOFs();

        BlkDiag = this.GetDiagonal([], viewLabel);

        if isempty(BlkDiag.GetApplyInverse())
          BlkDiag.FactorBlkDiag();
        end

        lambda = MaxEigenvalue(@this.DinvAv, n, this, BlkDiag);
        this.SetDinvALambda(lambda, viewLabel);
      end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [v] = DinvAv(this, y, Amat, BlkDiag)
      % Internal function used for via matrix-vector product form of eigs
      applyFcn = Amat.GetApply();
      z = applyFcn(Amat, y, []);
      applyInvFcn = BlkDiag.GetApplyInverse();
      v = applyInvFcn(BlkDiag, z , []);
    end

  end %methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

end %Operator classdef
