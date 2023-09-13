classdef CoalesceDropFactory < VerboseObject
% Factory which constructs graphs using matrices
% These graphs are then
% typically passed to an aggregation or F/C point factory which will be used
% to coarsen the fine level space.
%
% Specifically, we take a matrix and drop entries using the
% PreDropFunc_(). The block matrix is then collapsed to a point matrix using
% AmalgFunc_() to define the values in the collapsed point matrix.  Finally,
% another round of dropping entries occuring using PostDropFunc_().
%
% PreDropFunc and PostDropFunc must take a matrix and return the dropped matrix.
% Function prototypes : function [matrix] = MyDropFunc(matrix, DropData);
%
% Note: The blocking associated with the amalgamation is assumed to be
%       defined within the matrix.
%

  properties (Access = private)
      PreDropFunc_     % User supplied function which takes the inital matrix A and drops entries in it
      PreDropData_     % User supplied data passed to PreDropFunc_()
      AmalgFunc_       % User supplied function which takes a dense block and replaces it by a single number
      AmalgData_       % User supplied data passed to AmalgFunc_()
      PostDropFunc_    % User supplied function which drops entries in Ahat where Ahat is the amalgamated matrix
      PostDropData_    % User supplied data passed to PostDropData_()
      OutputType_      % if 'Binary', all edge weights = 1.

      reUseGraph_ = false

      AName_ = 'A';
      AFilteredName_ = 'Afiltered';
   end
   methods
      function [this] = CoalesceDropFactory(arg)
         % Copy constructor
         if nargin == 1 && isa(arg, class(this)), this.Copy_(arg,[]); return; end
         %

         % constructor based on default values.
         this.OutputType_   = 'Binary';
         this.PreDropData_  = 0.;
         this.PostDropData_ = 0.;
         this.AmalgData_    = 0.;
         this.PreDropFunc_  = [];
         this.PostDropFunc_ = []; %@SAdrop
         this.AmalgFunc_    = @MatrixNorm;
      end

      % TODO: if ~varexist('PreDropData') => reinit PreDropData
      function SetPreDropSpecifications(this,PreDropFunc,PreDropData)
      % Assign predrop function used to decide matrix elements ignored before amalgamation occurs
         this.PreDropFunc_ = PreDropFunc;
         if varexist('PreDropData'), this.PreDropData_ = PreDropData; end
      end
      function SetPostDropSpecifications(this,PostDropFunc,PostDropData)
      % Assign postdrop function used to decide matrix elements ignored after amalgamation occurs
         this.PostDropFunc_ = PostDropFunc;
         if varexist('PostDropData'), this.PostDropData_     = PostDropData; end
      end
      function SetAmalgSpecifications(this,AmalgFunc,AmalgData)
      % Assign amalgamation function which maps block values into a scalar representing the block
         this.AmalgFunc_ = AmalgFunc;
         if varexist('AmalgData'), this.AmalgData_     = AmalgData; end
      end
      function SetNonBinary(this)
      % Indicates that constructed graph should include weighted edges
         this.OutputType_ = 'NotBinary';
      end
      function SetBinary(this)
      % Indicates that constructed graph should include nonweighted edges
         this.OutputType_ = 'Binary';
      end

      function [PreDropFunc, PreDropData] = GetPreDropSpecifications(this)
         PreDropFunc = this.PreDropFunc_;
         PreDropData = this.PreDropData_;
      end
      function [PostDropFunc, PostDropData] = GetPostDropSpecifications(this)
         PostDropFunc = this.PostDropFunc_;
         PostDropData = this.PostDropData_;
      end
      function [AmalgFunc, AmalgData] = GetAmalgSpecifications(this)
         AmalgFunc = this.AmalgFunc_;
         AmalgData = this.AmalgData_;
      end
      function [OutputType] = GetOutputType(this) %TODO: non symetric with 'Set' function
         OutputType = this.OutputType_;
      end

      function [ToF] = ReUseGraph(this, ToF)
        if varexist('ToF'),
          ToFold = this.reUseGraph_;
          this.reUseGraph_ = ToF;
          ToF = ToFold;
        else ToF = this.reUseGraph_;
        end
      end

      function SetAName(this, name)
        this.AName_ = name;
      end

      function [name]=GetAName(this)
        name=this.AName_;
      end

      function SetAFilteredName(this, name)
        this.AFilteredName_ = name;
      end
      function [name]=GetAFilteredName(this)
        name=this.AFilteredName_;
      end

      function SetNeeds(this, CurrentLevel)
          % Obtain any cross factory specifications
          %if ~this.ReUseGraph()

          %end
      end

      function flag = Build(this,CurrentLevel,Specs)
        if strcmp(class(CurrentLevel),'Operator'), fprintf('MueMat:Deprecated %s','Use Build_() instead of Build()'); flag = this.Build_(CurrentLevel); dbstack(1); keyboard; return; end

        if CurrentLevel.IsAvailable('Graph',this),
            fprintf('reuse graph\n');
            return;
        end; % todo: make use of factories
%         if this.ReUseGraph()
%           if ~CurrentLevel.IsAvailable('Graph')
%             error('ReUseGraph: Aggregates not available. Use CurrentLevel.Print() for debug.');
%           else
%             return; % do nothing and return directly
%           end
%         end

        if CurrentLevel.IsAvailable('AuxMatrix') && ~strcmp(this.AName_,'AuxMatrix')
          dbstack
          keyboard
          this.AName_ = 'AuxMatrix';
          warning('MueMat:Deprecated','Please add a call to CoalesceDropFact.SetAName(''AuxMatrix'') in your driver. Right now, CoalesceDropFact always use AuxMatrix when it is available but it will not be the case on future version of MueMat');
        end

	Temprowmap = CurrentLevel.Get(this.AName_).GetRowMap();
        if(this.AmalgData_==0.0 & isempty(this.PreDropFunc_) & isempty(this.PostDropFunc_) & Temprowmap.MaxBlkSize==1)
	  if(isempty(CurrentLevel.Get(this.AFilteredName_))==true)
	    Afiltered=CurrentLevel.Get(this.AName_);
            Graph=Afiltered;
	  else
	    Afiltered=CurrentLevel.Get(this.AFilteredName_);
            Graph=Afiltered;
          end
	else
            fprintf('RST Warning: I had to comment this MueMat code out because I could not get it to work!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
% I couldn't get this to work right ... so I'm commenting it out? RST
%	  if(isempty(CurrentLevel.Get(this.AFilteredName_))==true)
	    [Graph,Afiltered]=this.Build_(CurrentLevel.Get(this.AName_),CurrentLevel);
%	  else
%	    [Graph,Afiltered]=this.Build_(CurrentLevel.Get(this.AFilteredName_),CurrentLevel);
%          end
	end

        % Notes: if you use an AuxMatrix that is already a nodal
        % matrix, then
        % - you can use CoalesceFact to do dropping on it if it has
        % not been done on it yet.
        % - Afiltered is also a nodal matrix and cannot be used
        % directly as a filteredA for prolongator smoother

        if ~isempty(Afiltered)
          Afiltered = Operator(Afiltered, CurrentLevel.Get(this.AName_).GetRowMap(),CurrentLevel.Get(this.AName_).GetColMap(), @MatlabApply);
          CurrentLevel.Set(this.AFilteredName_, Afiltered);
        end

        CurrentLevel.Set('Graph', Graph, this);
      end

      function [Graph,Afiltered] = Build_(this,Amat, CurrentLevel)
      % Build a graph from the matrix based on factory specifications
      % A block matrix is collapsed to a point matrix and small nonzeros
      % are dropped when defining the graph.

        if ~varexist('CurrentLevel'), CurrentLevel=[]; end

        [Graph,Afiltered]=this.AmalgamateAndOrDrop(Amat, CurrentLevel,this.PreDropFunc_,  ...
                                  this.PreDropData_, ...
                                  this.AmalgFunc_,   this.AmalgData_, ...
                                  this.PostDropFunc_,this.PostDropData_, ...
                                  this.OutputType_);

      end
   end
   methods (Access = protected)
     function [newA,Afiltered] =AmalgamateAndOrDrop(this, A, CurrentLevel, PreDropFunc, PreDropData, ...
                           AmalgFunc, AmalgData, PostDropFunc, PostDropData,...
                           OutputType)
      % Drop entries from matrix, amalgamate, and then drop more entries.
      % Please note that any of the above phases can be omitted by simply mapping
      % in a null object for the function pointer.
      %
      %
      %          Input:
      %             A           :Block matrix to be amalgamated
      %
      %             PreDropFunc :User supplied function which
      %                          determine whether A_ij is dropped
      %
      %             PreDropData :Data passed to PreDropFunc
      %
      %             AmalgFunc   :User supplied function which takes
      %                          a block (a dense matrix) and replaces
      %                          it by a single number
      %
      %             AmalgData   :Data passed to AmalgFunc
      %
      %             PostDropFunc:User supplied function which
      %                          determine whether A_ij is dropped
      %                          where A is amalgamated matrix.
      %
      %             PostDropData:Data passed to PreDropFunc
         Afiltered = [];
         Arowmap = A.GetRowMap();
         Acolmap = A.GetColMap();
         Adata   = A.GetMatrixData();

         % Initialization of newA (output matrix) members
         if ~isempty(AmalgFunc)
            newRowMap           = Map(Arowmap.NNodes(), 1);
            newColMap           = Map(Acolmap.NNodes(), 1);
         else
            newRowMap = Arowmap;
            newColMap = Acolmap;
         end

         matrix = A.GetMatrixData();
         original = matrix;

         % Predropping
         if ~isempty(PreDropFunc)
           matrix = PreDropFunc(matrix,CurrentLevel,PreDropData);
           Afiltered = matrix; % Afiltered is an output argument of the method
           % The user wants to create a 2nd filtered matrix. The one
           % already created is used as the matrix graph for things like
           % aggregation. The 2nd filtered matrix is used for things like
           % prolongator smoothing (or emin).
           if isfield(PreDropData,'Afiltered')
             twidget = PreDropData.Afiltered;
             TheFunc = PreDropFunc;
             if isfield(twidget,'Func')
                TheFunc = twidget.Func;
             end
             Afiltered = TheFunc(original,CurrentLevel,twidget);
             % Just in case an AuxMatrix was used for dropping, 
             % we'll always go back a grab the original A from CurrentLevel
             % and .* this with spones(Afiltered)

             ActualA  = CurrentLevel.Get('A').GetMatrixData;
             Afiltered= ActualA .* spones(Afiltered);
             %
             % If requested, make sure that Rowsum(Afiltered) = RowSum(A)
             % 
             if isfield(twidget,'FixRowSums')
                nnn     = size(Afiltered,1);
                GoodRowSums = ActualA*ones(nnn,1);
                BadRowSums  = Afiltered*ones(nnn,1);
                % do it like this for rounding errors
                Afiltered = Afiltered - spdiags(BadRowSums,0,nnn,nnn);
                Afiltered = Afiltered + spdiags(GoodRowSums,0,nnn,nnn);
             end
           end
         end



         % Matrix format conversion
         [Arows, Acols, Avals] = find(matrix);
         NN = length(Arows);
         clear matrix;

         % Amalgamation
         if ~isempty(AmalgFunc)
            BlkdVals = zeros(NN,1);
            if Arowmap.HasConstBlkSize()
               RowDim = Arowmap.ConstBlkSize();
               BlkdRows = ceil(Arows/RowDim - .001);
               Arows = Arows - RowDim*(BlkdRows-1);
            elseif Arowmap.HasVariableBlkSize()
               VarBlkRPtr = Arowmap.Vptr();
               BlkRIndex = zeros(Arowmap.NDOFs(),1);
               j = VarBlkRPtr(1);
               for i=1:Arowmap.NDOFs()
                  while i >= VarBlkRPtr(j), j = j+1; end
                  BlkRIndex(i) = j-1;
               end
               BlkdRows = zeros(NN,1);
               for i=1:NN
                  BlkdRows(i) = BlkRIndex(Arows(i));
                  Arows(i) = Arows(i) - VarBlkRPtr(BlkdRows(i))+1;
               end
            else
               fprintf('AmalgamateAndOrDrop:: Unknown matrix type\n');
            end
            if Acolmap.HasConstBlkSize()
               ColDim = Acolmap.ConstBlkSize();
               BlkdCols = ceil(Acols/ColDim - .001);
               Acols = Acols - ColDim*(BlkdCols-1);
            elseif Acolmap.HasVariableBlkSize()
               VarBlkCPtr = Acolmap.Vptr();
               BlkCIndex  = zeros(Acolmap.NDOFs(),1);
               j = VarBlkCPtr(1);
               for i=1:Acolmap.NDOFs()
                  while i >= VarBlkCPtr(j), j = j+1; end
                  BlkCIndex(i) = j-1;
               end
               BlkdCols = zeros(NN,1);
               for i=1:NN
                  BlkdCols(i) = BlkCIndex(Acols(i));
                  Acols(i) = Acols(i) - VarBlkCPtr(BlkdCols(i))+1;
               end
            else
               fprintf('AmalgamateAndOrDrop:: Unknown matrix type\n');
            end

               %
               % Sort all the lists so that the ordering correspond to blocks.
               % In particular, if kth entry corresponds to the (i,j)th block,
               % we want the k+1 entry to correspond to (ii,jj)th block where
               % ii >= i and jj >= j. Done by first sorting on block rows
               % and then within a block row sorting on block columns.
               %
               [BlkdRows,iii] = sort(BlkdRows);
               BlkdCols = BlkdCols(iii); Arows   = Arows(iii);
               Acols    = Acols(iii);    Avals   = Avals(iii);

                StartRow = 1; EndRow=StartRow;
                while EndRow < NN,
                   while (EndRow<NN) && ...
                         (BlkdRows(EndRow)==BlkdRows(EndRow+1)),
                      EndRow = EndRow + 1;
                   end
                   [ttemp,jjj] = sort(BlkdCols(StartRow:EndRow));
                   BlkdCols(StartRow:EndRow) = ttemp;
                   ttemp = Arows(StartRow:EndRow); Arows(StartRow:EndRow)= ttemp(jjj);
                   ttemp = Acols(StartRow:EndRow); Acols(StartRow:EndRow)= ttemp(jjj);
                   ttemp = Avals(StartRow:EndRow); Avals(StartRow:EndRow)= ttemp(jjj);
                   StartRow = EndRow+1; EndRow = StartRow;
                end

                %
                % Now, go through the sorted lists, build block matrix for
                % each block entry, invoke the amalgmate function
                %
                j = 1;
                RowMapHasVariableBlkSize = Arowmap.HasVariableBlkSize();
                ColMapHasVariableBlkSize = Acolmap.HasVariableBlkSize();

                for BlkRow=1:Arowmap.NNodes()
                    while  (j <= NN) && (BlkdRows(j) <= BlkRow),
                       BlkCol = BlkdCols(j);
                       if RowMapHasVariableBlkSize,
                          RowDim = VarBlkRPtr(BlkRow+1)-VarBlkRPtr(BlkRow);
                       end
                       if ColMapHasVariableBlkSize,
                          ColDim = VarBlkCPtr(BlkCol+1)-VarBlkCPtr(BlkCol);
                       end
                       matrix = zeros(RowDim,ColDim);

                       while (j <= NN ) && (BlkdCols(j) == BlkCol) && (BlkdRows(j) == BlkRow),
                          matrix( Arows(j), Acols(j)) = Avals(j);
                          j = j+1;
                       end
                       BlkdVals(j-1) = AmalgFunc(matrix,AmalgData);
                    end
                end
                SubInds = find(BlkdVals ~= 0);
                Arows = BlkdRows(SubInds);
                Acols = BlkdCols(SubInds);
                Avals = BlkdVals(SubInds);
                %NN = length(Arows); % unused
         end

         % Matrix format conversion
         matrix = sparse(Arows,Acols,Avals,newRowMap.NDOFs(), newColMap.NDOFs());
         clear Arows Acols Avals NN;

         % PostDropping
         if ~isempty(PostDropFunc),
           matrix = PostDropFunc(matrix,CurrentLevel,PostDropData);
         end

         %
         if strcmp(OutputType,'Binary'),
           newMatrixData = spones(matrix);
         else
           newMatrixData = matrix;
         end

         newA = Operator(newMatrixData,newRowMap, newColMap,@MatlabApply);

      end %AmalgamateAndOrDrop()

   end %Methods

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

end %classdef CoalesceDropFactory

