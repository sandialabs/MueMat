%
% Extract block diagonal of Operator where blocks are given by Collection.
% In particular, the total number of blocks is Collection.NSubsets and the
% ith block of the block diagonal is made up of blocks in Operator
% corresponding to Collection.Subset(i). We still have the possibility
% to store this matrix in either a sparse or a dense format. Here,
% it might more sense to store it as sparse given that the diagonal matrices
% might vary a fair amount and might be pretty large.
%
%
%
function [BlkDiag] = ExtractNonContigBlkDiag(BlkDiag,Amat,Collection)

   MatrixData= Amat.GetMatrixData();
   VarBlkPtr = zeros(Collection.NSubsets+1,1);
   Bsizes    = zeros(Collection.NSubsets,1);
   VarBlkPtr(1) = 1;
   RowMap = Amat.GetRowMap();
   ColMap = Amat.GetColMap();
   if RowMap.HasConstBlkSize(),
      ConstBlkSize = RowMap.ConstBlkSize();
      NRegBlksPerBigBlk = zeros(Collection.NSubsets,1);
      for i=1:Collection.NSubsets
         NRegBlksPerBigBlk(i) = length(Collection.Subsets(i).BlkList);
      end
      NTotalRows = ConstBlkSize*sum(NRegBlksPerBigBlk);
      MaxBlkSize = ConstBlkSize*max(NRegBlksPerBigBlk);
   elseif RowMap.HasVariableBlkSize()
      NTotalRows = 0;
      MaxBlkSize = 0;
      for i=1:Collection.NSubsets
         cols = Node2DOF(Collection.Subsets(i).BlkList, ColMap);
         NTotalRows = NTotalRows + length(cols);
         if length(cols) > MaxBlkSize, MaxBlkSize = length(cols); end
      end
   else
      fprintf('ExtractNonContigBlkDiag:: Unrecognized matrix type\n');
      keyboard;
   end

   first = 1;
   numNz = 0;
   fprintf('extracting block diagonal\n');
   tic
   offset = zeros(Collection.NSubsets,1);
   offset(1) = 1;
   % Save each block diagonal submatrix in a separate variable called IJVk, where k is in 1..#blocks.
   % In this way, the exact number of nonzeros in the block diagonal is known, and memory for
   % it can be allocated all at once.  After the loop, put all the IJVk into a single
   % structure IJV, and convert it to a sparse matrix.  This is about 5x faster than growing
   % a single IJV structure in the loop.  The downside is that at its highpoint, this needs
   % twice the memory of storing just the block diagonal -- oh well.
   maxNnz=0; maxIdx = 1; maxRows = 0;
   for i=1:Collection.NSubsets
      cols = Node2DOF(Collection.Subsets(i).BlkList, ColMap);
      last  = first + length(cols)-1;
      [I,J,V] = find(MatrixData(cols,cols));
      eval(['IJV' num2str(i) '= [I+first-1 J V];']);
      numNz = numNz + length(I);
      if length(I) > maxNnz, maxNnz = length(I); maxIdx = i; maxRows = length(cols); end
      offset(i+1) = offset(i) + length(I);
      VarBlkPtr(i+1) = VarBlkPtr(i) + (last-first+1);
      Bsizes(i) = length(cols);
      first = last+1;
   end
   IJV = zeros(numNz,3);
   offset = 1;
   for i=1:Collection.NSubsets
      eval(['tt=IJV' num2str(i) ';']);
      if length(tt) > 0
        IJV(offset:offset+length(tt)-1,:) = tt;
      end
      offset = offset + length(tt);
      eval(['clear IJV' num2str(i) ';']);
   end


   SparseMat = spconvert(IJV);
   clear IJV
   toc
   fprintf('actual #nonzeros for block diagonal = %g\n',offset-1);
   fprintf('max #nnz in a block = %d (block %d, %d x %d)\n',maxNnz,maxIdx,maxRows,maxRows);
   NewMap = Map(Collection.NSubsets,Bsizes);

   % BlkDiag = OperatorView(NewMap,NewMap,@DBlkApply);
   BlkDiag.SetRowMap(NewMap);
   BlkDiag.SetColMap(NewMap);
   BlkDiag.SetApply(@DBlkApply);
   BlkDiag.SetMatrixData(SparseMat);
   %BlkDiag.SetApplyInverse(@DinvBlkApply);
%  FactorBlkDiag(BlkDiag);

