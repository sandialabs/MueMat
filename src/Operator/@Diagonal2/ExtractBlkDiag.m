%
% Extract the block diagonal of the Operator.
%

function ExtractBlkDiag(BlkDiag, Amat) % BlkDiag == this

   RowMap = Amat.GetRowMap();
   ColMap = Amat.GetColMap();
   if xor(HasConstBlkSize(ColMap),HasConstBlkSize(RowMap)),
      fprintf('ExtractBlkDiag:: Mismatched Row/Col Block types\n');
      keyboard;
   end

   MatrixData = Amat.GetMatrixData();
   NBlks = RowMap.NNodes();
   if HasConstBlkSize(RowMap),
      ConstBlkSize = RowMap.ConstBlkSize();
      Ntotal = NBlks*ConstBlkSize;

      % %Simplest way to extract the block diagonal
      % DiagonalData = spalloc(size(MatrixData, 1), size(MatrixData, 2), NBlks*ConstBlkSize*ConstBlkSize);
      % for i=1:NBlks
      %  range = (i-1)*ConstBlkSize + [1:ConstBlkSize];
      %  DiagonalData(range, range) = MatrixData(range, range);
      % end

      % % Faster way to extract the block diagonal
      % Creating a filter (DiagBlockPattern) to keep only the block diag values of MatrixData
      a = ones(ConstBlkSize, ConstBlkSize, NBlks);
      c = num2cell(a,[1 2]);
      sparsec = cellfun(@sparse, c, 'UniformOutput', false); % convert to sparse matrix before using blkdiag
      DiagBlockPattern = blkdiag(sparsec{:});
      DiagonalData = MatrixData .* DiagBlockPattern;

      %
      ApplyFunc       = @DBlkApply;
      ApplyInvFunc    = @DinvBlkApply;
   elseif RowMap.HasVariableBlkSize()
      fprintf('This case is not handle (yet) by Diagonal2. Use the class Diagonal instead\n');
      keyboard
   else
      fprintf('ExtractBlkDiag:: Unrecognized matrix type\n');
      keyboard;
   end

   BlkDiag.SetRowMap(RowMap);
   BlkDiag.SetColMap(ColMap);
   BlkDiag.SetApply(ApplyFunc);

   BlkDiag.SetMatrixData(DiagonalData);

end