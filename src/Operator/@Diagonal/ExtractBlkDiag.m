%
% Extract the block diagonal of Operator (if it has not already been done).
%
% Case 1: Operator has constant block size
%
%      The block diagonal matrix is stored as a dense array
%      of dimension NN x NN*RowMap.NNodes where NN is the constant block size.
%      For example,
%           BlkDiag.MatrixData:
%                             [   Blk 1  ]
%                             [   Blk 2  ]
%                                  ...
%                             [   Blk NN ]
%
% Case 2: Operator has a variable block size
%
%      This could still be done as a dense matrix where the width is
%      equal to the largest block or as a sparse matrix where the width
%      is equal to the largest block. The sparse way would be useful
%      if we are doing block Gauss-Seidel and the blocks are not the same
%      as how A is blocked. In this case, we might end up with some pretty
%      big block diagonal matrix and I want these to be sparse. I've kind
%      of implemented both here using the hardwired variable SparseWay.

function ExtractBlkDiag(BlkDiag,Amat)  % BlkDiag == this

   SparseWay=1;   % indicates dense or sparse storage of a variable block matrix
   % TODO : SparseWay=0 is not working (MueMat crashs outside of
   % this function bcse such storage for BlkDiag is not supported
   % elsewhere in the code). Extraction of BlkDiag must also be
   % optimized for SparseWay=0.

   % Diagonal already computed. To recompute block diagonal, one
   % must first delete it before invoking this function.
   %

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

      DenseMat = zeros(Ntotal,ConstBlkSize);
      first = 1;
      for i=1:NBlks
         last  = i*ConstBlkSize;
         [rows,cols,vals] = find(MatrixData(:,first:last));
         for j=1:length(rows)
            if (rows(j) >= first) && (rows(j) <= last)
               DenseMat(rows(j),cols(j)) = vals(j);
            end
         end
         first = last+1;
      end

      if ConstBlkSize == 1, MatrixData = sparse(DenseMat);
      else                 MatrixData = DenseMat;        end
      ApplyFunc       = @DBlkApply;
      ApplyInvFunc    = @DinvBlkApply;
   elseif RowMap.HasVariableBlkSize()
      VarBlkPtr = RowMap.Vptr();
      Ntotal = VarBlkPtr(NBlks+1)-1;
      Bsizes = zeros(NBlks,1);
      for i=1:NBlks
         Bsizes(i) = VarBlkPtr(i+1)-VarBlkPtr(i);
      end
      MaxBlkSize = max(Bsizes);
      SparseMat  = spalloc(Ntotal,MaxBlkSize,Ntotal*MaxBlkSize);
      first = VarBlkPtr(1);
      if SparseWay,
         for i=1:NBlks
            last  = VarBlkPtr(i+1)-1;
            [rows,cols,vals] = find(MatrixData(:,first:last));
            for j=1:length(rows)
              if (rows(j) >= first) && (rows(j) <= last)
                SparseMat(rows(j),cols(j)) = vals(j);
              end
            end
            first = last+1;
         end
         MatrixData  = SparseMat;
      else
         DenseMat = zeros(Ntotal,MaxBlkSize);
         first = VarBlkPtr(1);
         count = 1;
         for i=1:-NBlks
            last  = VarBlkPtr(i+1)-1;
            inds = (first:last);
            DenseMat(count:count+Bsizes(i)-1,1:Bsizes(i))= MatrixData(inds,inds);%slow
            first = last+1;
            count = count + Bsizes(i);
         end
         MatrixData  = DenseMat;
      end
      ApplyFunc       = @DBlkApply;
      ApplyInvFunc    = @DinvBlkApply;
    else
      fprintf('ExtractBlkDiag:: Unrecognized matrix type\n');
      keyboard;
   end
%   BlkDiag = OperatorView(RowMap,ColMap,ApplyFunc);
   BlkDiag.SetRowMap(RowMap);
   BlkDiag.SetColMap(ColMap);
   BlkDiag.SetApply(ApplyFunc);

   BlkDiag.SetMatrixData(MatrixData); %TODO: MatrixData is used for Amat and then for the Diagonal. Another name should be used (DiagonalData)
   %BlkDiag.SetApplyInverse(ApplyInvFunc);
   %FactorBlkDiag(BlkDiag);
%  Amat.SetDiagonal(BlkDiag);
