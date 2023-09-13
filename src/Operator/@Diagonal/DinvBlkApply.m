%
% This file is generated from the template DxBlkApply.tmpl.m
%
% This function is a matrix-solve with a block diagonal matrix where
% the blocks have a constant size (HasConstBlkSize) or have
% different sizes (HasVariableBlkSize).
%
% There is also a possibility of performing the matrix-solve only on a
% subset of block rows which can be contiguous or scattered
% (see CreateDOFSubset.m).
% The block diagonal matrix is either already factored or factors will
% be produced (via FactorBlkDiag.m).
% The block diagonal matrix LU is stored in a special way (i.e. different
% from a general matrix operator) in order to gain some efficiency (see
% ExtractBlkDiag.m and ExtractNonContigBlkDiag.m for details on the
% format).
%
% When a matrix-solve is performed only involving a subset
% of rows, there are two possibilities for the form of the input vector:
%     Case 1:  Standard case where the input vector is the entire
%              input vector for the matrix even though only a subset
%              of this data is used in performing the subset
%              matrix-solve.
%
%     Case 2:  Since this is a block diagonal matrix, we know that the
%              subset of rows for the output is in fact also the subset of
%              columns needed for the input. Thus, it is assumed that
%              we are only given this input data and it is stored contiguously.
%
% Case 1 and Case 2 are identified by the input parameter VecType.
%              (case 1: VecType=GLOBAL, case 2: VecType=LOCAL)
%
% Note: The output vector length is the number of rows in Subset stored
%       contiguously (regardless of whether the input corresponds to
%       Case 1 or Case 2).
%
% Optimization:
% --------------------------------------------------------------------------
%  We don't know the number of nz of OutVec in advance.
%  Submatrices of OutVec are stored separately and assembled at the end.
%  To store submatrices, we use a IJVCellMatrix object.
%
%  For (ConstBlkSize == 1), the simplier subfunction DinvNoBlkApply is called.
%
% Loop Variables:
% --------------------------------------------------------------------------
% * fGlobalRow lGlobalRow:
%   = first and last row for LU Blk (entire matrix).
% * fLocalRow lLocalRow:
%   = first and last row for OutVec Blk. Rows of OutVec are stored
%    contiguously and the output vector length is the number of rows in Subset
% * row indices of Vec depends on VecType (case 1=GLOBAL, case 2=LOCAL).
%
% The only difference between ConstBlkSize and VariableBlkSize is
% how are set the loop variables.

function [OutVec] = DinvBlkApply(A, Vec, Subset, VecType)
   GLOBAL = 1; % defined by mue_include. Inclusion of mue_include is slow
   if ~varexist('Subset'), Subset = []; end
   if ~varexist('VecType'), VecType = GLOBAL; end

   RowMap = A.GetRowMap();

   if RowMap.HasConstBlkSize() && (RowMap.ConstBlkSize() == 1)
     OutVec = DinvNoBlkApply(A, Vec, Subset, VecType); %optim
   else

     BlkList = SubsetBlkList(A, Subset);

     % input (LU)
     A.FactorBlkDiag(); % This do nothing if the factorization is already computed
     widget = A.GetDiagonalInvData();

     % output (OutVec)
     nRow = SubsetNRow(A, Subset);
     %TEST: OutVecTEST = zeros(nRow,size(Vec,2));
     OutVecIJV = IJVCellMatrix(length(BlkList)); % Unassembled matrix

     % computation
     if RowMap.HasConstBlkSize()
       BlkSize = RowMap.ConstBlkSize();

       fLocalRow = 1;
       for i=BlkList;
         fGlobalRow = (i-1)*BlkSize+1;
         lGlobalRow = fGlobalRow + BlkSize-1;
         lLocalRow = fLocalRow + BlkSize-1;

         %OutVecPart = DinvSubFuncBlkApply(widget, Vec, [fGlobalRow fLocalRow], [lGlobalRow lLocalRow], BlkSize, VecType);
         OutVecPart = DinvSubFuncBlkApply(widget, Vec, [fGlobalRow fLocalRow], [lGlobalRow lLocalRow], i, BlkSize, VecType);
         %TEST: OutVecTEST(fLocalRow:lLocalRow,:) = OutVecPart;
         [I,J,V] = find(OutVecPart);
         I = I + fLocalRow - 1; % shift == OutVec(fLocalRow:lLocalRow,:) = DSubFuncBlkApply(...)
         OutVecIJV.Add(i,I,J,V);

         fLocalRow = lLocalRow+1;
       end

     else
       VarBlkPtr = RowMap.Vptr();

       fLocalRow = 1;
       for i=BlkList;
         fGlobalRow = VarBlkPtr(i);
         lGlobalRow = VarBlkPtr(i+1)-1;
         BlkSize = lGlobalRow-fGlobalRow+1;
         lLocalRow = fLocalRow + BlkSize-1;
         %OutVecPart = DinvSubFuncBlkApply(widget, Vec, [fGlobalRow fLocalRow], [lGlobalRow lLocalRow], BlkSize, VecType);
         OutVecPart = DinvSubFuncBlkApply(widget, Vec, [fGlobalRow fLocalRow], [lGlobalRow lLocalRow], i, BlkSize, VecType);
         %TEST: OutVecTEST(fLocalRow:lLocalRow,:) = OutVecPart;
         [I,J,V] = find(OutVecPart);
         I = I + fLocalRow - 1; % shift == OutVec(fLocalRow:lLocalRow,:) = DSubFuncBlkApply(...)
         OutVecIJV.Add(i,I,J,V);

         fLocalRow = lLocalRow+1;
       end

     end

     [II, JJ, VV] = OutVecIJV.GetIJV();
     OutVec = sparse(II,JJ,VV,nRow,size(Vec,2)); % note: it's important to give the matrix size (if empty row/col)
     %TEST: if (nnz(OutVec-OutVecTEST) ~= 0), ST=dbstack; error('%s:: test error.\n', ST(1).name); end

   end

end


% Main kernel (to compute OutVec)
%
% fRow(GLOBAL) lRow(GLOBAL):
%  = first and last row for LU Blk (entire matrix).
% fRow(LOCAL) lRow(LOCAL):
%  = first and last row for OutVec Blk. Rows of OutVec are stored
%    contiguously and the output vector length is the number of rows in Subset
% row indices of Vec depends on VecType (case 1=GLOBAL, case 2=LOCAL).
%function [OutVecPart] = DinvSubFuncBlkApply(widget, Vec, fRow, lRow, BlkSize, VecType)
function [OutVecPart] = DinvSubFuncBlkApply(widget, Vec, fRow, lRow, BlkInd, BlkSize, VecType)
   GLOBAL = 1; % defined by mue_include. Inclusion of mue_include is slow

   Perm = widget.Perm; Qmat = widget.Qmat;
   if isfield(widget,'L') && isfield(widget,'U')
     L = widget.L{BlkInd};
     U = widget.U{BlkInd};
   else
     M = widget.LU;
     LU = M(fRow(GLOBAL):lRow(GLOBAL),1:BlkSize);
     L = tril(LU,-1) + speye(BlkSize);
     U = triu(LU);
   end

   rhs = Vec(fRow(VecType):lRow(VecType),:);
   pp = Perm(fRow(GLOBAL):lRow(GLOBAL));
   OutVecPart = U \ ( L \ rhs(pp,:));
   if ~isempty(widget.Qmat)
     qq = Qmat(fRow(GLOBAL):lRow(GLOBAL));
     OutVecPart = OutVecPart(qq,:);
   end
end


% if ConstBlkSize == 1
function [OutVec] = DinvNoBlkApply(A, Vec, Subset, VecType)
    GLOBAL = 1; % defined by mue_include. Inclusion of mue_include is slow
    if ~varexist('Subset'), Subset = []; end
    if ~varexist('VecType'), VecType = GLOBAL; end

    BlkList = SubsetBlkList(A, Subset);
    if (VecType == GLOBAL), VecInds = BlkList; else VecInds = 1:length(BlkList); end

    % input (LU)
     A.FactorBlkDiag(); % This do nothing if the factorization is already computed
     widget = A.GetDiagonalInvData();
     LU = widget.LU; Perm = widget.Perm;

    % output (OutVec)
    OutVec = spalloc(length(BlkList)*1,size(Vec,2),nnz(Vec(VecInds)));

    % computation
    D = spdiags(sparse(LU(BlkList)),0,length(BlkList),length(BlkList));
    % No permutation needed because we have a series of 1x1 block matrices.
    OutVec(1:length(BlkList),:) = D \ Vec(VecInds,:);
end
