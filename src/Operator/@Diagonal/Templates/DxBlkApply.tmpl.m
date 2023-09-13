#ifndef INV

#define FUNC D
#define FUNCCOMMENT matrix-vector product
#define AMAT Amat                // input name

#else /* ifndef INV */

#define FUNC Dinv
#define FUNCCOMMENT matrix-solve
#define AMAT LU                  // input name

#endif /* ifndef INV */

#define FUNCNAME2(func,type) func ## type ## BlkApply
#define FUNCNAME(func,type)  FUNCNAME2(func,type)

%
% This file is generated from the template DxBlkApply.tmpl.m
%
% This function is a FUNCCOMMENT with a block diagonal matrix where
% the blocks have a constant size (HasConstBlkSize) or have
% different sizes (HasVariableBlkSize).
%
% There is also a possibility of performing the FUNCCOMMENT only on a
% subset of block rows which can be contiguous or scattered
% (see CreateDOFSubset.m).
#ifdef INV
% The block diagonal matrix is either already factored or factors will
% be produced (via FactorBlkDiag.m).
#endif
% The block diagonal matrix AMAT is stored in a special way (i.e. different
% from a general matrix operator) in order to gain some efficiency (see
% ExtractBlkDiag.m and ExtractNonContigBlkDiag.m for details on the
% format).
%
% When a FUNCCOMMENT is performed only involving a subset
% of rows, there are two possibilities for the form of the input vector:
%     Case 1:  Standard case where the input vector is the entire
%              input vector for the matrix even though only a subset
%              of this data is used in performing the subset
%              FUNCCOMMENT.
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
%  For (ConstBlkSize == 1), the simplier subfunction FUNCNAME(FUNC,No) is called.
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

function [OutVec] = FUNCNAME(FUNC,)(A, Vec, Subset, VecType)
   GLOBAL = 1; % defined by mue_include. Inclusion of mue_include is slow
   if ~varexist('Subset'), Subset = []; end
   if ~varexist('VecType'), VecType = GLOBAL; end

   RowMap = A.GetRowMap();

   if RowMap.HasConstBlkSize() && (RowMap.ConstBlkSize() == 1)
     OutVec = FUNCNAME(FUNC,No)(A, Vec, Subset, VecType); %optim
   else

     BlkList = SubsetBlkList(A, Subset);

     % input (AMAT)
     #include "input.tmpl.m" /* setup Amat or LU */

     % output (OutVec)
     nRow = SubsetNRow(A, Subset);
     %TEST: OutVecTEST = zeros(nRow,size(Vec,2));
     OutVecIJV = IJVCellMatrix(length(BlkList));  % Unassembled matrix

     % computation
     if RowMap.HasConstBlkSize()
       BlkSize = RowMap.ConstBlkSize();

       fLocalRow = 1;
       for i=BlkList;
         fGlobalRow = (i-1)*BlkSize+1;
         lGlobalRow = fGlobalRow + BlkSize-1;
         lLocalRow  = fLocalRow + BlkSize-1;

         #include "main_loop.tmpl.m" /* fill OutVecIJV */

         fLocalRow = lLocalRow+1;
       end

     else
       VarBlkPtr = RowMap.Vptr();

       fLocalRow = 1;
       for i=BlkList;
         fGlobalRow = VarBlkPtr(i);
         lGlobalRow = VarBlkPtr(i+1)-1;
         BlkSize    = lGlobalRow-fGlobalRow+1;
         lLocalRow  = fLocalRow + BlkSize-1;

         #include "main_loop.tmpl.m" /* fill OutVecIJV */

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
function [OutVecPart] = FUNCNAME(FUNC,SubFunc)(M, Vec, fRow, lRow, BlkSize, VecType)
   GLOBAL = 1; % defined by mue_include. Inclusion of mue_include is slow

#ifdef INV
   LU = M(fRow(GLOBAL):lRow(GLOBAL),1:BlkSize);
   L  = tril(LU,-1) + speye(BlkSize);
   U  = triu(LU);

   OutVecPart = U \ ( L \ Vec(fRow(VecType):lRow(VecType),:));
#else
   OutVecPart = M(fRow(GLOBAL):lRow(GLOBAL)) * Vec(fRow(VecType):lRow(VecType),:);
#endif

end


% if ConstBlkSize == 1
function [OutVec] = FUNCNAME(FUNC,No)(A, Vec, Subset, VecType)
    GLOBAL = 1; % defined by mue_include. Inclusion of mue_include is slow
    if ~varexist('Subset'), Subset = []; end
    if ~varexist('VecType'), VecType = GLOBAL; end

    BlkList = SubsetBlkList(A, Subset);
    if (VecType == GLOBAL), VecInds = BlkList; else VecInds = 1:length(BlkList); end

    % input (AMAT)
    #include "input.tmpl.m" /* setup Amat or LU */

    % output (OutVec)
    OutVec = spalloc(length(BlkList)*1,size(Vec,2),nnz(Vec(VecInds)));

    % computation
    D = diag(sparse(AMAT(BlkList)));
#ifdef INV
    OutVec(1:length(BlkList),:) = D \ Vec(VecInds,:);
#else
    OutVec(1:length(BlkList),:) = D * Vec(VecInds,:);
#endif
end
