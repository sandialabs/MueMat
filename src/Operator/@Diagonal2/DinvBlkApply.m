% This function is a matrix-solve with a block diagonal matrix

function [OutVec] = DinvBlkApply(A, Vec, Subset, VecType)
   mue_include

   if varexist('Subset')  && ~isempty(Subset),  keyboard; end
   if varexist('VecType') && VecType ~= GLOBAL, keyboard; end

   if ~varexist('Subset'), Subset = []; end
   if ~varexist('VecType'), VecType = GLOBAL; end

   Amat = A.GetMatrixData();

   % without using precomputed LU: OutVec = Amat \ Vec;

   % using precomputed LU:
   A.FactorBlkDiag(); % This do nothing if the factorization is already computed
   widget = A.GetDiagonalInvData();

   OutVec = widget.Q_*(widget.U_\(widget.L_\(widget.P_*Vec)));
end