% This function is a matrix-vector product with a block diagonal matrix

function [OutVec] = DBlkApply(A, Vec, Subset, VecType)
   mue_include

   if varexist('Subset')  && ~isempty(Subset),  keyboard; end
   if varexist('VecType') && VecType ~= GLOBAL, keyboard; end

   if ~varexist('Subset'),  Subset = []; end
   if ~varexist('VecType'), VecType = GLOBAL; end

    OutVec = A*Vec;

end
