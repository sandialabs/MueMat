% Compute the number of rows of a Subset
% used by Diag/DBlkApply.m and Diag/DinvBlkApply.m

function [nRow] = SubsetNRow(A, Subset)
% TODO: Duplicated on DinvCblkApply and DCblkApply
   %if ~varexist('Subset'), Subset = []; end

   BlkList   = SubsetBlkList(A, Subset);

   if A.GetRowMap().HasConstBlkSize()
     BlkSize = A.GetRowMap().ConstBlkSize();
     nRow = length(BlkList)*BlkSize;
   else
     VarBlkPtr = A.GetRowMap().Vptr();

     if (~isempty(Subset) && strcmp(Subset.type,'Scattered'))
       nRow = 0;
       for i=BlkList
         nRow = nRow + VarBlkPtr(i+1) - VarBlkPtr(i);
       end
     else
       nRow = VarBlkPtr(BlkList(end)+1)-VarBlkPtr(BlkList(1));
     end

   end
 end

