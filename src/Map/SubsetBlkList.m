% Get the block list of a Subset
% used by Diag/DBlkApply.m and Diag/DinvBlkApply.m

function [BlkList] = SubsetBlkList(A, Subset)
   if (~varexist('Subset') || isempty(Subset))
      BlkList = 1:A.GetRowMap().NNodes();
   elseif strcmp(Subset.type,'Contiguous')
      BlkList = Subset.FirstBlk:Subset.LastBlk;
   else
      BlkList = Subset.BlkList;
   end
end