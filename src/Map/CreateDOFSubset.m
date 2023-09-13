%  Subset indicates which rows should be included within an operation such
%  as a matrix-vector product or a block diagonal solve.
%
%      Subset.type     : Contiguous or scatted set of rows ['Contiguous' 'Scattered']
%      Subset.First,   : If Subset.type=='Contiguous', rows correspond to
%      Subset.Last     :    Subset.First:Subset.Last
%      Subset.BlkList  : If Subset.type=='Scattered', the ith block row is BlkList[i]
%      Subset.FirstBlk : If Subset.type=='Contiguous', block rows correspond to
%      Subset.LastBlk  :   Subset.FirstBlk:Subset.LastBlk

function [Subset] = CreateDOFSubset(Map, AllOrSubset, FirstBlk, LastBlk, BlkList)

   if strcmp(AllOrSubset,'All'),
      FirstBlk = 1;
      LastBlk  = Map.NNodes();
      AllOrSubset = 'Contiguous';
   end

   if strcmp(AllOrSubset,'Contiguous'),
      Subset.type = 'Contiguous';
      Subset.FirstBlk = FirstBlk;
      Subset.LastBlk = LastBlk;
      if Map.HasConstBlkSize()
         constBlkSize = Map.ConstBlkSize();
         Subset.First    = (FirstBlk-1)*constBlkSize+1;
         Subset.Last     = LastBlk*constBlkSize;
      else
         VarBlkPtr = Map.Vptr();
         Subset.First = VarBlkPtr(FirstBlk);
         Subset.Last  = VarBlkPtr(LastBlk+1)-1;
      end
   elseif strcmp(AllOrSubset,'Scattered'),
      Subset.type = 'Scattered';
      Subset.BlkList = BlkList;
      Subset.First = -1;
      Subset.Last  = -1;
   else
      fprintf('CreateDOFSubset:: Unrecognized AllOrSubset = %s\n',AllOrSubset);
      keyboard;
   end
