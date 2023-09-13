%
% Take a  list of block rows and convert them to individual
% point rows using the blksize information that is contained
% in operator.
%
function [DOFList] = Node2DOF(BlkRows,Map)


  NBlks =  length(BlkRows);
  if Map.HasConstBlkSize()
     BlkSize = Map.ConstBlkSize();
     DOFList = zeros(NBlks*BlkSize,1);
     for i=1:NBlks
        DOFList((i-1)*BlkSize+1:i*BlkSize) = ...
               ((BlkRows(i)-1)*BlkSize+1:BlkRows(i)*BlkSize);
     end
  else
     VarBlkPtr = Map.Vptr();
     DOFList = zeros(NBlks*Map.MaxBlkSize(),1);
     count = 0;
     for i=1:NBlks,
        Bsize = VarBlkPtr(BlkRows(i)+1)-VarBlkPtr(BlkRows(i));
        DOFList(count+1:count+Bsize) = ...
                  (VarBlkPtr(BlkRows(i)):VarBlkPtr(BlkRows(i)+1)-1);
        count = count+Bsize;
     end
     DOFList = DOFList(1:count);
  end

