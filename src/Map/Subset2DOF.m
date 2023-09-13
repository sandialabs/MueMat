%
%  Given a Subset object, convert this to a list of point rows.
%
function [indices] = Subset2DOF(Subset,Map)

   if  isempty(Subset),
      indices=(1:Map.NDOFs())';
   elseif strcmp(Subset.type,'Contiguous'),
         indices = (Subset.First:Subset.Last)';
   elseif strcmp(Subset.type,'Scattered'),
         indices = Node2DOF(Subset.BlkList, Map);
   else
         fprintf('Unknown Subset type:: %s\n',Subset.type);
         keyboard;
   end



