function [nullspace] = BuildNullSpace(Amat, outputlevel)
% build nullspace by using identity matrices
% warning:  Does not include rigid body rotations as needed for elasticity

   if nargin < 2, outputlevel = 1; end

   RowMap = Amat.GetRowMap();
   nullspace = sparse(RowMap.NDOFs(),RowMap.MaxBlkSize());
   if outputlevel > 0
      fprintf('Nullspace: Using default nullspace size=%d\n',RowMap.MaxBlkSize());
   end
   if(RowMap.MaxBlkSize()==1)
     nullspace(:,1)=1;
   else
     for i=1:RowMap.NNodes()
       DOFs = Node2DOF(i,RowMap);
       bsize = length(DOFs);
       nullspace(DOFs,1:bsize) = speye(bsize,bsize);
     end
   end


