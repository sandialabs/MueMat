% Build the block matrix of a 1D Laplace problem
% with either Constant or Variable sized blocks.
%
% The matrix represents a finite difference discretization of the
% poisson problem on a uniform grid and with Dirichlet boundary
% conditions.
%
% In the variable block size case, the kth DOF within a node has a
% nonzero connection to the kth DOF of any adjacent nodes. If
% adjacent nodes don't have a kth DOF, then this DOF is obviously
% not linked to this adjacent node.
%
% Note: picking epsi ~= 0 puts some small crud into the matrix blocks
%       so that they are not all diagonal.
%
% See also: BuildLaplace1D, TutorialOperator

function [Amat] = BuildLaplace1DBlk(ConstBlkSize, VarBlkPtr, NBlks)

   epsi = 1.e-4;
   if (ConstBlkSize > 0)
      Perturb = epsi*ones(ConstBlkSize,ConstBlkSize);
      fprintf('\n*****************************************************************************\n\n');
      fprintf('Building %d x %d block matrix with constant block size = %d\n',NBlks,NBlks,ConstBlkSize);
      Ntotal = NBlks*ConstBlkSize;
      newMap = Map(NBlks, ConstBlkSize);
      MatrixData = spalloc(Ntotal,Ntotal,3*ConstBlkSize*Ntotal);
      for i=1:NBlks
         FirstRow = (i-1)*ConstBlkSize+1;
         if epsi ~= 0,
            ii= FirstRow; kk=i*ConstBlkSize;
            MatrixData(ii:kk,ii:kk) = Perturb;
            if i ~= NBlks, MatrixData(ii:kk,ii+ConstBlkSize:kk+ConstBlkSize)=Perturb; end
            if i ~= 1,     MatrixData(ii:kk,ii-ConstBlkSize:kk-ConstBlkSize)=Perturb; end
         end
         for j=1:ConstBlkSize
            MatrixData(FirstRow+j-1,FirstRow+j-1) = 2;
            if i ~= 1,     MatrixData(FirstRow+j-1,FirstRow+j-1-ConstBlkSize)= -1; end
            if i ~= NBlks, MatrixData(FirstRow+j-1,FirstRow+j-1+ConstBlkSize)= -1; end
         end
      end
   else
      Ntotal = VarBlkPtr(NBlks+1)-1;
      if (Ntotal < NBlks)
         fprintf('BuildPois:: VarBlkPtr not set properly\n');
         keyboard;
      end
      Bsizes = zeros(NBlks,1);
      for i=1:NBlks
         Bsizes(i) = VarBlkPtr(i+1)-VarBlkPtr(i);
      end
      newMap = Map(NBlks, Bsizes);
      fprintf('\n*****************************************************************************\n\n');
      fprintf('Building %d x %d block matrix with variable block size = [%d,%d]\n',NBlks,NBlks,min(Bsizes),max(Bsizes));
      MatrixData = spalloc(Ntotal,Ntotal,3*max(Bsizes)*Ntotal);
      for i=1:NBlks
         if (epsi ~= 0)
            ii= VarBlkPtr(i); kk = VarBlkPtr(i+1)-1; bb=Bsizes(i);
            MatrixData(ii:kk,ii:kk) = epsi*ones(bb,bb);
            if i ~=NBlks,cc=Bsizes(i+1); MatrixData(ii:kk,kk+1:kk+cc)=epsi*ones(bb,cc);end
            if i ~=1,    cc=Bsizes(i-1); MatrixData(ii:kk,ii-cc:ii-1)=epsi*ones(bb,cc);end
         end
         for j=1:Bsizes(i)
            MatrixData(VarBlkPtr(i)+j-1,VarBlkPtr(i)+j-1) = 2;
            if (i~= 1) && (Bsizes(i-1) >= j),
                MatrixData(VarBlkPtr(i)+j-1,VarBlkPtr(i-1)+j-1) = -1;
            end
            if (i~= NBlks) && (Bsizes(i+1) >= j),
               MatrixData(VarBlkPtr(i)+j-1,VarBlkPtr(i+1)+j-1) = -1;
            end
         end
      end
   end
   %Amat.RowMap  = newMap;
   %Amat.ColMap  = newMap;
   %Amat.Apply      = @MatlabApply;
   %Amat.MatrixData = MatrixData;
   Amat = Operator(MatrixData,newMap,newMap,@MatlabApply, 'Laplace1DBlk');
