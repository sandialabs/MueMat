%
% Factor the block diagonal matrix BlkDiag (if not already been done).
% Store the result back with the diagonal matrix in the field
% MatrixInvData.
%
% It would definitely be nicer if we have a factor function so that
% we could also do ilu or anything we want here.
%
function FactorBlkDiag(BlkDiag)

   % factorization already computed. To recompute, must
   % delete the old factors before invoking this function.
   %
   if ~ isempty(BlkDiag.GetDiagonalInvData()), return; end

   RowMap = BlkDiag.GetRowMap();

   MatrixData= BlkDiag.GetMatrixData();
   NBlks     = RowMap.NNodes();

   if RowMap.HasConstBlkSize()
      ConstBlkSize = RowMap.ConstBlkSize();
      if (NBlks*ConstBlkSize ~= RowMap.NDOFs())
          fprintf('FactorBlkDiag:: Size mismatch\n');
          keyboard;
      end
      Ntotal = RowMap.NDOFs();

      DenseMat = zeros(Ntotal,ConstBlkSize);
      Perm     = zeros(Ntotal,1);
      if ConstBlkSize == 1,  DenseMat = MatrixData;
      else
         first = 1;
         for i=1:NBlks
            last  = i*ConstBlkSize;
            tt = MatrixData(first:last,1:ConstBlkSize);
            [L,U,P] = lu(tt);
            DenseMat(first:last,1:ConstBlkSize) = U + L - speye(ConstBlkSize,ConstBlkSize);
            [aaa,bbb,ccc] = find(P); [bbb,iii] = sort(bbb); aaa = aaa(iii);
            Perm(aaa+first-1) = (1:ConstBlkSize);
            first = last+1;
         end;
      end
      widget.LU = DenseMat; widget.Perm = Perm; widget.Qmat = [];
      BlkDiag.SetDiagonalInvData(widget);
      BlkDiag.SetApplyInverse(@DinvBlkApply);

   elseif RowMap.HasVariableBlkSize()
      fprintf('factoring block diagonal\n');
      tic;
      VarBlkPtr = RowMap.Vptr();
      Ntotal = VarBlkPtr(NBlks+1)-1;
      Bsizes = zeros(NBlks,1);
      for i=1:NBlks
         Bsizes(i) = VarBlkPtr(i+1)-VarBlkPtr(i);
      end
      MaxBlkSize = max(Bsizes);
      Perm     = zeros(Ntotal,1);
      Qmat     = zeros(Ntotal,1);
      first = VarBlkPtr(1);
      % The main idea here is to iterate over the diagonal blocks:
      %   1) calculate the inverse of each block
      %   2) save L+U-I in row/col/val vectors
      % Once all inverses are done, convert global matrix to sparse form
      for i=1:NBlks
         last  = VarBlkPtr(i+1)-1;
         tt = MatrixData(first:last,1:Bsizes(i));
         if last > first
           % In this case, P*tt*Q = LU.  The Q reduces tt's bandwidth prior to factorization.
           [L,U,P,Q] = lu(tt);
           vv = rand(size(tt,1),1);
         else
           [L,U,P] = lu(tt);
           Q = [];
         end %if last > first
         [I,J,V] = find(U);
         eval(['UIJV' num2str(i) '= [I J V];']);
         [I,J,V] = find(L);
         eval(['LIJV' num2str(i) '= [I J V];']);
         [aaa,bbb,ccc] = find(P); [bbb,iii] = sort(bbb); aaa = aaa(iii);
         Perm(aaa+first-1) = (1:last-first+1);
         [aaa,bbb,ccc] = find(Q); [bbb,iii] = sort(bbb); aaa = aaa(iii);
         Qmat(aaa+first-1) = (1:last-first+1);
         first = last+1;
      end %for i=1:NBlks
      offset = 1;
      for i=1:NBlks
        eval(['tt=UIJV' num2str(i) ';']);
        if length(tt) > 0
          widget.U(offset) = {spconvert(tt)};
        else
          widget.U(offset) = {[]};
        end
        offset = offset+1;
        eval(['clear UIJV' num2str(i) ';']);
      end

      offset = 1;
      for i=1:NBlks
        eval(['tt=LIJV' num2str(i) ';']);
        if length(tt) > 0
          widget.L(offset) = {spconvert(tt)};
        else
          widget.L(offset) = {[]};
        end
        %disp(widget.L(offset))
        offset = offset + 1;
        eval(['clear LIJV' num2str(i) ';']);
      end
      widget.Perm = Perm; widget.Qmat = Qmat;
      BlkDiag.SetDiagonalInvData(widget);
      BlkDiag.SetApplyInverse(@DinvBlkApply);
      toc
   else
      fprintf('FactorBlkDiag:: Unrecognized matrix type\n');
      keyboard;
   end
