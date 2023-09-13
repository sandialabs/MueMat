function FactorBlkDiag(BlkDiag)

% factorization already computed. To recompute, must
   % delete the old factors before invoking this function.
   %
   if ~isempty(BlkDiag.GetDiagonalInvData()), return; end

   MatrixData = BlkDiag.GetMatrixData();

   [widget.L_,widget.U_,widget.P_,widget.Q_] = lu(MatrixData);

   BlkDiag.SetDiagonalInvData(widget);
   BlkDiag.SetApplyInverse(@DinvBlkApply);

end