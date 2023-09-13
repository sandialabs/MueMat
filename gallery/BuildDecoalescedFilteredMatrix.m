function Matrix2=BuildDecoalescedFilteredMatrix(FullMatrix,CFMatrix)
% BuildDecoalescedFilteredMatrix
%
% SYNTAX M=BuildDecoalescedFilteredMatrix(FullMatrix,CFMatrix)
%
%   FullMatrix  - the full matrix we want to filter
%   CFMatrix    - the coalesced and filtered matrix
%   Rb      - number of rows per block
%   Cb      - number of columns per block
%
% Builds a version of FullMatrix, which has the sparsity pattern
% of a de-coalesced version of CFMatrix, but has the same rowsums
% as FullMatrix.  This is useful in smoothed aggregation for block
% PDEs, so you can use the filtered A for the AP.
%
[M,N]=size(FullMatrix.GetMatrixData());
[M2,N2]=size(CFMatrix.GetMatrixData());
Nr=M/M2; Nc=N/N2;
[R,C,V]=find(CFMatrix.GetMatrixData());
Blk=Nr*Nc;
NNZ=length(R);
Rblk=reshape(repmat(1:Nr,Nc,1),Blk,1);
Cblk=reshape(repmat(1:Nc,Nr,1)',Blk,1);

NewR=reshape(repmat((R-1)*Nr,1,Blk)',NNZ*Blk,1)+reshape(repmat(Rblk',NNZ,1)',NNZ*Blk,1);
NewC=reshape(repmat((C-1)*Nc,1,Blk)',NNZ*Blk,1)+reshape(repmat(Cblk',NNZ,1)',NNZ*Blk,1);
V=ones(NNZ*Blk,1);

% Build the matrix of all 1's
Matrix2=sparse(NewR,NewC,V);

% Get all of the non-dropped entries.
Matrix2=Matrix2.*FullMatrix.GetMatrixData();

% Now probe to get the dropped contributions to the rowsum and
% add those to the diagonal
one_v=ones(N,1);
D=FullMatrix.GetMatrixData()*one_v - Matrix2*one_v;
Matrix2=Matrix2 + spdiags(D,0,N,N);

Matrix2 = Operator(Matrix2,FullMatrix.GetRowMap(), FullMatrix.GetColMap(),@MatlabApply,'Decoalesced Filtered Matrix');
