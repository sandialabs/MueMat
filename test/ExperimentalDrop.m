% Drop entries from a matrix.
%
% This dropping function can be used as a PreDropFunc or
% PostDropFunc in CoalesceDropFactory.
%
% Foreach a_ij, we keep a_ij if 
%
%  (a_ij^2 / (a_ii*a_jj)) > tol^2 or 
%  (a_ji^2 / (a_ii*a_jj)) > tol^2
%
% If, however, we dropped all the off-diagonal nonzeros in
% a row, we put back one nonzero in the row corresponding 
% to the largest off-diagonal.
%
% Finally, we modify the matrix diagonal so that rowsums
% are preserved after dropping.
%
% See also: SAdropExample
%
function [matrix] = ExperimentalDrop(matrix, level, tol)
  n = size(matrix,2);
  
  invDiag = 1./diag(matrix);
  invDiag = sparse(1:n,1:n,invDiag);
  
  pattern = matrix;
  pattern = pattern.^2;          % square entries
  pattern = pattern*invDiag;     % col scaling
  pattern = (pattern'*invDiag)'; % row scaling

  % At this point pattern_ij should be equal to a_ij^2/(a_ii*a_jj)
  % Now compute the largest off-diag in each row.

  [dummy,LargestOffDiag] = max( pattern - speye(n,n));
  
  toltol  = tol*tol;
  pattern = (pattern>toltol);

  % For any pattern row which contains only 1 nonzero (the diagonal)
  % when the original matrix row contained more than 1 nonzero,
  % add one more nonzero to the pattern corresponding to the 
  % largest offdiagonal entry.

  nnzs     = pattern*ones(n,1);
  orignnzs = (matrix~=0)*ones(n,1);
  nnzs(orignnzs == 1) = 0;
  AllOffDiagsDropped = find(nnzs == 1);

  for k=1:length(AllOffDiagsDropped)
    pattern(AllOffDiagsDropped(k),LargestOffDiag(AllOffDiagsDropped(k))) = 1;
  end;

  pattern = spones(pattern+pattern');
  
  % drop entries from matrix according to pattern and change diagonal
  % so that row sums are the same as before dropping.

  oldsums = matrix*ones(n,1);
  matrix = matrix .* pattern;
  newsums = matrix*ones(n,1);
  matrix = matrix + sparse(1:n,1:n,oldsums-newsums);
